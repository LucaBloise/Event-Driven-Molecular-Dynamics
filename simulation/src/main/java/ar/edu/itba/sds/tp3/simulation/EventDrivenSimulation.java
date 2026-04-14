package ar.edu.itba.sds.tp3.simulation;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;

public final class EventDrivenSimulation {
    private final SimulationConfig config;
    private final Random random;
    private final List<Particle> particles;
    private final PriorityQueue<Event> priorityQueue;

    private double currentTime;
    private long processedEvents;
    private long writtenFrames;

    public EventDrivenSimulation(final SimulationConfig config) {
        this.config = config;
        this.random = new Random(config.getSeed());
        this.particles = new ArrayList<>(config.getParticleCount());
        this.priorityQueue = new PriorityQueue<>();
        this.currentTime = 0.0;
        this.processedEvents = 0L;
        this.writtenFrames = 0L;
    }

    public SimulationResult run(final Path outputFile) throws IOException {
        initializeParticles();
        scheduleInitialEvents();

        final long startNanos = System.nanoTime();
        double lastWrittenTime = Double.NEGATIVE_INFINITY;

        try (SimulationOutputWriter writer = new SimulationOutputWriter(outputFile)) {
            writer.writeHeader();
            writer.writeFrame(writtenFrames++, processedEvents, currentTime, "INITIAL", -1, -1, particles);
            lastWrittenTime = currentTime;

            while (currentTime < config.getEndTime()) {
                final Event next = pollNextValidEvent();
                if (next == null) {
                    advanceAll(config.getEndTime() - currentTime);
                    currentTime = config.getEndTime();
                    break;
                }

                if (next.getTime() > config.getEndTime()) {
                    advanceAll(config.getEndTime() - currentTime);
                    currentTime = config.getEndTime();
                    break;
                }

                advanceAll(next.getTime() - currentTime);
                currentTime = next.getTime();

                final ProcessedEvent processedEvent = processEvent(next);
                processedEvents++;

                if (processedEvents % config.getSnapshotEveryEvents() == 0) {
                    writer.writeFrame(
                            writtenFrames++,
                            processedEvents,
                            currentTime,
                            processedEvent.label,
                            processedEvent.particleA,
                            processedEvent.particleB,
                            particles
                    );
                    lastWrittenTime = currentTime;
                }

                scheduleAfterEvent(processedEvent);
            }

            if (Math.abs(lastWrittenTime - currentTime) > Particle.EPS) {
                writer.writeFrame(writtenFrames++, processedEvents, currentTime, "FINAL", -1, -1, particles);
            }
        }

        final long endNanos = System.nanoTime();
        final double executionSeconds = (endNanos - startNanos) / 1.0e9;
        return new SimulationResult(currentTime, processedEvents, writtenFrames, executionSeconds);
    }

    private void initializeParticles() {
        particles.clear();

        final double innerRadius = config.getInnerCollisionRadius();
        final double outerRadius = config.getOuterCollisionRadius();
        final double minCenterDistance = 2.0 * config.getParticleRadius();
        final double innerSq = innerRadius * innerRadius;
        final double outerSq = outerRadius * outerRadius;

        for (int id = 0; id < config.getParticleCount(); id++) {
            boolean placed = false;

            for (int attempt = 0; attempt < config.getMaxPlacementAttemptsPerParticle(); attempt++) {
                final double theta = random.nextDouble() * 2.0 * Math.PI;
                final double radial = Math.sqrt(random.nextDouble() * (outerSq - innerSq) + innerSq);

                final double x = radial * Math.cos(theta);
                final double y = radial * Math.sin(theta);

                if (overlapsExisting(x, y, minCenterDistance)) {
                    continue;
                }

                final double direction = random.nextDouble() * 2.0 * Math.PI;
                final double vx = config.getParticleSpeed() * Math.cos(direction);
                final double vy = config.getParticleSpeed() * Math.sin(direction);

                particles.add(new Particle(
                        id,
                        x,
                        y,
                        vx,
                        vy,
                        config.getParticleMass(),
                        config.getParticleRadius(),
                        ParticleState.FRESH
                ));
                placed = true;
                break;
            }

            if (!placed) {
                throw new IllegalStateException(
                        "Could not place all particles without overlap. Try lowering --n or increasing --max-placement-attempts"
                );
            }
        }
    }

    private boolean overlapsExisting(final double x, final double y, final double minCenterDistance) {
        final double minDistSquared = minCenterDistance * minCenterDistance;
        for (final Particle particle : particles) {
            final double dx = x - particle.getX();
            final double dy = y - particle.getY();
            final double distSq = dx * dx + dy * dy;
            if (distSq < minDistSquared) {
                return true;
            }
        }
        return false;
    }

    private void scheduleInitialEvents() {
        priorityQueue.clear();
        for (int i = 0; i < particles.size(); i++) {
            scheduleParticleEvents(i);
        }
    }

    private void scheduleAfterEvent(final ProcessedEvent event) {
        scheduleParticleEvents(event.particleA);
        if (event.particleB >= 0) {
            scheduleParticleEvents(event.particleB);
        }
    }

    private void scheduleParticleEvents(final int particleIndex) {
        if (particleIndex < 0 || particleIndex >= particles.size()) {
            return;
        }

        final Particle particle = particles.get(particleIndex);

        final double dtOuter = particle.timeToHitCircle(config.getOuterCollisionRadius());
        scheduleWallEvent(dtOuter, EventType.OUTER_WALL_COLLISION, particleIndex);

        final double dtInner = particle.timeToHitCircle(config.getInnerCollisionRadius());
        scheduleWallEvent(dtInner, EventType.INNER_WALL_COLLISION, particleIndex);

        for (int otherIndex = 0; otherIndex < particles.size(); otherIndex++) {
            if (otherIndex == particleIndex) {
                continue;
            }

            final Particle other = particles.get(otherIndex);
            final double dt = particle.timeToHitParticle(other);
            if (!isFiniteFutureTime(dt)) {
                continue;
            }

            final double eventTime = currentTime + dt;
            if (eventTime > config.getEndTime() + Particle.EPS) {
                continue;
            }

            priorityQueue.add(new Event(
                    eventTime,
                    EventType.PARTICLE_COLLISION,
                    particleIndex,
                    otherIndex,
                    particle.getCollisionCount(),
                    other.getCollisionCount()
            ));
        }
    }

    private void scheduleWallEvent(final double dt, final EventType type, final int particleIndex) {
        if (!isFiniteFutureTime(dt)) {
            return;
        }

        final double eventTime = currentTime + dt;
        if (eventTime > config.getEndTime() + Particle.EPS) {
            return;
        }

        final Particle particle = particles.get(particleIndex);
        priorityQueue.add(new Event(
                eventTime,
                type,
                particleIndex,
                -1,
                particle.getCollisionCount(),
                -1
        ));
    }

    private Event pollNextValidEvent() {
        while (!priorityQueue.isEmpty()) {
            final Event event = priorityQueue.poll();
            if (event.isValid(particles)) {
                return event;
            }
        }
        return null;
    }

    private ProcessedEvent processEvent(final Event event) {
        switch (event.getType()) {
            case PARTICLE_COLLISION:
                return processParticleCollision(event.getAIndex(), event.getBIndex());
            case OUTER_WALL_COLLISION:
                return processOuterWallCollision(event.getAIndex());
            case INNER_WALL_COLLISION:
                return processInnerWallCollision(event.getAIndex());
            default:
                throw new IllegalStateException("Unknown event type: " + event.getType());
        }
    }

    private ProcessedEvent processParticleCollision(final int aIndex, final int bIndex) {
        final Particle a = particles.get(aIndex);
        final Particle b = particles.get(bIndex);
        a.bounceOffParticle(b);
        return new ProcessedEvent("PARTICLE_COLLISION", aIndex, bIndex);
    }

    private ProcessedEvent processOuterWallCollision(final int index) {
        final Particle particle = particles.get(index);
        particle.snapToCircle(config.getOuterCollisionRadius());
        particle.bounceOffCircularBoundary();

        if (particle.getState() == ParticleState.USED) {
            particle.setState(ParticleState.FRESH);
        }

        return new ProcessedEvent("OUTER_WALL_COLLISION", index, -1);
    }

    private ProcessedEvent processInnerWallCollision(final int index) {
        final Particle particle = particles.get(index);
        particle.snapToCircle(config.getInnerCollisionRadius());
        particle.bounceOffCircularBoundary();

        if (particle.getState() == ParticleState.FRESH) {
            particle.setState(ParticleState.USED);
        }

        return new ProcessedEvent("INNER_WALL_COLLISION", index, -1);
    }

    private void advanceAll(double dt) {
        if (dt < -Particle.EPS) {
            throw new IllegalStateException("Negative time step detected: " + dt);
        }

        if (dt < 0.0) {
            dt = 0.0;
        }

        for (final Particle particle : particles) {
            particle.advance(dt);
        }
    }

    private static boolean isFiniteFutureTime(final double dt) {
        return Double.isFinite(dt) && dt > Particle.EPS;
    }

    private static final class ProcessedEvent {
        private final String label;
        private final int particleA;
        private final int particleB;

        private ProcessedEvent(final String label, final int particleA, final int particleB) {
            this.label = label;
            this.particleA = particleA;
            this.particleB = particleB;
        }
    }
}