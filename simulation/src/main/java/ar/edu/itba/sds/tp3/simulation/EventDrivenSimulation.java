package ar.edu.itba.sds.tp3.simulation;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Random;

public final class EventDrivenSimulation {
    private final SimulationConfig config;
    private final Random random;
    private final List<Particle> particles;
    private final PriorityQueue<Event> priorityQueue;
    private final Map<Long, Event> activePairCollisionEvents;

    private double currentTime;
    private long processedEvents;
    private long writtenFrames;

    public EventDrivenSimulation(final SimulationConfig config) {
        this.config = config;
        this.random = new Random(config.getSeed());
        this.particles = new ArrayList<>(config.getParticleCount());
        this.priorityQueue = new PriorityQueue<>();
        this.activePairCollisionEvents = new HashMap<>();
        this.currentTime = 0.0;
        this.processedEvents = 0L;
        this.writtenFrames = 0L;
    }

    public SimulationResult run(final Path outputFile) throws IOException {
        initializeParticles();
        scheduleInitialEvents();

        final long startNanos = System.nanoTime();
        final boolean writeFrames = outputFile != null;
        final boolean writeDeltaEvents = outputFile != null && config.shouldWriteDeltaOutputEvents();
        double lastWrittenTime = Double.NEGATIVE_INFINITY;

        SimulationOutputWriter writer = null;
        SimulationDeltaOutputWriter deltaWriter = null;
        try {
            if (writeFrames && !writeDeltaEvents) {
                writer = new SimulationOutputWriter(outputFile);
                writer.writeHeader();
                writer.writeFrame(writtenFrames++, processedEvents, currentTime, "INITIAL", -1, -1, particles);
                lastWrittenTime = currentTime;
            }

            if (writeDeltaEvents) {
                deltaWriter = new SimulationDeltaOutputWriter(outputFile);
                deltaWriter.writeHeader();
                deltaWriter.writeInitialState(particles);
            }

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

                if (writeFrames && !writeDeltaEvents && processedEvents % config.getSnapshotEveryEvents() == 0) {
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

                if (writeDeltaEvents) {
                    deltaWriter.writeEvent(
                            processedEvents,
                            currentTime,
                            processedEvent.label,
                            processedEvent.particleA,
                            processedEvent.particleB,
                            processedEvent.changedParticles
                    );
                }

                scheduleAfterEvent(processedEvent);
            }

            if (writeFrames && !writeDeltaEvents && Math.abs(lastWrittenTime - currentTime) > Particle.EPS) {
                writer.writeFrame(writtenFrames++, processedEvents, currentTime, "FINAL", -1, -1, particles);
            }

            if (writeDeltaEvents) {
                deltaWriter.writeFinal(processedEvents, currentTime);
            }
        } finally {
            if (writer != null) {
                writer.close();
            }
            if (deltaWriter != null) {
                deltaWriter.close();
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

        // Strategy: try random placement first; if it fails, use hexagonal lattice.
        boolean useRandomFirst = config.getParticleCount() <= 700;
        boolean placedRandomly = useRandomFirst && placeRandomly(innerRadius, outerRadius, minCenterDistance);
        
        if (!placedRandomly) {
            placeOnLattice(innerRadius, outerRadius, minCenterDistance);
        }

        // Set velocities for all particles
        for (int id = 0; id < particles.size(); id++) {
            final Particle p = particles.get(id);
            final double direction = random.nextDouble() * 2.0 * Math.PI;
            final double vx = config.getParticleSpeed() * Math.cos(direction);
            final double vy = config.getParticleSpeed() * Math.sin(direction);
            p.setVx(vx);
            p.setVy(vy);
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

    private boolean placeRandomly(final double innerRadius, final double outerRadius, final double minCenterDistance) {
        final int n = config.getParticleCount();
        final double innerSq = innerRadius * innerRadius;
        final double outerSq = outerRadius * outerRadius;
        final int maxAttempts = Math.max(20000, 50 * n);

        for (int id = 0; id < n; id++) {
            boolean placed = false;
            for (int attempt = 0; attempt < maxAttempts; attempt++) {
                final double theta = random.nextDouble() * 2.0 * Math.PI;
                final double radial = Math.sqrt(random.nextDouble() * (outerSq - innerSq) + innerSq);
                final double x = radial * Math.cos(theta);
                final double y = radial * Math.sin(theta);

                if (overlapsExisting(x, y, minCenterDistance)) {
                    continue;
                }

                particles.add(new Particle(
                        id,
                        x,
                        y,
                        0.0,
                        0.0,
                        config.getParticleMass(),
                        config.getParticleRadius(),
                        ParticleState.FRESH
                ));
                placed = true;
                break;
            }
            if (!placed) {
                return false;
            }
        }
        return true;
    }

    private void placeOnLattice(final double innerRadius, final double outerRadius, final double minCenterDistance) {
        final double r = config.getParticleRadius();
        final double dx = 2.0 * r;
        final double dy = Math.sqrt(3.0) * r;

        List<Position> bestPositions = buildHexLatticePositions(outerRadius, innerRadius, dx, dy, 0.0, 0.0);

        // Try multiple lattice offsets and keep the densest valid arrangement.
        for (int attempt = 0; attempt < 32; attempt++) {
            final double shiftX = random.nextDouble() * dx;
            final double shiftY = random.nextDouble() * dy;
            final List<Position> candidate = buildHexLatticePositions(outerRadius, innerRadius, dx, dy, shiftX, shiftY);
            if (candidate.size() > bestPositions.size()) {
                bestPositions = candidate;
            }
        }

        if (bestPositions.size() < config.getParticleCount()) {
            throw new IllegalStateException(
                    "Not enough lattice positions for N=" + config.getParticleCount()
                            + " (max=" + bestPositions.size() + "). Increase domain or decrease particle count."
            );
        }

        Collections.shuffle(bestPositions, random);
        for (int id = 0; id < config.getParticleCount(); id++) {
            final Position p = bestPositions.get(id);
            particles.add(new Particle(
                    id,
                    p.x,
                    p.y,
                    0.0,
                    0.0,
                    config.getParticleMass(),
                    config.getParticleRadius(),
                    ParticleState.FRESH
            ));
        }
    }

    private List<Position> buildHexLatticePositions(final double outerRadius,
                                                     final double innerRadius,
                                                     final double dx,
                                                     final double dy,
                                                     final double shiftX,
                                                     final double shiftY) {
        final List<Position> positions = new ArrayList<>();
        int row = 0;

        for (double yy = -outerRadius + shiftY; yy <= outerRadius + 1e-12; yy += dy) {
            final double rowOffset = (row % 2 == 0) ? 0.0 : dx / 2.0;
            for (double xx = -outerRadius + rowOffset + shiftX; xx <= outerRadius + 1e-12; xx += dx) {
                final double dist = Math.sqrt(xx * xx + yy * yy);
                if (dist >= innerRadius && dist <= outerRadius) {
                    positions.add(new Position(xx, yy));
                }
            }
            row++;
        }

        return positions;
    }

    private void scheduleInitialEvents() {
        priorityQueue.clear();
        activePairCollisionEvents.clear();
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

            schedulePairCollisionEvent(new Event(
                    eventTime,
                    EventType.PARTICLE_COLLISION,
                    particleIndex,
                    otherIndex,
                    particle.getCollisionCount(),
                    other.getCollisionCount()
            ));
        }
    }

    private void schedulePairCollisionEvent(final Event candidate) {
        final long key = pairKey(candidate.getAIndex(), candidate.getBIndex());
        final Event current = activePairCollisionEvents.get(key);

        if (current != null && current.isValid(particles) && current.getTime() <= candidate.getTime() + Particle.EPS) {
            return;
        }

        activePairCollisionEvents.put(key, candidate);
        priorityQueue.add(candidate);
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
            if (event.getType() == EventType.PARTICLE_COLLISION) {
                final long key = pairKey(event.getAIndex(), event.getBIndex());
                final Event current = activePairCollisionEvents.get(key);

                if (current != event) {
                    continue;
                }

                if (!event.isValid(particles)) {
                    activePairCollisionEvents.remove(key);
                    continue;
                }

                activePairCollisionEvents.remove(key);
                return event;
            }

            if (event.isValid(particles)) {
                return event;
            }
        }
        return null;
    }

    private static long pairKey(final int indexA, final int indexB) {
        final int min = Math.min(indexA, indexB);
        final int max = Math.max(indexA, indexB);
        return ((long) min << 32) | (max & 0xffffffffL);
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
        return new ProcessedEvent("PARTICLE_COLLISION", aIndex, bIndex, List.of(a, b));
    }

    private ProcessedEvent processOuterWallCollision(final int index) {
        final Particle particle = particles.get(index);
        particle.snapToCircle(config.getOuterCollisionRadius());
        particle.bounceOffCircularBoundary();

        if (particle.getState() == ParticleState.USED) {
            particle.setState(ParticleState.FRESH);
        }

        return new ProcessedEvent("OUTER_WALL_COLLISION", index, -1, List.of(particle));
    }

    private ProcessedEvent processInnerWallCollision(final int index) {
        final Particle particle = particles.get(index);
        particle.snapToCircle(config.getInnerCollisionRadius());
        particle.bounceOffCircularBoundary();

        if (particle.getState() == ParticleState.FRESH) {
            particle.setState(ParticleState.USED);
        }

        return new ProcessedEvent("INNER_WALL_COLLISION", index, -1, List.of(particle));
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
        private final List<Particle> changedParticles;

        private ProcessedEvent(
                final String label,
                final int particleA,
                final int particleB,
                final List<Particle> changedParticles
        ) {
            this.label = label;
            this.particleA = particleA;
            this.particleB = particleB;
            this.changedParticles = changedParticles;
        }
    }

    private static final class Position {
        final double x;
        final double y;

        Position(final double x, final double y) {
            this.x = x;
            this.y = y;
        }
    }
}