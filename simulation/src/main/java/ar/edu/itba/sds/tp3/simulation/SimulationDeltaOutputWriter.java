package ar.edu.itba.sds.tp3.simulation;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.List;
import java.util.Locale;

public final class SimulationDeltaOutputWriter implements Closeable {
    private final BufferedWriter writer;

    public SimulationDeltaOutputWriter(final Path outputFile) throws IOException {
        if (outputFile.getParent() != null) {
            Files.createDirectories(outputFile.getParent());
        }

        this.writer = Files.newBufferedWriter(
                outputFile,
                StandardCharsets.UTF_8,
                StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING,
                StandardOpenOption.WRITE
        );
    }

    public void writeHeader() throws IOException {
        writer.write("# Event-Driven Molecular Dynamics output delta v1");
        writer.newLine();
        writer.write("# INITIAL_PARTICLE <id> <x_m> <y_m> <vx_m_s> <vy_m_s> <state> <color_r> <color_g> <color_b>");
        writer.newLine();
        writer.write("# EVENT <event_index> <time_s> <event_type> <particle_a> <particle_b>");
        writer.newLine();
        writer.write("# CHANGED_PARTICLE <id> <x_m> <y_m> <vx_m_s> <vy_m_s> <state> <color_r> <color_g> <color_b>");
        writer.newLine();
    }

    public void writeInitialState(final List<Particle> particles) throws IOException {
        writer.write("BEGIN_INITIAL_STATE");
        writer.newLine();

        for (final Particle particle : particles) {
            final ParticleState state = particle.getState();
            writer.write(String.format(
                    Locale.US,
                    "INITIAL_PARTICLE %d %.10f %.10f %.10f %.10f %s %d %d %d",
                    particle.getId(),
                    particle.getX(),
                    particle.getY(),
                    particle.getVx(),
                    particle.getVy(),
                    state.name(),
                    state.getRed(),
                    state.getGreen(),
                    state.getBlue()
            ));
            writer.newLine();
        }

        writer.write("END_INITIAL_STATE");
        writer.newLine();
    }

    public void writeEvent(
            final long eventIndex,
            final double time,
            final String eventType,
            final int particleA,
            final int particleB,
            final List<Particle> changedParticles
    ) throws IOException {
        writer.write(String.format(
                Locale.US,
                "EVENT %d %.10f %s %d %d",
                eventIndex,
                time,
                eventType,
                particleA,
                particleB
        ));
        writer.newLine();

        for (final Particle particle : changedParticles) {
            final ParticleState state = particle.getState();
            writer.write(String.format(
                    Locale.US,
                    "CHANGED_PARTICLE %d %.10f %.10f %.10f %.10f %s %d %d %d",
                    particle.getId(),
                    particle.getX(),
                    particle.getY(),
                    particle.getVx(),
                    particle.getVy(),
                    state.name(),
                    state.getRed(),
                    state.getGreen(),
                    state.getBlue()
            ));
            writer.newLine();
        }

        writer.write("END_EVENT");
        writer.newLine();
    }

    public void writeFinal(final long eventIndex, final double time) throws IOException {
        writer.write(String.format(Locale.US, "FINAL %d %.10f", eventIndex, time));
        writer.newLine();
    }

    @Override
    public void close() throws IOException {
        writer.close();
    }
}
