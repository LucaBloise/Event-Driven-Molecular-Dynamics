package ar.edu.itba.sds.tp3.simulation;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

public final class SimulationPropertiesWriter {
    private SimulationPropertiesWriter() {
    }

    public static void write(
            final Path propertiesFile,
            final SimulationConfig config,
            final SimulationResult result,
            final Path outputFile
    ) throws IOException {
        final Map<String, String> values = new LinkedHashMap<>();
        values.put("format_version", "1");
        values.put("output_format", config.shouldWriteDeltaOutputEvents() ? "event-delta-v1" : "none");
        values.put("output_file", outputFile != null ? outputFile.getFileName().toString() : "none");

        values.put("n_particles", Integer.toString(config.getParticleCount()));
        values.put("tf_seconds", format(config.getEndTime()));
        values.put("domain_diameter_m", format(config.getDomainDiameter()));
        values.put("obstacle_radius_m", format(config.getObstacleRadius()));
        values.put("particle_radius_m", format(config.getParticleRadius()));
        values.put("particle_mass_kg", format(config.getParticleMass()));
        values.put("particle_speed_m_s", format(config.getParticleSpeed()));
        values.put("inner_collision_radius_m", format(config.getInnerCollisionRadius()));
        values.put("outer_collision_radius_m", format(config.getOuterCollisionRadius()));
        values.put("seed", Long.toString(config.getSeed()));
        values.put("snapshot_every_events", Integer.toString(config.getSnapshotEveryEvents()));
        values.put("write_output_frames", "false");
        values.put("write_delta_output_events", Boolean.toString(config.shouldWriteDeltaOutputEvents()));
        values.put("max_placement_attempts_per_particle", Integer.toString(config.getMaxPlacementAttemptsPerParticle()));

        values.put("final_simulated_time_s", format(result.getFinalSimulatedTime()));
        values.put("processed_events", Long.toString(result.getProcessedEvents()));
        values.put("written_frames", Long.toString(result.getWrittenFrames()));
        values.put("execution_time_s", format(result.getExecutionSeconds()));

        if (propertiesFile.getParent() != null) {
            Files.createDirectories(propertiesFile.getParent());
        }

        try (BufferedWriter writer = Files.newBufferedWriter(
                propertiesFile,
                StandardCharsets.UTF_8,
                StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING,
                StandardOpenOption.WRITE
        )) {
            for (final Map.Entry<String, String> entry : values.entrySet()) {
                writer.write(entry.getKey());
                writer.write('=');
                writer.write(entry.getValue());
                writer.newLine();
            }
        }
    }

    private static String format(final double value) {
        return String.format(Locale.US, "%.10f", value);
    }
}