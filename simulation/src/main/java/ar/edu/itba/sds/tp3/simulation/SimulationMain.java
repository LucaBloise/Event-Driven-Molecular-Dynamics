package ar.edu.itba.sds.tp3.simulation;

import java.nio.file.Path;
import java.util.Locale;

public final class SimulationMain {

    public static void main(final String[] args) {
        if (SimulationConfig.hasHelpFlag(args)) {
            SimulationConfig.printUsage();
            return;
        }

        try {
            final SimulationConfig config = SimulationConfig.fromArgs(args);
            final Path runDirectory = config.resolveRunDirectory();
            final Path outputFile = config.shouldWriteDeltaOutputEvents()
                    ? runDirectory.resolve("output.txt")
                    : null;
            final Path propertiesFile = runDirectory.resolve("properties.txt");

            final EventDrivenSimulation simulation = new EventDrivenSimulation(config);
            final SimulationResult result = simulation.run(outputFile);

            SimulationPropertiesWriter.write(propertiesFile, config, result, outputFile);

            System.out.println("Simulation finished successfully.");
            System.out.println("Run directory: " + runDirectory.toAbsolutePath());
            if (outputFile != null) {
                System.out.println("Output file: " + outputFile.toAbsolutePath() + " (mode: delta)");
            } else {
                System.out.println("Output file: <disabled by --no-output>");
            }
            System.out.println("Properties file: " + propertiesFile.toAbsolutePath());
            System.out.printf(Locale.US, "Processed events: %d%n", result.getProcessedEvents());
            System.out.printf(Locale.US, "Written frames: %d%n", result.getWrittenFrames());
            System.out.printf(Locale.US, "Final simulated time (s): %.6f%n", result.getFinalSimulatedTime());
            System.out.printf(Locale.US, "Execution time (s): %.6f%n", result.getExecutionSeconds());
        } catch (final Exception exception) {
            System.err.println("Simulation failed: " + exception.getMessage());
            exception.printStackTrace(System.err);
            System.exit(1);
        }
    }
}