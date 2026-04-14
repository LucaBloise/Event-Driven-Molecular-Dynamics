package ar.edu.itba.sds.tp3.simulation;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

public final class SimulationConfig {
    private static final DateTimeFormatter RUN_TIMESTAMP = DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss");

    private final int particleCount;
    private final double endTime;
    private final double domainDiameter;
    private final double obstacleRadius;
    private final double particleRadius;
    private final double particleMass;
    private final double particleSpeed;
    private final long seed;
    private final int snapshotEveryEvents;
    private final int maxPlacementAttemptsPerParticle;
    private final String outputBaseDirectory;
    private final String outputDirectory;
    private final String runName;

    private SimulationConfig(
            final int particleCount,
            final double endTime,
            final double domainDiameter,
            final double obstacleRadius,
            final double particleRadius,
            final double particleMass,
            final double particleSpeed,
            final long seed,
            final int snapshotEveryEvents,
            final int maxPlacementAttemptsPerParticle,
            final String outputBaseDirectory,
            final String outputDirectory,
            final String runName
    ) {
        this.particleCount = particleCount;
        this.endTime = endTime;
        this.domainDiameter = domainDiameter;
        this.obstacleRadius = obstacleRadius;
        this.particleRadius = particleRadius;
        this.particleMass = particleMass;
        this.particleSpeed = particleSpeed;
        this.seed = seed;
        this.snapshotEveryEvents = snapshotEveryEvents;
        this.maxPlacementAttemptsPerParticle = maxPlacementAttemptsPerParticle;
        this.outputBaseDirectory = outputBaseDirectory;
        this.outputDirectory = outputDirectory;
        this.runName = runName;
    }

    public static SimulationConfig fromArgs(final String[] args) {
        final Map<String, String> options = parseOptions(args);

        final int particleCount = parseInt(options, "n", 100);
        final double endTime = parseDouble(options, "tf", 5.0);
        final double domainDiameter = parseDouble(options, "diameter", 80.0);
        final double obstacleRadius = parseDouble(options, "obstacle-radius", 1.0);
        final double particleRadius = parseDouble(options, "particle-radius", 1.0);
        final double particleMass = parseDouble(options, "mass", 1.0);
        final double particleSpeed = parseDouble(options, "speed", 1.0);
        final long seed = parseLong(options, "seed", System.currentTimeMillis());
        final int snapshotEveryEvents = parseInt(options, "snapshot-every", 1);
        final int maxPlacementAttemptsPerParticle = parseInt(options, "max-placement-attempts", 10_000);
        final String outputBaseDirectory = options.getOrDefault("output-base-dir", "simulation/outputs");
        final String outputDirectory = options.get("output-dir");
        final String runName = options.get("run-name");

        final SimulationConfig config = new SimulationConfig(
                particleCount,
                endTime,
                domainDiameter,
                obstacleRadius,
                particleRadius,
                particleMass,
                particleSpeed,
                seed,
                snapshotEveryEvents,
                maxPlacementAttemptsPerParticle,
                outputBaseDirectory,
                outputDirectory,
                runName
        );
        config.validate();
        return config;
    }

    public static boolean hasHelpFlag(final String[] args) {
        for (final String arg : args) {
            if ("--help".equals(arg) || "-h".equals(arg)) {
                return true;
            }
        }
        return false;
    }

    public static void printUsage() {
        System.out.println("Event-Driven Molecular Dynamics Simulation (System 1)");
        System.out.println();
        System.out.println("Usage:");
        System.out.println("  java ...SimulationMain [options]");
        System.out.println();
        System.out.println("Options:");
        System.out.println("  --n=<int>                      Number of particles (default: 100)");
        System.out.println("  --tf=<double>                  End time in simulated seconds (default: 5.0)");
        System.out.println("  --diameter=<double>            Outer domain diameter in meters (default: 80.0)");
        System.out.println("  --obstacle-radius=<double>     Central obstacle radius in meters (default: 1.0)");
        System.out.println("  --particle-radius=<double>     Mobile particle radius in meters (default: 1.0)");
        System.out.println("  --mass=<double>                Mobile particle mass in kg (default: 1.0)");
        System.out.println("  --speed=<double>               Initial speed modulus in m/s (default: 1.0)");
        System.out.println("  --seed=<long>                  RNG seed (default: current time)");
        System.out.println("  --snapshot-every=<int>         Write state every K events (default: 1)");
        System.out.println("  --max-placement-attempts=<int> Max random placement retries per particle (default: 10000)");
        System.out.println("  --output-base-dir=<path>       Base directory for generated runs (default: simulation/outputs)");
        System.out.println("  --output-dir=<path>            Exact run directory path (overrides output-base-dir/run-name)");
        System.out.println("  --run-name=<name>              Run folder name when output-dir is not provided");
        System.out.println("  --help, -h                     Print this message");
    }

    public int getParticleCount() {
        return particleCount;
    }

    public double getEndTime() {
        return endTime;
    }

    public double getDomainDiameter() {
        return domainDiameter;
    }

    public double getObstacleRadius() {
        return obstacleRadius;
    }

    public double getParticleRadius() {
        return particleRadius;
    }

    public double getParticleMass() {
        return particleMass;
    }

    public double getParticleSpeed() {
        return particleSpeed;
    }

    public long getSeed() {
        return seed;
    }

    public int getSnapshotEveryEvents() {
        return snapshotEveryEvents;
    }

    public int getMaxPlacementAttemptsPerParticle() {
        return maxPlacementAttemptsPerParticle;
    }

    public String getOutputBaseDirectory() {
        return outputBaseDirectory;
    }

    public String getOutputDirectory() {
        return outputDirectory;
    }

    public String getRunName() {
        return runName;
    }

    public double getOuterCollisionRadius() {
        return domainDiameter / 2.0 - particleRadius;
    }

    public double getInnerCollisionRadius() {
        return obstacleRadius + particleRadius;
    }

    public Path resolveRunDirectory() throws IOException {
        final Path runDirectory;

        if (outputDirectory != null && !outputDirectory.isBlank()) {
            runDirectory = Paths.get(outputDirectory);
        } else {
            final String resolvedRunName = runName != null && !runName.isBlank()
                    ? runName
                    : defaultRunName();
            runDirectory = Paths.get(outputBaseDirectory).resolve(resolvedRunName);
        }

        Files.createDirectories(runDirectory);
        return runDirectory;
    }

    public String defaultRunName() {
        final String ts = LocalDateTime.now().format(RUN_TIMESTAMP);
        return String.format(Locale.US, "run_n%d_seed%d_%s", particleCount, seed, ts);
    }

    private void validate() {
        if (particleCount <= 0) {
            throw new IllegalArgumentException("--n must be > 0");
        }
        if (endTime <= 0.0) {
            throw new IllegalArgumentException("--tf must be > 0");
        }
        if (domainDiameter <= 0.0) {
            throw new IllegalArgumentException("--diameter must be > 0");
        }
        if (obstacleRadius < 0.0) {
            throw new IllegalArgumentException("--obstacle-radius must be >= 0");
        }
        if (particleRadius <= 0.0) {
            throw new IllegalArgumentException("--particle-radius must be > 0");
        }
        if (particleMass <= 0.0) {
            throw new IllegalArgumentException("--mass must be > 0");
        }
        if (particleSpeed < 0.0) {
            throw new IllegalArgumentException("--speed must be >= 0");
        }
        if (snapshotEveryEvents <= 0) {
            throw new IllegalArgumentException("--snapshot-every must be > 0");
        }
        if (maxPlacementAttemptsPerParticle <= 0) {
            throw new IllegalArgumentException("--max-placement-attempts must be > 0");
        }

        final double inner = getInnerCollisionRadius();
        final double outer = getOuterCollisionRadius();
        if (inner >= outer) {
            throw new IllegalArgumentException(
                    "Invalid geometry: obstacleRadius + particleRadius must be smaller than diameter/2 - particleRadius"
            );
        }
    }

    private static Map<String, String> parseOptions(final String[] args) {
        final Map<String, String> options = new HashMap<>();

        for (int i = 0; i < args.length; i++) {
            final String arg = args[i];
            if (!arg.startsWith("--")) {
                throw new IllegalArgumentException("Unexpected argument: " + arg);
            }

            final String key;
            final String value;

            final int eq = arg.indexOf('=');
            if (eq >= 0) {
                key = arg.substring(2, eq);
                value = arg.substring(eq + 1);
            } else {
                key = arg.substring(2);
                if (i + 1 < args.length && !args[i + 1].startsWith("--")) {
                    value = args[++i];
                } else {
                    value = "true";
                }
            }

            options.put(key, value);
        }

        return options;
    }

    private static int parseInt(final Map<String, String> options, final String key, final int fallback) {
        if (!options.containsKey(key)) {
            return fallback;
        }
        return Integer.parseInt(options.get(key));
    }

    private static long parseLong(final Map<String, String> options, final String key, final long fallback) {
        if (!options.containsKey(key)) {
            return fallback;
        }
        return Long.parseLong(options.get(key));
    }

    private static double parseDouble(final Map<String, String> options, final String key, final double fallback) {
        if (!options.containsKey(key)) {
            return fallback;
        }
        return Double.parseDouble(options.get(key));
    }
}