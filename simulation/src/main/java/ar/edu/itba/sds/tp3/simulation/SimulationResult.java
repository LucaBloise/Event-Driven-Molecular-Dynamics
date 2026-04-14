package ar.edu.itba.sds.tp3.simulation;

public final class SimulationResult {
    private final double finalSimulatedTime;
    private final long processedEvents;
    private final long writtenFrames;
    private final double executionSeconds;

    public SimulationResult(
            final double finalSimulatedTime,
            final long processedEvents,
            final long writtenFrames,
            final double executionSeconds
    ) {
        this.finalSimulatedTime = finalSimulatedTime;
        this.processedEvents = processedEvents;
        this.writtenFrames = writtenFrames;
        this.executionSeconds = executionSeconds;
    }

    public double getFinalSimulatedTime() {
        return finalSimulatedTime;
    }

    public long getProcessedEvents() {
        return processedEvents;
    }

    public long getWrittenFrames() {
        return writtenFrames;
    }

    public double getExecutionSeconds() {
        return executionSeconds;
    }
}