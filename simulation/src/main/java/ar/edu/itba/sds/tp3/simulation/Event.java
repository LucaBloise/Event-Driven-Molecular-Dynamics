package ar.edu.itba.sds.tp3.simulation;

import java.util.List;

public final class Event implements Comparable<Event> {
    private final double time;
    private final EventType type;
    private final int aIndex;
    private final int bIndex;
    private final int aCollisionCount;
    private final int bCollisionCount;

    public Event(
            final double time,
            final EventType type,
            final int aIndex,
            final int bIndex,
            final int aCollisionCount,
            final int bCollisionCount
    ) {
        this.time = time;
        this.type = type;
        this.aIndex = aIndex;
        this.bIndex = bIndex;
        this.aCollisionCount = aCollisionCount;
        this.bCollisionCount = bCollisionCount;
    }

    public double getTime() {
        return time;
    }

    public EventType getType() {
        return type;
    }

    public int getAIndex() {
        return aIndex;
    }

    public int getBIndex() {
        return bIndex;
    }

    public boolean isValid(final List<Particle> particles) {
        if (aIndex >= 0 && particles.get(aIndex).getCollisionCount() != aCollisionCount) {
            return false;
        }
        if (bIndex >= 0 && particles.get(bIndex).getCollisionCount() != bCollisionCount) {
            return false;
        }
        return true;
    }

    @Override
    public int compareTo(final Event other) {
        final int byTime = Double.compare(this.time, other.time);
        if (byTime != 0) {
            return byTime;
        }

        final int byType = this.type.compareTo(other.type);
        if (byType != 0) {
            return byType;
        }

        final int byA = Integer.compare(this.aIndex, other.aIndex);
        if (byA != 0) {
            return byA;
        }

        return Integer.compare(this.bIndex, other.bIndex);
    }
}