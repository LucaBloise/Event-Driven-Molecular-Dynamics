package ar.edu.itba.sds.tp3.simulation;

public enum ParticleState {
    FRESH(0, 200, 0),
    USED(148, 0, 211);

    private final int red;
    private final int green;
    private final int blue;

    ParticleState(final int red, final int green, final int blue) {
        this.red = red;
        this.green = green;
        this.blue = blue;
    }

    public int getRed() {
        return red;
    }

    public int getGreen() {
        return green;
    }

    public int getBlue() {
        return blue;
    }
}