package ar.edu.itba.sds.tp3.simulation;

public final class Particle {
    public static final double EPS = 1.0e-9;

    private final int id;
    private final double mass;
    private final double radius;

    private double x;
    private double y;
    private double vx;
    private double vy;
    private int collisionCount;
    private ParticleState state;

    public Particle(
            final int id,
            final double x,
            final double y,
            final double vx,
            final double vy,
            final double mass,
            final double radius,
            final ParticleState state
    ) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.vx = vx;
        this.vy = vy;
        this.mass = mass;
        this.radius = radius;
        this.state = state;
        this.collisionCount = 0;
    }

    public int getId() {
        return id;
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public double getVx() {
        return vx;
    }

    public double getVy() {
        return vy;
    }

    public double getMass() {
        return mass;
    }

    public double getRadius() {
        return radius;
    }

    public int getCollisionCount() {
        return collisionCount;
    }

    public ParticleState getState() {
        return state;
    }

    public void setState(final ParticleState state) {
        this.state = state;
    }

    public void advance(final double dt) {
        x += vx * dt;
        y += vy * dt;
    }

    public void snapToCircle(final double circleRadius) {
        final double norm = Math.hypot(x, y);
        if (norm <= EPS) {
            return;
        }
        final double scale = circleRadius / norm;
        x *= scale;
        y *= scale;
    }

    public double timeToHitParticle(final Particle other) {
        if (this == other) {
            return Double.POSITIVE_INFINITY;
        }

        final double dx = other.x - this.x;
        final double dy = other.y - this.y;
        final double dvx = other.vx - this.vx;
        final double dvy = other.vy - this.vy;

        final double dvdr = dx * dvx + dy * dvy;
        if (dvdr >= 0.0) {
            return Double.POSITIVE_INFINITY;
        }

        final double dvdv = dvx * dvx + dvy * dvy;
        if (dvdv <= EPS) {
            return Double.POSITIVE_INFINITY;
        }

        final double drdr = dx * dx + dy * dy;
        final double sigma = this.radius + other.radius;
        final double discriminant = dvdr * dvdr - dvdv * (drdr - sigma * sigma);
        if (discriminant < 0.0) {
            return Double.POSITIVE_INFINITY;
        }

        final double dt = -(dvdr + Math.sqrt(discriminant)) / dvdv;
        return dt > EPS ? dt : Double.POSITIVE_INFINITY;
    }

    public double timeToHitCircle(final double circleRadius) {
        final double a = vx * vx + vy * vy;
        if (a <= EPS) {
            return Double.POSITIVE_INFINITY;
        }

        final double b = 2.0 * (x * vx + y * vy);
        final double c = x * x + y * y - circleRadius * circleRadius;
        final double discriminant = b * b - 4.0 * a * c;
        if (discriminant < 0.0) {
            return Double.POSITIVE_INFINITY;
        }

        final double sqrt = Math.sqrt(discriminant);
        final double t1 = (-b - sqrt) / (2.0 * a);
        final double t2 = (-b + sqrt) / (2.0 * a);

        double dt = Double.POSITIVE_INFINITY;
        if (t1 > EPS) {
            dt = t1;
        }
        if (t2 > EPS && t2 < dt) {
            dt = t2;
        }
        return dt;
    }

    public void bounceOffParticle(final Particle other) {
        final double dx = other.x - this.x;
        final double dy = other.y - this.y;
        final double dvx = other.vx - this.vx;
        final double dvy = other.vy - this.vy;
        final double dvdr = dx * dvx + dy * dvy;
        final double dist = this.radius + other.radius;

        final double magnitude =
                (2.0 * this.mass * other.mass * dvdr) / ((this.mass + other.mass) * dist);
        final double fx = magnitude * dx / dist;
        final double fy = magnitude * dy / dist;

        this.vx += fx / this.mass;
        this.vy += fy / this.mass;
        other.vx -= fx / other.mass;
        other.vy -= fy / other.mass;

        this.collisionCount++;
        other.collisionCount++;
    }

    public void bounceOffCircularBoundary() {
        final double norm = Math.hypot(x, y);
        if (norm <= EPS) {
            return;
        }

        final double nx = x / norm;
        final double ny = y / norm;
        final double dot = vx * nx + vy * ny;

        vx -= 2.0 * dot * nx;
        vy -= 2.0 * dot * ny;
        collisionCount++;
    }
}