import java.awt.*;

public class Torus extends SceneObject{

    // torus parameters
    private Vector3 position;
    private double radius_torus;
    private double radius_tube;

    // surface coefficients
    private final double TORUS_KD = 0.8;
    private final double TORUS_KS = 1.2;
    private final double TORUS_ALPHA = 10;
    private final double TORUS_REFLECTIVITY = 0.3;

    public Torus(Vector3 position, double radius_torus, double radius_tube, ColorRGB colour){
        this.position = position;
        this.radius_torus = radius_torus;
        this.radius_tube = radius_tube;
        this.colour = colour;
        this.phong_kD = TORUS_KD;
        this.phong_kS = TORUS_KS;
        this.phong_alpha = TORUS_ALPHA;
        this.reflectivity = TORUS_REFLECTIVITY;
    }

    public Torus(Vector3 position, double radius_torus, double radius_tube, ColorRGB colour,
                 double kD, double kS, double alphaS, double reflectivity) {
        this.position = position;
        this.radius_torus = radius_torus;
        this.radius_tube = radius_tube;
        this.colour = colour;
        this.phong_kD = kD;
        this.phong_kS = kS;
        this.phong_alpha = alphaS;
        this.reflectivity = reflectivity;
    }

    @Override
    public RaycastHit intersectionWith(Ray ray) {
        // get torus parameters
        double R = this.radius_torus;
        double r = this.radius_tube;
        Vector3 C = this.position;

        // get ray parameters
        Vector3 E = ray.getOrigin().subtract(C);
        Vector3 D = ray.getDirection();

        double xe = E.x;
        double ye = E.y;
        double ze = E.z;
        double xd = D.x;
        double yd = D.y;
        double zd = D.z;

        // calculate coefficients for t: c4(t^4) + c3(t^3) + c2(t^2) + c1(t) + c0 = 0
        double[] c = new double[5];
        double base1 = xd * xd + yd * yd + zd * zd;
        double base2 = xe * xd + ye * yd + ze * zd;
        double base3 = xe * xe + ye * ye + ze * ze - (r * r + R * R);
        c[4] = base1 * base1;
        c[3] = 4.0 * base1 * base2;
        c[2] = 2.0 * base1 * base3 + 4.0 * base2 * base2 + 4.0 * R * R * yd * yd;
        c[1] = 4.0 * base3 * base2 + 8.0 * R * R * ye * yd;
        c[0] = base3 * base3 - 4.0 * R * R * (r * r - ye * ye);

        // calculate the value(s) of t using QuarticRootFinder
        QuarticRootFinder qr = new QuarticRootFinder();
        double[] t = qr.QuarticRoot(c);

        // find the nearest (smallest), non-negative t
        double min_t = Double.MAX_VALUE;
        for (int i = 0; i < t.length; i ++) {
            if (t[i] > 0) min_t = Math.min(min_t, t[i]);
        }

        // calculate distance, intersection location, normal
        RaycastHit rc = new RaycastHit();
        if (min_t != Double.MAX_VALUE){
            Vector3 location = ray.evaluateAt(min_t);
            rc = new RaycastHit(this, min_t, location, this.getNormalAt(location));
        }

        return rc;
    }

    @Override
    public Vector3 getNormalAt(Vector3 position) {
        // calculate the normal of the torus surface at a given position
        double x = position.x;
        double y = position.y;
        double z = position.z;
        double base = (x * x + y * y + z * z) - (radius_torus * radius_torus + radius_tube * radius_torus);
        return new Vector3(4.0 * x * base,
                           4.0 * y * (base + 2.0 * radius_torus * radius_tube),
                           4.0 * z * base).normalised();
    }
}
