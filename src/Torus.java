import java.awt.*;

public class Torus extends SceneObject{
    public static final double EPSILON = 1.0E-7f;

    // torus parameters
    private Vector3 position;
    private Vector3 axis;
    private double radius_torus;
    private double radius_tube;

    // surface coefficients
    private final double TORUS_KD = 0.8;
    private final double TORUS_KS = 1.2;
    private final double TORUS_ALPHA = 10;
    private final double TORUS_REFLECTIVITY = 0.3;

    public Torus(Vector3 position, Vector3 axis, double radius_torus, double radius_tube, ColorRGB colour){
        this.position = position;
        this.axis = axis.normalised();
        this.radius_torus = radius_torus;
        this.radius_tube = radius_tube;
        this.colour = colour;
        this.phong_kD = TORUS_KD;
        this.phong_kS = TORUS_KS;
        this.phong_alpha = TORUS_ALPHA;
        this.reflectivity = TORUS_REFLECTIVITY;
    }

    public Torus(Vector3 position, Vector3 axis, double radius_torus, double radius_tube, ColorRGB colour,
                 double kD, double kS, double alphaS, double reflectivity) {
        this.position = position;
        this.axis = axis.normalised();
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
        Vector3 A = this.axis;

        // get ray parameters
        Vector3 E = ray.getOrigin();
        Vector3 D = ray.getDirection();

        // transform the ray (because equation only works for torus centering at origin)
        Vector3 centerToRayOrigin = E.subtract(C);
        double centerToRayOriginDotDirection = D.dot(centerToRayOrigin);
        double	centerToRayOriginDotDirectionSquared = centerToRayOrigin.dot(centerToRayOrigin);
        double innerRadiusSquared = r * r;
        double outerRadiusSquared = R * R;
        double axisDotCenterToRayOrigin	= A.dot(centerToRayOrigin);
        double axisDotRayDirection = A.dot(D);
        double a = 1 - axisDotRayDirection * axisDotRayDirection;
        double b = 2 * (centerToRayOrigin.dot(D) - axisDotCenterToRayOrigin * axisDotRayDirection);
        double c = centerToRayOriginDotDirectionSquared - axisDotCenterToRayOrigin * axisDotCenterToRayOrigin;
        double d = centerToRayOriginDotDirectionSquared + outerRadiusSquared - innerRadiusSquared;

        // calculate the coefficients for c4(x^4) + c3(c^3) + c2(c^2) + c1(c) + c0 = 0
        double[] coefficients = new double[5];
        coefficients[4] = 1;
        coefficients[3] = 4 * centerToRayOriginDotDirection;
        coefficients[2] = 2 * d + coefficients[3] * coefficients[3] * 0.25f - 4 * outerRadiusSquared * a;
        coefficients[1] = coefficients[3] * d - 4 * outerRadiusSquared * b;
        coefficients[0] = d * d - 4 * outerRadiusSquared * c;

        // calculate the value(s) of t using QuarticRootFinder
        QuarticRootFinder qr = new QuarticRootFinder();
        double[] t = qr.QuarticRoot(coefficients);

        // find the nearest (smallest), non-negative t
        double min_t = Double.MAX_VALUE;
        for (int i = 0; i < t.length; i ++) {
            if (t[i] > EPSILON && min_t - t[i] > EPSILON) min_t = t[i];
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
    public Vector3 getNormalAt(Vector3 point) {
        // calculate the normal of the torus surface at a given position
        Vector3 centerToPoint = point.subtract(this.position);
        double centerToPointDotAxis = centerToPoint.dot(this.axis);
        Vector3 direction = centerToPoint.subtract(this.axis.scale(centerToPointDotAxis));
        direction = direction.normalised();
        Vector3 normal = point.subtract(this.position).add(direction.scale(this.radius_torus));
        return normal.normalised();
    }
}
