public class QuarticRootFinder {

    /*** Ferrari Method ***/
    private final double EQN_EPS = 1e-9;

    private boolean IsZero(double x){
        return ((x) > -EQN_EPS && (x) < EQN_EPS);
    }

    public double[] QuarticRoot(double c[]){
        double[] coeffs = new double[4];
        double[] solution;

        // normalised form: x^4 + Ax^3 + Bx^2 + Cx + D = 0
        double A = c[3] / c[4];
        double B = c[2] / c[4];
        double C = c[1] / c[4];
        double D = c[0] / c[4];

        // substitute x = y - A/4 to eliminate cubic term: x^4 + px^2 + qx + r = 0
        double squareA = A * A;
        double p = - 3.0 / 8 * squareA + B;
        double q = 1.0 / 8.0 * squareA * A - 1.0 / 2.0 * A * B + C;
        double r = - 3.0 / 256.0 * squareA * squareA + 1.0 / 16.0 * squareA * B - 1.0 / 4.0 * A * C + D;

        if (IsZero(r)){
            // no absolute term: y(y^3 + py + q) = 0
            coeffs[0] = q;
            coeffs[1] = p;
            coeffs[2] = 0;
            coeffs[3] = 1;
            double[] ans = CubicRoot(coeffs);
            solution = new double[ans.length+1];
            for (int i = 0; i < ans.length; i ++) solution[i] = ans[i];
            solution[ans.length] = 0;

        } else{
            // solve the resolvent cubic ...
            coeffs[ 0 ] = 1.0 / 2.0 * r * p - 1.0 / 8.0 * q * q;
            coeffs[ 1 ] = - r;
            coeffs[ 2 ] = - 1.0 / 2.0 * p;
            coeffs[ 3 ] = 1;

            double[] psolution;
            psolution = CubicRoot(coeffs);

            // ... and take the one real solution ...
            double z,u,v;
            z = psolution[ 0 ];

            // ... to build two quadric equations
            u = z * z - r;
            v = 2 * z - p;

            if (IsZero(u)) { u = 0; }
            else if (u > 0) { u = Math.sqrt(u); }
            else{ return new double[]{}; }

            if (IsZero(v)) { v = 0; }
            else if (v > 0) { v = Math.sqrt(v); }
            else { return new double[]{}; }

            coeffs[ 0 ] = z - u;
            coeffs[ 1 ] = q < 0 ? -v : v;
            coeffs[ 2 ] = 1;

            psolution = QuadricRoot(coeffs);

            coeffs[ 0 ]= z + u;
            coeffs[ 1 ] = q < 0 ? v : -v;
            coeffs[ 2 ] = 1;

            double[] qsolution = QuadricRoot(coeffs);

            // concatenate psolution + qsolution
            solution = new double[psolution.length + qsolution.length];
            for (int i = 0; i < psolution.length; i ++) solution[i] = psolution[i];
            for (int i = 0; i < qsolution.length; i ++) solution[i+psolution.length] = qsolution[i];
        }

        /* resubstitute */

        double sub = 1.0/4 * A;
        for (int i = 0; i < solution.length; i ++)
            solution[i] -= sub;

        return solution;
    }


    private double[] CubicRoot(double c[]){
        double[] solution;

        // normalised form: x^3 + Ax^2 + Bx + C = 0
        double A = c[2] / c[3];
        double B = c[1] / c[3];
        double C = c[0] / c[3];

        // substitute x = y - A/3 to eliminate quadric term: x^3 +px + q = 0
        double squareA = A * A;
        double p = 1.0 / 3 * (- 1.0 / 3 * squareA + B);
        double q = 1.0 / 2 * (2.0 / 27 * A * squareA - 1.0 / 3 * A * B + C);

        /* use Cardano's formula */
        double cubicp = p * p * p;
        double D = q * q + cubicp;

        if (IsZero(D)){
            if (IsZero(q)){       // one triple solution
                solution = new double[]{0};
            }
            else{                 // one single and one double solution
                double u = Math.cbrt(-q);
                solution = new double[]{2*u, -u};
            }
        }
        else if (D < 0){          // Casus irreducibilis: three real solutions
            double phi = 1.0/3 * Math.acos(-q / Math.sqrt(-cubicp));
            double t = 2 * Math.sqrt(-p);
            solution = new double[]{t * Math.cos(phi), - t * Math.cos(phi + Math.PI / 3), - t * Math.cos(phi - Math.PI / 3)};
        }
        else{                     // one real solution
            double sqrtD = Math.sqrt(D);
            double u = Math.cbrt(sqrtD - q);
            double v = - Math.cbrt(sqrtD + q);
            solution = new double[]{u + v};
        }

        /* resubstitute */
        double sub = 1.0 / 3 * A;
        for (int i = 0; i < solution.length; ++i)
            solution[i] -= sub;

        return solution;
    }


    private double[] QuadricRoot(double c[]){
        double[] solution;

        // normalised form: x^2 + px + q = 0
        double p = c[1] / (2 * c[2]);
        double q = c[0] / c[2];

        // calculate determinant
        double D = p * p - q;

        if (IsZero(D)){
            solution = new double[]{-p};
        }
        else if (D < 0){
            solution = new double[]{};
        }
        else{
            double sqrtD = Math.sqrt(D);
            solution = new double[]{sqrtD-p, -sqrtD-p};
        }

        return solution;
    }

}







