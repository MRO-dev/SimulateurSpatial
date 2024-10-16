package org.example.leoplaw.tolerance;

public class ToleranceBuilder {

    private final double firstParam;
    private final double secondParam;
    private final double thirdParam;
    private final double fourthParam;
    private final double fifthParam;

    public ToleranceBuilder() {
        this.firstParam  = -1;
        this.secondParam = -1;
        this.thirdParam  = -1;
        this.fourthParam = -1;
        this.fifthParam  = -1;
    }

    private ToleranceBuilder(final double firstParam, final double secondParam, final double thirdParam,
                             final double fourthParam, final double fifthParam) {
        this.firstParam  = firstParam;
        this.secondParam = secondParam;
        this.thirdParam  = thirdParam;
        this.fourthParam = fourthParam;
        this.fifthParam  = fifthParam;
    }

    public ToleranceBuilder withFirstParam(final double firstParam) {
        return new ToleranceBuilder(firstParam, secondParam, thirdParam, fourthParam, fifthParam);
    }

    public ToleranceBuilder withSecondParam(final double secondParam) {
        return new ToleranceBuilder(firstParam, secondParam, thirdParam, fourthParam, fifthParam);
    }

    public ToleranceBuilder withThirdParam(final double thirdParam) {
        return new ToleranceBuilder(firstParam, secondParam, thirdParam, fourthParam, fifthParam);
    }

    public ToleranceBuilder withFourthParam(final double fourthParam) {
        return new ToleranceBuilder(firstParam, secondParam, thirdParam, fourthParam, fifthParam);
    }

    public ToleranceBuilder withFifthParam(final double fifthParam) {
        return new ToleranceBuilder(firstParam, secondParam, thirdParam, fourthParam, fifthParam);
    }

    public KeplerianTolerance buildKeplerianTolerance() {
        double[] tolerance = { firstParam, secondParam, thirdParam, fourthParam, fifthParam };

        // Default Semi Axis Tolerance
        if (firstParam <= 0) {
            tolerance[0] = 10e5;

        }
        // Default eccentricity tolerance
        if (secondParam <= 0 || secondParam >= 1) {
            tolerance[1] = 1.5;
        }

        // Default inclination tolerance
        if (thirdParam <= 0 || thirdParam > 180) {
            tolerance[2] = 360;
        }

        // Default w tolerance
        if (fourthParam <= 0 || fourthParam > 360) {
            tolerance[3] = 360;
        }

        // Default W tolerance
        if (fifthParam <= 0 || fifthParam > 360) {
            tolerance[4] = 360;
        }

        return new KeplerianTolerance(tolerance);
    }

    public EquinoctialTolerance buildEquinoctialTolerance() {
        double[] tolerance = { firstParam, secondParam, thirdParam, fourthParam, fifthParam };

        // Default Semi Axis Tolerance
        if (firstParam <= 0) {
            tolerance[0] = 10e5;

        }
        // Default ex tolerance
        if (secondParam <= 0 || secondParam >= 1) {
            tolerance[1] = 1.5;
        }

        // Default ey tolerance
        if (thirdParam <= 0 || thirdParam >= 1) {
            tolerance[2] = 1.5;
        }

        // Default hx
        if (fourthParam <= 0 || fourthParam > 100) {
            tolerance[3] = 100;
        }

        // Default W tolerance
        if (fifthParam <= 0 || fourthParam > 100) {
            tolerance[4] = 100;
        }

        return new EquinoctialTolerance(tolerance);
    }

    public abstract class Tolerance {

        final double[] tolerances;

        private Tolerance(final double[] tolerances) {
            this.tolerances = tolerances;
        }

        public double[] getTolerances() {
            return tolerances;
        }

        public int isFirstParameterCriterion() {

            if (tolerances[0] < 5e6) {
                return 1;
            }
            else {
                return 0;
            }
        }

        public abstract int isSecondParameterCriterion();

        public abstract int isThirdParameterCriterion();

        public abstract int isFourthParameterCriterion();

        public abstract int isFifthParameterCriterion();
    }

    public class KeplerianTolerance extends Tolerance {

        protected KeplerianTolerance(final double[] tolerances) {
            super(tolerances);
        }

        // Return 0 if there is no criterion convergence for eccentricity , 1 otherwise
        public int isSecondParameterCriterion() {

            if (tolerances[1] < 1) {
                return 1;
            }
            else {
                return 0;
            }

        }

        // Return 0 if there is no criterion convergence for inclination , 1 otherwise
        public int isThirdParameterCriterion() {
            if (tolerances[2] < 180) {
                return 1;
            }
            else {
                return 0;
            }
        }

        // Return 0 if there is no criterion convergence for perigee argument , 1 otherwise
        public int isFourthParameterCriterion() {
            if (tolerances[3] < 360) {
                return 1;
            }
            else {
                return 0;
            }
        }

        // Return 0 if there is no criterion convergence for ascending node argument, 1 otherwise
        public int isFifthParameterCriterion() {
            if (tolerances[4] < 360) {
                return 1;
            }
            else {
                return 0;
            }
        }
    }

    public class EquinoctialTolerance extends Tolerance {

        private EquinoctialTolerance(final double[] tolerances) {
            super(tolerances);
        }

        // Return 0 if there is no criterion convergence for eccentricity , 1 otherwise
        public int isSecondParameterCriterion() {

            if (tolerances[1] < 1) {
                return 1;
            }
            else {
                return 0;
            }
        }

        // Return 0 if there is no criterion convergence for inclination , 1 otherwise
        public int isThirdParameterCriterion() {
            if (tolerances[2] < 1) {
                return 1;
            }
            else {
                return 0;
            }
        }

        // Return 0 if there is no criterion convergence for perigee argument , 1 otherwise
        public int isFourthParameterCriterion() {
            if (tolerances[3] <= 50) {
                return 1;
            }
            else {
                return 0;
            }
        }

        // Return 0 if there is no criterion convergence for ascending node argument, 1 otherwise
        public int isFifthParameterCriterion() {
            if (tolerances[4] <= 50) {
                return 1;
            }
            else {
                return 0;
            }

        }

    }
}



