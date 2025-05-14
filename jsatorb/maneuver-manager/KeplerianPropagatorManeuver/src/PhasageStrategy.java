import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.orekit.data.DirectoryCrawler;
import org.orekit.forces.maneuvers.ImpulseManeuver;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.orekit.propagation.events.DateDetector;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;

import java.io.IOException;
import java.text.ParseException;
import java.util.Map;

/**
 * This class implements a phasing maneuver strategy.
 * It performs two impulses to modify a spacecraft's phase angle within the same orbit.
 */
public class PhasageStrategy extends AbstractManeuverStrategy {
    private double DELTA_THETA;  // desired phase change (degrees)

    /**
     * Constructor for PhasageStrategy
     *
     * @param sma Semi-major axis of the initial orbit (m)
     * @param ecc Eccentricity of the initial orbit
     * @param incDeg Inclination of the initial orbit (degrees)
     * @param raanDeg Right ascension of ascending node (degrees)
     * @param aopDeg Argument of perigee (degrees)
     * @param meanAnomDeg Mean anomaly (degrees)
     * @param dryMass Dry mass of the spacecraft (kg)
     * @param ergol Propellant mass (kg)
     * @param isp Specific impulse (s)
     * @param delta_theta Desired phase angle change (degrees)
     * @param dateStr Initial date
     * @param modeParameter Mode parameter (e.g., "blue1")
     * @param commandParameter Command parameter
     * @throws IOException If file operations fail
     */
    public PhasageStrategy(
            double sma, double ecc, double incDeg, double raanDeg, double aopDeg, double meanAnomDeg,
            double dryMass, double ergol, double isp,
            double delta_theta,
            String dateStr, String modeParameter, String commandParameter) throws IOException {
        super(modeParameter, commandParameter);
        this.SMA    = sma;
        this.ECC    = ecc;
        this.INC    = FastMath.toRadians(incDeg);
        this.RAAN   = FastMath.toRadians(raanDeg);
        this.PA     = FastMath.toRadians(aopDeg);
        this.ANO    = FastMath.toRadians(meanAnomDeg);
        this.DRYMASS = dryMass;
        this.ERGOL   = ergol;
        this.ISP     = isp;
        this.DELTA_THETA = delta_theta;
        this.DATE    = dateStr;
        System.out.println("Phase angle target: " + delta_theta + "°");
    }

    /**
     * Calculate final mass after a ΔV using the Tsiolkovsky rocket equation.
     */
    public static double calculateFinalMass(double initialMass, double deltaV, double ISP, double g0) {
        return initialMass * FastMath.exp(-Math.abs(deltaV) / (ISP * g0));
    }

    /**
     * Validate orbit and mass parameters before the maneuver.
     */
    private void validateOrbitParameters() {
        if (SMA <= 0) {
            throw new IllegalArgumentException("Semi-major axis must be positive");
        }
        if (DRYMASS <= 0 || ERGOL < 0) {
            throw new IllegalArgumentException("Invalid mass parameters");
        }
    }

    /**
     * Validate maneuver feasibility (time order and propellant).
     */
    private void validateManeuverParameters(AbsoluteDate initialDate, AbsoluteDate apsideDate, double finalMass) {
        if (!isEqualOrAfterWithTolerance(initialDate, apsideDate, TIME_TOLERANCE_SECONDS)) {
            throw new IllegalArgumentException("Initial date must be >= reference apside date");
        }
        if (finalMass < DRYMASS || ERGOL < 1e-3) {
            throw new IllegalArgumentException("Insufficient propellant for phasing maneuver");
        }
        if (ECC > 0.05) {
            throw new IllegalArgumentException("Eccentricity too high for accurate phasing calculation");
        }
    }

    /**
     * Normalize angle to [0, 2π)
     */
    private static double normalizeAngle0to2pi(double angle) {
        angle = angle % (2*Math.PI);
        if (angle < 0) angle += 2*Math.PI;
        return angle;
    }

    /**
     * Get the mean longitude (ω + M)
     */
    private static double meanLongitude(KeplerianOrbit o) {
        return normalizeAngle0to2pi(o.getPerigeeArgument() + o.getMeanAnomaly());
    }

    /**
     * Print detailed state information for debugging
     */
    private static void logState(final String label, final SpacecraftState s) {
        KeplerianOrbit o = new KeplerianOrbit(s.getOrbit());

        // Normalize mean and true anomaly to [0, 2π)
        double meanAnomaly = normalizeAngle0to2pi(o.getMeanAnomaly());
        double trueAnomaly = normalizeAngle0to2pi(o.getTrueAnomaly());

        System.out.printf(
                "\n[%s]\n" +
                        "  Epoch   : %s\n" +
                        "  Mass    : %.3f kg\n" +
                        "  a       : %.3f km\n" +
                        "  e       : %.8f\n" +
                        "  i       : %.3f °\n" +
                        "  RAAN    : %.3f °\n" +
                        "  ω (PA)  : %.3f °\n" +
                        "  M       : %.3f °\n" +
                        "  ν (TA)  : %.3f °\n",
                label,
                s.getDate(),
                s.getMass(),
                o.getA() / 1000.0,
                o.getE(),
                FastMath.toDegrees(o.getI()),
                FastMath.toDegrees(o.getRightAscensionOfAscendingNode()),
                FastMath.toDegrees(o.getPerigeeArgument()),
                FastMath.toDegrees(meanAnomaly),
                FastMath.toDegrees(trueAnomaly)
        );
    }

    @Override
    public void loadMassData(Boolean isMassCalculation) throws IOException {
        super.loadMassData(isMassCalculation);
        if (isMassCalculation) {
            if (cmdData.containsKey("DELTA_THETA")) {
                this.DELTA_THETA = Double.parseDouble(cmdData.get("DELTA_THETA"));
            }
        }
    }

    @Override
    public void loadTimeData(Boolean isTimeCalculation) throws IOException {
        super.loadTimeData(isTimeCalculation);
        if (tipData.containsKey("DELTA_THETA")) {
            this.DELTA_THETA = Double.parseDouble(tipData.get("DELTA_THETA"));
        }
        System.out.println(" - timeIntermediateParametersFile => ManeuverType: " + MANEUV_TYPE
                + ", PhaseAngle: " + DELTA_THETA + "°");
    }

    @Override
    public void computeAndExecute() throws IOException, ParseException {
        /* 0) --- House-keeping ------------------------------------------------ */
        final long wallStart = System.currentTimeMillis();
        manager.addProvider(new DirectoryCrawler(orekitData));
        validateOrbitParameters();

        /* 1) --- Frames & dates ---------------------------------------------- */
        final Frame inertial = FramesFactory.getEME2000();
        final MqttService mqttService = new MqttService();
        final AbsoluteDate apsisDate = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
        final AbsoluteDate tleDate = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
        final AbsoluteDate initialDate = tleDate.shiftedBy(manoeuverRelativeDate);
        final AbsoluteDate endHorizonDate = new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());
        final double g0 = 9.80665;

        final double m0 = DRYMASS + ERGOL;
        System.out.printf("%n=== PHASING MANOEUVRE  ==================================================%n");
        System.out.printf("Dry mass              : %,.3f kg%n", DRYMASS);
        System.out.printf("Propellant on‑board   : %,.3f kg%n", ERGOL);
        System.out.printf("ISP                   : %,.1f s%n", ISP);
        System.out.printf("Initial wet mass m0   : %,.3f kg%n", m0);
        System.out.printf("Reference apsis date  : %s%n", apsisDate);
        System.out.printf("Burn‑1 epoch (t0)     : %s%n", initialDate);

        /* 2) --- Initial circular state -------------------------------------- */
        final KeplerianOrbit orbitAtApsis = new KeplerianOrbit(
                SMA, ECC, INC, PA, RAAN, ANO, PositionAngle.MEAN,
                inertial, apsisDate, MU);
        SpacecraftState stateApsis = new SpacecraftState(orbitAtApsis, m0);
        logState("state@apsis", stateApsis);

        final KeplerianPropagator prepProp = new KeplerianPropagator(orbitAtApsis);
        final SpacecraftState stateInit = new SpacecraftState(
                prepProp.propagate(initialDate).getOrbit(), m0);
        logState("state@initial", stateInit);

        /* 3) --- Phase still to perform -------------------------------------- */
        final double nTgt = FastMath.sqrt(MU / FastMath.pow(SMA, 3));
        final double drift = nTgt * initialDate.durationFrom(apsisDate);

        final double dThetaCmd = FastMath.toRadians(DELTA_THETA);
        final double dTheta = MathUtils.normalizeAngle(dThetaCmd, 0);
        System.out.printf("Δθ commanded          : %+.4f deg%n", FastMath.toDegrees(dThetaCmd));
        System.out.printf("Natural drift         : %+.4f deg%n", FastMath.toDegrees(drift));
        System.out.printf("Δθ still to perform   : %+.4f deg%n", FastMath.toDegrees(dTheta));

        /* 4) --- First-guess phasing period ---------------------------------- */
        final int kInt =  1;
        double Tph = (kInt * 2 * FastMath.PI - dTheta) / nTgt;   // s
        System.out.printf("k_int                 : %d%n", kInt);
        System.out.printf("T_ph (initial guess)  : %.3f min%n", Tph/60);

        /* 5) --- Phasing-orbit geometry -------------------------------------- */
        double aPh = FastMath.cbrt(MU * FastMath.pow(Tph/(2*FastMath.PI), 2));
        boolean outer = aPh > SMA;
        double rp = outer ? SMA : 2*aPh - SMA;
        double ra = outer ? 2*aPh - SMA : SMA;
        double ePh = (ra - rp) / (ra + rp);

        System.out.printf("Phasing orbit (first) : a=%.3f km  e=%.6f  rp=%.3f km  ra=%.3f km%n",
                aPh/1e3, ePh, rp/1e3, ra/1e3);

        /* 6) --- ΔV budget ---------------------------------------------------- */
        double vCirc = FastMath.sqrt(MU / SMA);
        double vPeriPh = FastMath.sqrt(MU*(2/SMA - 1/aPh));
        double dv1 = vPeriPh - vCirc;
        double dv2 = -dv1;
        double m1 = calculateFinalMass(m0, dv1, ISP, g0);
        double m2 = calculateFinalMass(m1, dv2, ISP, g0);
        System.out.printf("ΔV1 / ΔV2             : %+.3f / %+.3f m/s%n", dv1, dv2);
        System.out.printf("Propellant expected   : %.3f kg%n", m0-m2);

        /* 7) --- Burn-1 ------------------------------------------------------- */
        final Vector3D dir1 = stateInit.getPVCoordinates().getVelocity()
                .normalize().scalarMultiply(dv1);
        final DateDetector firstBurnTrigger = new DateDetector(initialDate.shiftedBy(0.001));
        final ImpulseManeuver firstBurn = new ImpulseManeuver(firstBurnTrigger, dir1, ISP);

        final KeplerianPropagator prop1 = new KeplerianPropagator(stateInit.getOrbit());
        prop1.addEventDetector(firstBurn);

        SpacecraftState stateAfterBurn1 = prop1.propagate(initialDate.shiftedBy(0.002));
        stateAfterBurn1 = new SpacecraftState(stateAfterBurn1.getOrbit(), m1);
        logState("state@afterBurn1", stateAfterBurn1);

        final Vector3D rB1 = stateAfterBurn1.getPVCoordinates().getPosition();
        final Vector3D h = Vector3D.crossProduct(rB1,
                stateAfterBurn1.getPVCoordinates().getVelocity()).normalize(); // plane normal

        /* 8) --- Newton refinement with signed angular error ------------------ */
        final int MAX_IT = 8;
        final double ANG_TOL_DEG = 1.0e-4;                  // 0.36 arc-sec
        double Tcorr = Tph;

        for (int it = 0; it < MAX_IT; ++it) {
            AbsoluteDate guess = initialDate.shiftedBy(Tcorr);
            Vector3D rGuess = prop1.propagate(guess).getPVCoordinates().getPosition();

            double unsigned = Vector3D.angle(rB1, rGuess);           // (0,π)
            double signed = FastMath.signum(
                    Vector3D.dotProduct(Vector3D.crossProduct(rB1, rGuess), h)) * unsigned;
            double angErrDeg = FastMath.toDegrees(signed);

            if (FastMath.abs(angErrDeg) < ANG_TOL_DEG) {
                System.out.printf("Newton refine  it=%d  |Δθ|=%.6f deg  ✓%n", it, FastMath.abs(angErrDeg));
                break;
            }

            double dt = signed / nTgt;      // 1-st order (keeps the sign!)
            Tcorr -= dt;

            System.out.printf("Newton refine  it=%d  Δθ=%+.6f deg  dt=%+.1f s  → Tph=%.3f min%n",
                    it, angErrDeg, dt, Tcorr/60);
            if (it == MAX_IT-1)
                System.out.println("  (stopped -- reached max iterations)");
        }

        /* 9) --- Coast until Burn-2 ------------------------------------------ */
        final AbsoluteDate tB2 = initialDate.shiftedBy(Tcorr);
        SpacecraftState stateBeforeBurn2 = prop1.propagate(tB2);
        stateBeforeBurn2 = new SpacecraftState(stateBeforeBurn2.getOrbit(), m1);
        logState("state@beforeBurn2", stateBeforeBurn2);

        // Write first maneuver results and wait for trigger
        Map<String, String> payload1 = Utils.createDatePayload(initialDate, tB2);
        Utils.writeJsonPayload(payload1);
        Utils.logResults(resultFileName, stateAfterBurn1,
                new KeplerianOrbit(stateAfterBurn1.getOrbit()), m0, DRYMASS);
        writeManeuverTimestamp(postManeuverDateFileName, tB2);
        mqttService.sendFileAndWaitForTrigger(resultFileName, publishTopic, triggerTopic);

        // Reference orbit (if we hadn't performed any maneuvers)
        KeplerianPropagator refProp = new KeplerianPropagator(stateInit.getOrbit());
        SpacecraftState refState = refProp.propagate(tB2);
        refState = new SpacecraftState(refState.getOrbit(), m0);
        logState("Reference state (no maneuver)", refState);

        /* 10) --- Burn-2 ------------------------------------------------------ */
        final Vector3D dir2 = stateBeforeBurn2.getPVCoordinates().getVelocity()
                .normalize().scalarMultiply(dv2);
        final DateDetector secondBurnTrigger = new DateDetector(tB2.shiftedBy(0.001));
        final ImpulseManeuver secondBurn = new ImpulseManeuver(secondBurnTrigger, dir2, ISP);

        final KeplerianPropagator prop2 = new KeplerianPropagator(stateBeforeBurn2.getOrbit());
        prop2.addEventDetector(secondBurn);

        SpacecraftState finalState = prop2.propagate(tB2.shiftedBy(0.002));
        finalState = new SpacecraftState(finalState.getOrbit(), m2);
        logState("state@final", finalState);
        final Vector3D rB2 = finalState.getPVCoordinates().getPosition();

        /* 11) --- Co-location & phase check ---------------------------------- */
        double dR = rB2.subtract(rB1).getNorm();
        double dAng = FastMath.toDegrees(Vector3D.angle(rB1, rB2));
        System.out.printf("Co-location check     : Δr = %.3f m   Δθ = %.6f deg%n", dR, dAng);

        if (dAng > 1.0e-3)      // 0.001 deg  ≈ 3.6″
            throw new IllegalStateException(
                    String.format("Burn-2 not co-located with Burn-1 (Δθ = %.6f deg)", dAng));

        double phaseAchieved = MathUtils.normalizeAngle(
                meanLongitude(new KeplerianOrbit(finalState.getOrbit()))
                        - meanLongitude(new KeplerianOrbit(refState.getOrbit())), 0);
        System.out.printf("Phase achieved vs ref : %+.4f deg%n",
                FastMath.toDegrees(phaseAchieved));

        /* 12) --- Write final results ---------------------------------------- */
        Map<String, String> payload2 = Utils.createDatePayload(tB2, endHorizonDate);
        Utils.writeJsonPayload(payload2);
        Utils.logResults(resultFileName, finalState,
                new KeplerianOrbit(finalState.getOrbit()), m0, DRYMASS);
        writeManeuverTimestamp(postManeuverDateFileName, finalState.getDate());
        writeManeuverTimestamp(lastManeuverDateFile, finalState.getDate());
        mqttService.sendFileViaMQTT(resultFileName, publishTopic);

        System.out.printf("Phasing maneuver complete. Execution time: %.3f s%n",
                (System.currentTimeMillis() - wallStart)/1000.0);
    }

    @Override
    public void calculateErgolConsumption() throws IOException {
        System.out.println("Calculating propellant consumption for phase angle change of " + DELTA_THETA + "°");
        double start = System.currentTimeMillis();
        manager.addProvider(new DirectoryCrawler(orekitData));
        /* 1) --- Frames & dates ---------------------------------------------- */
        final Frame inertial = FramesFactory.getEME2000();
        final MqttService mqttService = new MqttService();
        final AbsoluteDate apsisDate = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
        final AbsoluteDate tleDate = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
        final AbsoluteDate initialDate = tleDate.shiftedBy(manoeuverRelativeDate);
        final AbsoluteDate endHorizonDate = new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());
        final double g0 = 9.80665;

        final double m0 = DRYMASS + ERGOL;
        System.out.printf("%n=== PHASING MANOEUVRE  ==================================================%n");
        System.out.printf("Dry mass              : %,.3f kg%n", DRYMASS);
        System.out.printf("Propellant on‑board   : %,.3f kg%n", ERGOL);
        System.out.printf("ISP                   : %,.1f s%n", ISP);
        System.out.printf("Initial wet mass m0   : %,.3f kg%n", m0);
        System.out.printf("Reference apsis date  : %s%n", apsisDate);
        System.out.printf("Burn‑1 epoch (t0)     : %s%n", initialDate);

        /* 2) --- Initial circular state -------------------------------------- */
        final KeplerianOrbit orbitAtApsis = new KeplerianOrbit(
                SMA, ECC, INC, PA, RAAN, ANO, PositionAngle.MEAN,
                inertial, apsisDate, MU);
        SpacecraftState stateApsis = new SpacecraftState(orbitAtApsis, m0);
        logState("state@apsis", stateApsis);

        final KeplerianPropagator prepProp = new KeplerianPropagator(orbitAtApsis);
        final SpacecraftState stateInit = new SpacecraftState(
                prepProp.propagate(initialDate).getOrbit(), m0);
        logState("state@initial", stateInit);

        /* 3) --- Phase still to perform -------------------------------------- */
        final double nTgt = FastMath.sqrt(MU / FastMath.pow(SMA, 3));
        final double drift = nTgt * initialDate.durationFrom(apsisDate);

        final double dThetaCmd = FastMath.toRadians(DELTA_THETA);
        final double dTheta = MathUtils.normalizeAngle(dThetaCmd, 0);
        System.out.printf("Δθ commanded          : %+.4f deg%n", FastMath.toDegrees(dThetaCmd));
        System.out.printf("Natural drift         : %+.4f deg%n", FastMath.toDegrees(drift));
        System.out.printf("Δθ still to perform   : %+.4f deg%n", FastMath.toDegrees(dTheta));

        /* 4) --- First-guess phasing period ---------------------------------- */
        final int kInt =  1;
        double Tph = (kInt * 2 * FastMath.PI - dTheta) / nTgt;   // s
        System.out.printf("k_int                 : %d%n", kInt);
        System.out.printf("T_ph (initial guess)  : %.3f min%n", Tph/60);

        /* 5) --- Phasing-orbit geometry -------------------------------------- */
        double aPh = FastMath.cbrt(MU * FastMath.pow(Tph/(2*FastMath.PI), 2));
        boolean outer = aPh > SMA;
        double rp = outer ? SMA : 2*aPh - SMA;
        double ra = outer ? 2*aPh - SMA : SMA;
        double ePh = (ra - rp) / (ra + rp);

        System.out.printf("Phasing orbit (first) : a=%.3f km  e=%.6f  rp=%.3f km  ra=%.3f km%n",
                aPh/1e3, ePh, rp/1e3, ra/1e3);

        /* 6) --- ΔV budget ---------------------------------------------------- */
        double vCirc = FastMath.sqrt(MU / SMA);
        double vPeriPh = FastMath.sqrt(MU*(2/SMA - 1/aPh));
        double dv1 = vPeriPh - vCirc;
        double dv2 = -dv1;
        double m1 = calculateFinalMass(m0, dv1, ISP, g0);
        double m2 = calculateFinalMass(m1, dv2, ISP, g0);
        System.out.printf("ΔV1 / ΔV2             : %+.3f / %+.3f m/s%n", dv1, dv2);
        System.out.printf("Propellant expected   : %.3f kg%n", m0-m2);

        /* 7) --- Burn-1 ------------------------------------------------------- */
        final Vector3D dir1 = stateInit.getPVCoordinates().getVelocity()
                .normalize().scalarMultiply(dv1);
        final DateDetector firstBurnTrigger = new DateDetector(initialDate.shiftedBy(0.001));
        final ImpulseManeuver firstBurn = new ImpulseManeuver(firstBurnTrigger, dir1, ISP);

        final KeplerianPropagator prop1 = new KeplerianPropagator(stateInit.getOrbit());
        prop1.addEventDetector(firstBurn);

        SpacecraftState stateAfterBurn1 = prop1.propagate(initialDate.shiftedBy(0.002));
        stateAfterBurn1 = new SpacecraftState(stateAfterBurn1.getOrbit(), m1);
        logState("state@afterBurn1", stateAfterBurn1);

        final Vector3D rB1 = stateAfterBurn1.getPVCoordinates().getPosition();
        final Vector3D h = Vector3D.crossProduct(rB1,
                stateAfterBurn1.getPVCoordinates().getVelocity()).normalize(); // plane normal

        /* 8) --- Newton refinement with signed angular error ------------------ */
        final int MAX_IT = 8;
        final double ANG_TOL_DEG = 1.0e-4;                  // 0.36 arc-sec
        double Tcorr = Tph;

        for (int it = 0; it < MAX_IT; ++it) {
            AbsoluteDate guess = initialDate.shiftedBy(Tcorr);
            Vector3D rGuess = prop1.propagate(guess).getPVCoordinates().getPosition();

            double unsigned = Vector3D.angle(rB1, rGuess);           // (0,π)
            double signed = FastMath.signum(
                    Vector3D.dotProduct(Vector3D.crossProduct(rB1, rGuess), h)) * unsigned;
            double angErrDeg = FastMath.toDegrees(signed);

            if (FastMath.abs(angErrDeg) < ANG_TOL_DEG) {
                System.out.printf("Newton refine  it=%d  |Δθ|=%.6f deg  ✓%n", it, FastMath.abs(angErrDeg));
                break;
            }

            double dt = signed / nTgt;      // 1-st order (keeps the sign!)
            Tcorr -= dt;

            System.out.printf("Newton refine  it=%d  Δθ=%+.6f deg  dt=%+.1f s  → Tph=%.3f min%n",
                    it, angErrDeg, dt, Tcorr/60);
            if (it == MAX_IT-1)
                System.out.println("  (stopped -- reached max iterations)");
        }

        /* 9) --- Coast until Burn-2 ------------------------------------------ */
        final AbsoluteDate tB2 = initialDate.shiftedBy(Tcorr);
        SpacecraftState stateBeforeBurn2 = prop1.propagate(tB2);
        stateBeforeBurn2 = new SpacecraftState(stateBeforeBurn2.getOrbit(), m1);
        logState("state@beforeBurn2", stateBeforeBurn2);

        // Write first maneuver results and wait for trigger
//        Map<String, String> payload1 = Utils.createDatePayload(initialDate, tB2);
//        Utils.writeJsonPayload(payload1);
//        Utils.logResults(resultFileName, stateAfterBurn1,
//                new KeplerianOrbit(stateAfterBurn1.getOrbit()), m0, DRYMASS);
//        writeManeuverTimestamp(postManeuverDateFileName, tB2);
//        mqttService.sendFileAndWaitForTrigger(resultFileName, publishTopic, triggerTopic);

        // Reference orbit (if we hadn't performed any maneuvers)
        KeplerianPropagator refProp = new KeplerianPropagator(stateInit.getOrbit());
        SpacecraftState refState = refProp.propagate(tB2);
        refState = new SpacecraftState(refState.getOrbit(), m0);
        logState("Reference state (no maneuver)", refState);

        /* 10) --- Burn-2 ------------------------------------------------------ */
        final Vector3D dir2 = stateBeforeBurn2.getPVCoordinates().getVelocity()
                .normalize().scalarMultiply(dv2);
        final DateDetector secondBurnTrigger = new DateDetector(tB2.shiftedBy(0.001));
        final ImpulseManeuver secondBurn = new ImpulseManeuver(secondBurnTrigger, dir2, ISP);

        final KeplerianPropagator prop2 = new KeplerianPropagator(stateBeforeBurn2.getOrbit());
        prop2.addEventDetector(secondBurn);

        SpacecraftState finalState = prop2.propagate(tB2.shiftedBy(0.002));
        finalState = new SpacecraftState(finalState.getOrbit(), m2);
        logState("state@final", finalState);
        final Vector3D rB2 = finalState.getPVCoordinates().getPosition();

        /* 11) --- Co-location & phase check ---------------------------------- */
        double dR = rB2.subtract(rB1).getNorm();
        double dAng = FastMath.toDegrees(Vector3D.angle(rB1, rB2));
        System.out.printf("Co-location check     : Δr = %.3f m   Δθ = %.6f deg%n", dR, dAng);

        if (dAng > 1.0e-3)      // 0.001 deg  ≈ 3.6″
            throw new IllegalStateException(
                    String.format("Burn-2 not co-located with Burn-1 (Δθ = %.6f deg)", dAng));

        double phaseAchieved = MathUtils.normalizeAngle(
                meanLongitude(new KeplerianOrbit(finalState.getOrbit()))
                        - meanLongitude(new KeplerianOrbit(refState.getOrbit())), 0);
        System.out.printf("Phase achieved vs ref : %+.4f deg%n",
                FastMath.toDegrees(phaseAchieved));
        double initialMass = DRYMASS + ERGOL;
        double fuelUsed = initialMass - m2;

        System.out.println("Fuel used " + fuelUsed + " kg");
        Utils.appendLineToFile(consommationErgolsFile, String.valueOf(fuelUsed));

        /* 12) --- Write final results ---------------------------------------- */
//        Map<String, String> payload2 = Utils.createDatePayload(tB2, endHorizonDate);
//        Utils.writeJsonPayload(payload2);
//        Utils.logResults(resultFileName, finalState,
//                new KeplerianOrbit(finalState.getOrbit()), m0, DRYMASS);
//        writeManeuverTimestamp(postManeuverDateFileName, finalState.getDate());
//        writeManeuverTimestamp(lastManeuverDateFile, finalState.getDate());
//        mqttService.sendFileViaMQTT(resultFileName, publishTopic);

    }

    @Override
    public void processReachOrbitTime() throws IOException {
        System.out.println("Processing orbit reach time for phase angle change of " + DELTA_THETA + "°");
        printAttributes();
        validateOrbitParameters();
        double start = System.currentTimeMillis();
        AbsoluteDate endHorizonDate = new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());

        try {
            // 1) Setup frames & dates
            Frame eme2000 = FramesFactory.getEME2000();
            AbsoluteDate apsideDate = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
            AbsoluteDate initialDate = Utils.parseDateFromTimestamp(DATE);

            // 2) Create initial orbit
            KeplerianOrbit orbitAtApside = new KeplerianOrbit(
                    SMA, ECC, INC, PA, RAAN, ANO, PositionAngle.MEAN,
                    eme2000, apsideDate, MU
            );

            // Validate parameters
            validateOrbitParameters();
            if (!isEqualOrAfterWithTolerance(initialDate, apsideDate, TIME_TOLERANCE_SECONDS)) {
                Utils.logApsideDate(apsideFile, null);
                throw new IllegalArgumentException(
                        "Initial date must be equal to or later than the apside date within the allowed tolerance.");
            }
            if (!isEqualOrAfterWithTolerance(endHorizonDate, initialDate, TIME_TOLERANCE_SECONDS)) {
                Utils.logApsideDate(apsideFile, null);
                throw new IllegalArgumentException("Initial date must be before the end of horizon time!");
            }

            // 3) Calculate phasing time
            double deltaThetaRad = FastMath.toRadians(DELTA_THETA);
            int kint = (FastMath.abs(deltaThetaRad) > FastMath.toRadians(140)) ? 2 : 1;
            double omega = FastMath.sqrt(MU / (SMA*SMA*SMA));
            double Tph = (kint * 2 * FastMath.PI - deltaThetaRad) / omega;

            // 4) Calculate the exact time when the phasing will be complete
            AbsoluteDate phasingEndDate = initialDate.shiftedBy(Tph);

            // 5) Log the time to the apsideFile
            Utils.logApsideDate(apsideFile, phasingEndDate);

            System.out.printf("Phasing period: %.2f min (%.2f s)%n", Tph/60, Tph);
            System.out.printf("Phasing completion date: %s%n", phasingEndDate);
            System.out.printf("Process completed in %.3f s%n",
                    (System.currentTimeMillis() - start) / 1000.0);
        } catch (Exception e) {
            System.err.println("Error processing orbit reach time: " + e.getMessage());
            e.printStackTrace();
            throw e;
        }
    }
}