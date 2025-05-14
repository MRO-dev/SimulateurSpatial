import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.orekit.data.DirectoryCrawler;
import org.orekit.forces.maneuvers.ImpulseManeuver;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.KeplerianOrbit;
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
 * It adjusts the spacecraft's mean anomaly by creating a desired phase angle separation.
 */
public class PhasageStrategyBackup extends AbstractManeuverStrategy {

    /** Desired mean-anomaly shift (phase angle), in radians. */
    private double MA_2;

    /**
     * Constructor for Phasage (phasing) maneuver.
     *
     * @param sma         initial SMA (meters)
     * @param ecc         initial eccentricity
     * @param incDeg      inclination (degrees)
     * @param raanDeg     RAAN (degrees)
     * @param aopDeg      argument of perigee (degrees)
     * @param meanAnomDeg initial mean anomaly (degrees)
     * @param dryMass     spacecraft dry mass
     * @param ergol       available propellant mass
     * @param isp         engine specific impulse
     * @param ma_2        desired mean anomaly shift (in degrees, we convert to radians)
     * @param dateStr     reference date (TLE-based or user-provided)
     * @param modeParameter    e.g. "blue1", "red2", etc. for your file modes
     * @param commandParameter e.g. "compute", "calculateMass", etc.
     */
    public PhasageStrategyBackup(
            double sma, double ecc, double incDeg, double raanDeg, double aopDeg, double meanAnomDeg,
            double dryMass, double ergol, double isp,
            double ma_2, // in degrees; we store as radians
            String dateStr, String modeParameter, String commandParameter
    ) throws IOException {
        super(modeParameter, commandParameter);

        // Basic orbital parameters (convert angles from degrees to radians).
        this.SMA     = sma;
        this.ECC     = ecc;
        this.INC     = FastMath.toRadians(incDeg);
        this.RAAN    = FastMath.toRadians(raanDeg);
        this.PA      = FastMath.toRadians(aopDeg);
        this.ANO     = FastMath.toRadians(meanAnomDeg);

        // Spacecraft mass parameters
        this.DRYMASS = dryMass;
        this.ERGOL   = ergol;
        this.ISP     = isp;

        // Convert the user-provided ma_2 (degrees) to radians
        this.MA_2    = FastMath.toRadians(ma_2);

        // Date from TLE or from the user
        this.DATE    = dateStr;
    }

    /** Rocket equation: final mass after a ΔV maneuver. */
    public static double calculateFinalMass(double initialMass, double deltaV, double ISP, double g0) {
        return initialMass * FastMath.exp(-Math.abs(deltaV) / (ISP * g0));
    }

    /** Check for valid orbit parameters and mass. */
    private void validateOrbitParameters() {
        if (SMA <= 0) {
            throw new IllegalArgumentException("Semi-major axis must be positive");
        }
        if (DRYMASS <= 0 || ERGOL < 0) {
            throw new IllegalArgumentException("Invalid mass parameters");
        }
    }

    /**
     * Validate that we have enough time/fuel and that the orbit is near-circular.
     */
    private void validateManeuverParameters(AbsoluteDate initialDate,
                                            AbsoluteDate apsideDate,
                                            double finalMass) {
        if (!isEqualOrAfterWithTolerance(initialDate, apsideDate, TIME_TOLERANCE_SECONDS)) {
            throw new IllegalArgumentException("Initial date must be >= reference apside date");
        }
        if (finalMass < DRYMASS || ERGOL < 1e-3) {
            throw new IllegalArgumentException("Not enough propellant for phasing maneuver");
        }
        // For a simpler phasing calculation, we assume near-circular orbits
        if (ECC > 0.05) {
            throw new IllegalArgumentException("Orbit eccentricity too high for accurate phasing calculation");
        }
    }

    /**
     * Load mass data from files. If isMassCalculation == true, we also parse "MA_2" from cmdData.
     */
    @Override
    public void loadMassData(Boolean isMassCalculation) throws IOException {
        super.loadMassData(isMassCalculation);
        if (isMassCalculation && cmdData.containsKey("MA_2")) {
            this.MA_2 = FastMath.toRadians(Double.parseDouble(cmdData.get("MA_2")));
        }
    }

    /**
     * Load time data from files. If isTimeCalculation == true, we also parse "MA_2" from tipData.
     */
    @Override
    public void loadTimeData(Boolean isTimeCalculation) throws IOException {
        super.loadTimeData(isTimeCalculation);
        if (tipData.containsKey("MA_2")) {
            this.MA_2 = FastMath.toRadians(Double.parseDouble(tipData.get("MA_2")));
        }
        System.out.println(" - timeIntermediateParametersFile => ManeuverType: " + MANEUV_TYPE
                + ", PhaseAngle: " + FastMath.toDegrees(MA_2) + "°");
    }

    /**
     * Compute and execute the phasing maneuver:
     *   1) Insert into lower or higher phasing orbit.
     *   2) Wait one (or possibly multiple) phasing orbit periods.
     *   3) Return to the original semi-major axis, offset by target mean anomaly.
     */
    @Override
    public void computeAndExecute() throws IOException, ParseException {

        double startTime = System.currentTimeMillis();
        manager.addProvider(new DirectoryCrawler(orekitData));

        // Setup frames, times
        Frame eme2000 = FramesFactory.getEME2000();
        MqttService mqttService = new MqttService();  // re-enable if you want real MQTT usage

        AbsoluteDate dateTLE       = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
        AbsoluteDate initialDate   = dateTLE.shiftedBy(manoeuverRelativeDate);
        AbsoluteDate apsideDate    = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
        AbsoluteDate endHorizonDate= new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());

        // Build orbit at apside date, then propagate to initial date
        KeplerianOrbit orbitAtApside = new KeplerianOrbit(
                SMA, ECC, INC, PA, RAAN, ANO,
                PositionAngle.MEAN, eme2000, apsideDate, MU);

        double initialMass = DRYMASS + ERGOL;
        validateOrbitParameters();

        KeplerianPropagator initialProp = new KeplerianPropagator(orbitAtApside);
        SpacecraftState stateAtInitial = initialProp.propagate(initialDate);

        KeplerianOrbit initialOrbit = new KeplerianOrbit(stateAtInitial.getOrbit());
        SpacecraftState initialState = new SpacecraftState(initialOrbit, initialMass);

        double T0 = initialOrbit.getKeplerianPeriod();

        // Normalize the phase angle to [-π, +π]
        double deltaTheta = FastMath.atan2(FastMath.sin(MA_2), FastMath.cos(MA_2));
        double fraction   = deltaTheta / (2 * FastMath.PI);

        // Basic single-orbit approach to T_phase
        double T_phase;
        if (deltaTheta >= 0) {
            // "catch up" => shorter period
            T_phase = T0 * (1.0 - fraction);
        } else {
            // "lag behind" => longer period
            T_phase = T0 * (1.0 + FastMath.abs(fraction));
        }
        if (T_phase <= 0) {
            throw new IllegalArgumentException("Non-physical phasing orbit period (check phase angle).");
        }

        // Compute phasing orbit's semi-major axis
        double a_phase = FastMath.cbrt(MU * (T_phase / (2.0 * FastMath.PI))
                * (T_phase / (2.0 * FastMath.PI)));

        boolean lowerOrbit = (T_phase < T0);
        double r0 = initialOrbit.getA();
        double r_other = 2.0 * a_phase - r0;  // perigee if lowerOrbit

        // -----------------------------------------------------------------
        //   ADDED LOOP: if perigee is below Earth, add extra orbits
        // -----------------------------------------------------------------
        int extraOrbits = 0;
        while (lowerOrbit && r_other < 6371000.0) {
            extraOrbits++;
            System.out.println("Phasing orbit perigee below Earth! Adding another full orbit (k="
                    + extraOrbits + ") to T_phase.");

            // We add 1 full original period to T_phase
            T_phase += T0;

            // Recompute a_phase, check again
            a_phase = FastMath.cbrt(MU * (T_phase / (2.0 * FastMath.PI))
                    * (T_phase / (2.0 * FastMath.PI)));
            r_other = 2.0 * a_phase - r0;

            // If T_phase is now bigger than T0, we might not be "lowerOrbit" anymore
            lowerOrbit = (T_phase < (T0 * (extraOrbits + 1.0)));
        }
        // If STILL too low after multiple orbits, we give up
        if (lowerOrbit && r_other < 6371000.0) {
            throw new IllegalArgumentException("Even after adding extra orbits, perigee is below Earth surface!");
        }

        // Now recalc velocities
        double v_circ          = FastMath.sqrt(MU / r0);
        double v_phasing_at_r0 = FastMath.sqrt(MU * (2.0 / r0 - 1.0 / a_phase));
        double deltaV1         = v_phasing_at_r0 - v_circ;
        double deltaV2         = v_circ - v_phasing_at_r0;

        // Check final mass after both burns
        double massAfter1 = calculateFinalMass(initialMass, deltaV1, ISP, g0);
        double finalMass  = calculateFinalMass(massAfter1, deltaV2, ISP, g0);

        // Validate feasibility: time, propellant, near-circular orbit
        validateManeuverParameters(initialDate, apsideDate, finalMass);

        // ================= 1) First Burn: Enter phasing orbit =================
        Vector3D burnDirection1 = initialState.getPVCoordinates()
                .getVelocity()
                .normalize()
                .scalarMultiply(deltaV1);

        DateDetector firstBurnTrigger = new DateDetector(initialDate.shiftedBy(0.001));
        ImpulseManeuver firstBurn     = new ImpulseManeuver(firstBurnTrigger, burnDirection1, ISP);

        KeplerianPropagator propagator = new KeplerianPropagator(initialOrbit);
        propagator.addEventDetector(firstBurn);
        SpacecraftState stateAfterFirst = propagator.propagate(firstBurnTrigger.getDate());

        // Jump forward by T_phase
        AbsoluteDate firstBurnEndDate = stateAfterFirst.getDate().shiftedBy(T_phase);
        stateAfterFirst = new SpacecraftState(stateAfterFirst.getOrbit(), massAfter1);

        KeplerianOrbit orbitAfterFirst = new KeplerianOrbit(stateAfterFirst.getOrbit());

        // Log and wait for trigger
        Map<String, String> payload1 = Utils.createDatePayload(initialDate, firstBurnEndDate);
        Utils.writeJsonPayload(payload1);
        Utils.logResults(resultFileName, stateAfterFirst, orbitAfterFirst, initialMass, DRYMASS);

        writeManeuverTimestamp(postManeuverDateFileName, firstBurnEndDate);
        //mqttService.sendFileAndWaitForTrigger(resultFileName, publishTopic, triggerTopic);
        // ^ re-enable if you want to block/wait for user signal

        // ============== 2) Second Burn: Return to original orbit ===============
        propagator = new KeplerianPropagator(stateAfterFirst.getOrbit());
        SpacecraftState stateAtReturn = propagator.propagate(firstBurnEndDate);
        stateAtReturn = new SpacecraftState(stateAtReturn.getOrbit(), massAfter1);

        Vector3D burnDirection2 = stateAtReturn.getPVCoordinates()
                .getVelocity()
                .normalize()
                .scalarMultiply(deltaV2);

        DateDetector secondBurnTrigger = new DateDetector(firstBurnEndDate.shiftedBy(0.001));
        ImpulseManeuver secondBurn     = new ImpulseManeuver(secondBurnTrigger, burnDirection2, ISP);

        propagator = new KeplerianPropagator(stateAtReturn.getOrbit());
        propagator.addEventDetector(secondBurn);
        SpacecraftState finalStateLocal = propagator.propagate(secondBurnTrigger.getDate());
        finalStateLocal = new SpacecraftState(finalStateLocal.getOrbit(), finalMass);

        KeplerianOrbit finalOrbit = new KeplerianOrbit(finalStateLocal.getOrbit());

        // Log final results
        Map<String, String> payload2 = Utils.createDatePayload(firstBurnEndDate, endHorizonDate);
        Utils.writeJsonPayload(payload2);
        Utils.logResults(resultFileName, finalStateLocal, finalOrbit, initialMass, DRYMASS);

        writeManeuverTimestamp(postManeuverDateFileName, finalStateLocal.getDate());
        writeManeuverTimestamp(lastManeuverDateFile, finalStateLocal.getDate());
        //mqttService.sendFileViaMQTT(resultFileName, publishTopic);

        double elapsedSec = (System.currentTimeMillis() - startTime) / 1000.0;
        System.out.println("Phasing maneuver complete. Execution time: " + elapsedSec + " s");
    }

    /**
     * Calculate total ergol consumption for the phasing maneuver without
     * actually executing the burn/propagation steps.
     */
    @Override
    public void calculateErgolConsumption() throws IOException {
        double initialMass = DRYMASS + ERGOL;
        double deltaTheta  = MA_2;

        System.out.println("Initial mass: " + initialMass + " kg");
        System.out.println("Desired phase angle: " + FastMath.toDegrees(deltaTheta) + "°");

        // Original orbit period
        double T0 = 2.0 * FastMath.PI * FastMath.sqrt(Math.pow(SMA, 3) / MU);

        // Normalize to [-π, +π]
        deltaTheta = FastMath.atan2(FastMath.sin(deltaTheta), FastMath.cos(deltaTheta));
        double fraction = deltaTheta / (2.0 * FastMath.PI);

        double T_phase;
        if (deltaTheta >= 0) {
            T_phase = T0 * (1.0 - fraction);
        } else {
            T_phase = T0 * (1.0 + FastMath.abs(fraction));
        }

        // We might do the same 'multi-orbit' logic here, if you want to approximate consumption
        // with extra orbits. But for now, the code is the simpler single-orbit approach.

        double a_phase = FastMath.cbrt(MU * (T_phase / (2.0 * FastMath.PI))
                * (T_phase / (2.0 * FastMath.PI)));
        double r0     = SMA;
        double v_circ = FastMath.sqrt(MU / r0);
        double v_phase_at_r0 = FastMath.sqrt(MU * (2.0 / r0 - 1.0 / a_phase));

        double deltaV1 = v_phase_at_r0 - v_circ;
        double deltaV2 = v_circ - v_phase_at_r0;

        System.out.println("ΔV1 (enter phasing): " + deltaV1 + " m/s");
        System.out.println("ΔV2 (exit phasing):  " + deltaV2 + " m/s");

        double massAfter1 = calculateFinalMass(initialMass, deltaV1, ISP, g0);
        double massAfter2 = calculateFinalMass(massAfter1, deltaV2, ISP, g0);

        double fuel1 = (initialMass - massAfter1);
        double fuel2 = (massAfter1 - massAfter2);
        double totalFuel = (initialMass - massAfter2);

        System.out.println("Fuel for first burn:  " + fuel1 + " kg");
        System.out.println("Fuel for second burn: " + fuel2 + " kg");
        System.out.println("Total fuel consumed:  " + totalFuel + " kg");
    }

    /**
     * For time-based calculations only: figure out when the second burn
     * happens based on the desired phase angle shift (single-orbit approach).
     */
    @Override
    public void processReachOrbitTime() throws IOException {
        printAttributes();
        AbsoluteDate initialDate = Utils.parseDateFromTimestamp(DATE);
        System.out.println("Initial date: " + initialDate);

        double T0 = 2.0 * FastMath.PI * FastMath.sqrt(Math.pow(SMA, 3) / MU);
        double deltaTheta = FastMath.atan2(FastMath.sin(MA_2), FastMath.cos(MA_2));
        double fraction = deltaTheta / (2.0 * FastMath.PI);

        double T_phase = (deltaTheta >= 0) ? (T0 * (1.0 - fraction))
                : (T0 * (1.0 + FastMath.abs(fraction)));

        AbsoluteDate secondBurnDate = initialDate.shiftedBy(T_phase);
        Utils.logApsideDate(apsideFile, secondBurnDate);

        System.out.println("Logged return-to-orbit burn time: " + secondBurnDate);
    }
}
