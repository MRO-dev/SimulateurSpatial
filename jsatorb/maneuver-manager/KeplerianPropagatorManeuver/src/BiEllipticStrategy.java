import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
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

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.util.Map;

/**
 * This class implements a bi-elliptic transfer maneuver strategy.
 * It performs three impulses to transfer from SMA1 to SMA2 via a high apogee orbit (SMA3 = 2*SMA2).
 */
public class BiEllipticStrategy extends AbstractManeuverStrategy {
    private double SMA_2;  // target semi-major axis (m)
    private double SMA_3;  // intermediate high-apogee semi-major axis (m)

    // Constructor
    public BiEllipticStrategy(
            double sma1, double ecc1, double incDeg, double raanDeg, double aopDeg, double meanAnomDeg,
            double dryMass, double ergol, double isp,
            double sma2,
            String dateStr, String modeParameter, String commandParameter) throws IOException {
        super(modeParameter, commandParameter);
        this.SMA    = sma1;
        this.ECC    = ecc1;
        this.INC    = FastMath.toRadians(incDeg);
        this.RAAN   = FastMath.toRadians(raanDeg);
        this.PA     = FastMath.toRadians(aopDeg);
        this.ANO    = FastMath.toRadians(meanAnomDeg);
        this.DRYMASS = dryMass;
        this.ERGOL   = ergol;
        this.ISP     = isp;
        this.SMA_2   = sma2;
        this.DATE    = dateStr;
        // Set intermediate apogee orbit SMA (twice the target SMA)
        this.SMA_3   = 2 * sma2;
    }

    /** Calculate final mass after a ΔV using the Tsiolkovsky rocket equation. */
    public static double calculateFinalMass(double initialMass, double deltaV, double ISP, double g0) {
        return initialMass * FastMath.exp(-Math.abs(deltaV) / (ISP * g0));
    }

    /** Validate orbit and mass parameters (SMA, masses) before the maneuver. */
    private void validateOrbitParameters() {
        if (SMA <= 0 || SMA_2 <= 0) {
            throw new IllegalArgumentException("Semi-major axes must be positive");
        }
        if (DRYMASS <= 0 || ERGOL < 0) {
            throw new IllegalArgumentException("Invalid mass parameters");
        }
    }

    /** Validate maneuver feasibility (time order and propellant). */
    private void validateManeuverParameters(AbsoluteDate initialDate, AbsoluteDate apsideDate, double finalMass) {
        if (!isEqualOrAfterWithTolerance(initialDate, apsideDate, TIME_TOLERANCE_SECONDS)) {
            throw new IllegalArgumentException("Initial date must be on/after reference apside date");
        }
        if (finalMass < DRYMASS || ERGOL < 1e-3) {
            throw new IllegalArgumentException("Insufficient propellant for bi-elliptic transfer");
        }
        if (ECC > 0.05) {
            throw new IllegalArgumentException("Eccentricity too high for bi-elliptic transfer assumption");
        }
    }

    @Override
    public void loadMassData(Boolean isMassCalculation) throws IOException {
        super.loadMassData(isMassCalculation);
        if (isMassCalculation) {
            // Load target SMA (km) from command data and compute SMA_3
            if (cmdData.containsKey("SMA_2")) {
                this.SMA_2 = Double.parseDouble(cmdData.get("SMA_2")) * 1000.0;
                this.SMA_3 = 2 * this.SMA_2;
            }
        }
    }

    @Override
    public void loadTimeData(Boolean isTimeCalculation) throws IOException {
        super.loadTimeData(isTimeCalculation);
        // Load target SMA (km) from intermediate parameters and set SMA_3
        if (tipData.containsKey("SMA_2")) {
            this.SMA_2 = Double.parseDouble(tipData.get("SMA_2")) * 1000.0;
            this.SMA_3 = 2 * this.SMA_2;
        }
        System.out.println(" - timeIntermediateParametersFile => ManeuverType: " + MANEUV_TYPE
                + ", SMA_2: " + SMA_2 + " m, SMA_3: " + SMA_3 + " m");
    }

    @Override
    public void computeAndExecute() throws IOException, java.text.ParseException {
        double startTime = System.currentTimeMillis();
        manager.addProvider(new DirectoryCrawler(orekitData));

        // Setup frames and dates
        Frame eme2000 = FramesFactory.getEME2000();
        MqttService mqttService = new MqttService();
        AbsoluteDate dateTLE = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
        AbsoluteDate initialDate = dateTLE.shiftedBy(manoeuverRelativeDate);
        AbsoluteDate apsideDate  = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
        AbsoluteDate endHorizonDate = new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());

        // Initialize orbit at apside and propagate to initial date
        KeplerianOrbit orbitAtApside = new KeplerianOrbit(
                SMA, ECC, INC, PA, RAAN, ANO, PositionAngle.MEAN, eme2000, apsideDate, MU);
        double initialMass = DRYMASS + ERGOL;
        KeplerianPropagator propagator = new KeplerianPropagator(orbitAtApside);
        SpacecraftState stateAtInitial = propagator.propagate(initialDate);
        KeplerianOrbit initialOrbit = new KeplerianOrbit(stateAtInitial.getOrbit());
        SpacecraftState initialState = new SpacecraftState(initialOrbit, initialMass);

        // Calculate ΔV for each of the three impulses
        double r1 = initialOrbit.getA();
        double r2 = SMA_2;
        double r3 = SMA_3;
        // Velocities on initial and final circular orbits
        double v_circ1 = FastMath.sqrt(MU / r1);
        double v_circ2 = FastMath.sqrt(MU / r2);
        // First transfer orbit (perigee=r1, apogee=r3)
        double a_trans1 = (r1 + r3) / 2.0;
        double v_trans1_perigee = FastMath.sqrt(MU * (2/r1 - 1/a_trans1));
        double v_trans1_apogee  = FastMath.sqrt(MU * (2/r3 - 1/a_trans1));
        // Second transfer orbit (perigee=r2, apogee=r3)
        double a_trans2 = (r2 + r3) / 2.0;
        double v_trans2_perigee = FastMath.sqrt(MU * (2/r2 - 1/a_trans2));
        double v_trans2_apogee  = FastMath.sqrt(MU * (2/r3 - 1/a_trans2));
        // ΔV calculations
        double deltaV1 = v_trans1_perigee - v_circ1;      // burn at r1 to raise apogee to r3
        double deltaV2 = v_trans2_apogee - v_trans1_apogee; // burn at r3 to drop perigee to r2
        double deltaV3 = v_circ2 - v_trans2_perigee;      // burn at r2 to circularize

        // Compute mass after each maneuver
        double massAfter1 = calculateFinalMass(initialMass, deltaV1, ISP, g0);
        double massAfter2 = calculateFinalMass(massAfter1, deltaV2, ISP, g0);
        double finalMass  = calculateFinalMass(massAfter2, deltaV3, ISP, g0);
        validateManeuverParameters(initialDate, apsideDate, finalMass);

        // ---------- First Maneuver (at initial orbit) ----------
        // Impulse at initialDate to raise apogee to r3
        Vector3D burnDirection1 = initialState.getPVCoordinates().getVelocity().normalize().scalarMultiply(deltaV1);
        DateDetector firstBurnTrigger = new DateDetector(initialDate.shiftedBy(0.001));
        ImpulseManeuver firstBurn = new ImpulseManeuver(firstBurnTrigger, burnDirection1, ISP);

        propagator = new KeplerianPropagator(initialOrbit);
        propagator.addEventDetector(firstBurn);
        SpacecraftState stateAfterFirst = propagator.propagate(firstBurnTrigger.getDate());
        AbsoluteDate firstBurnEndDate = stateAfterFirst.getDate().shiftedBy(
                Math.PI * FastMath.sqrt(Math.pow(a_trans1, 3) / MU)  // time to reach apogee r3 (half period of first transfer)
        );
        // Update mass after first burn
        stateAfterFirst = new SpacecraftState(stateAfterFirst.getOrbit(), massAfter1);
        KeplerianOrbit orbitAfterFirst = new KeplerianOrbit(stateAfterFirst.getOrbit());

        // Log after first maneuver
        Map<String, String> payload1 = Utils.createDatePayload(initialDate, firstBurnEndDate);
        Utils.writeJsonPayload(payload1);
        Utils.logResults(resultFileName, stateAfterFirst, orbitAfterFirst, initialMass, DRYMASS);
        writeManeuverTimestamp(postManeuverDateFileName, firstBurnEndDate);
        mqttService.sendFileAndWaitForTrigger(resultFileName, publishTopic, triggerTopic);

        // ---------- Second Maneuver (at apogee r3) ----------
        // Propagate to apogee time for second burn
        propagator = new KeplerianPropagator(stateAfterFirst.getOrbit());
        SpacecraftState stateAtApogee = propagator.propagate(firstBurnEndDate);
        // Impulse at apogee to lower perigee to r2
        Vector3D burnDirection2 = stateAtApogee.getPVCoordinates().getVelocity().normalize().scalarMultiply(deltaV2);
        DateDetector secondBurnTrigger = new DateDetector(firstBurnEndDate.shiftedBy(0.001));
        ImpulseManeuver secondBurn = new ImpulseManeuver(secondBurnTrigger, burnDirection2, ISP);

        propagator = new KeplerianPropagator(stateAtApogee.getOrbit());
        propagator.addEventDetector(secondBurn);
        SpacecraftState stateAfterSecond = propagator.propagate(secondBurnTrigger.getDate());
        AbsoluteDate secondBurnEndDate = stateAfterSecond.getDate().shiftedBy(
                Math.PI * FastMath.sqrt(Math.pow(a_trans2, 3) / MU)  // time to reach perigee r2 (half period of second transfer)
        );
        // Update mass after second burn
        stateAfterSecond = new SpacecraftState(stateAfterSecond.getOrbit(), massAfter2);
        KeplerianOrbit orbitAfterSecond = new KeplerianOrbit(stateAfterSecond.getOrbit());

        // Log after second maneuver
        Map<String, String> payload2 = Utils.createDatePayload(firstBurnEndDate, secondBurnEndDate);
        Utils.writeJsonPayload(payload2);
        Utils.logResults(resultFileName, stateAfterSecond, orbitAfterSecond, initialMass, DRYMASS);
        writeManeuverTimestamp(postManeuverDateFileName, secondBurnEndDate);
        mqttService.sendFileAndWaitForTrigger(resultFileName, publishTopic, triggerTopic);

        // ---------- Third Maneuver (at perigee r2) ----------
        // Propagate to perigee time for third burn
        propagator = new KeplerianPropagator(stateAfterSecond.getOrbit());
        SpacecraftState stateAtPerigee = propagator.propagate(secondBurnEndDate);
        // Impulse at perigee to circularize at r2
        Vector3D burnDirection3 = stateAtPerigee.getPVCoordinates().getVelocity().normalize().scalarMultiply(deltaV3);
        DateDetector thirdBurnTrigger = new DateDetector(secondBurnEndDate.shiftedBy(0.001));
        ImpulseManeuver thirdBurn = new ImpulseManeuver(thirdBurnTrigger, burnDirection3, ISP);

        propagator = new KeplerianPropagator(stateAtPerigee.getOrbit());
        propagator.addEventDetector(thirdBurn);
        SpacecraftState finalState = propagator.propagate(thirdBurnTrigger.getDate());
        finalState = new SpacecraftState(finalState.getOrbit(), finalMass);  // update mass after third burn
        KeplerianOrbit finalOrbit = new KeplerianOrbit(finalState.getOrbit());

        // Log after third (final) maneuver
        Map<String, String> payload3 = Utils.createDatePayload(secondBurnEndDate, endHorizonDate);
        Utils.writeJsonPayload(payload3);
        Utils.logResults(resultFileName, finalState, finalOrbit, initialMass, DRYMASS);
        writeManeuverTimestamp(postManeuverDateFileName, finalState.getDate());
        writeManeuverTimestamp(lastManeuverDateFile, finalState.getDate());
        mqttService.sendFileViaMQTT(resultFileName, publishTopic);

        System.out.println("Bi-elliptic transfer complete. Execution time: "
                + (System.currentTimeMillis() - startTime) / 1000.0 + " s");
    }

    @Override
    public void calculateErgolConsumption() throws IOException {
        double initialMass = DRYMASS + ERGOL;
        System.out.println("Initial mass: " + initialMass + " kg");
        // Recompute ΔVs for transparency
        double r1 = SMA;
        double r2 = SMA_2;
        double r3 = 2 * r2;
        double v_circ1 = FastMath.sqrt(MU / r1);
        double v_circ2 = FastMath.sqrt(MU / r2);
        double a_trans1 = (r1 + r3) / 2.0;
        double a_trans2 = (r2 + r3) / 2.0;
        double v_trans1_perigee = FastMath.sqrt(MU * (2/r1 - 1/a_trans1));
        double v_trans1_apogee  = FastMath.sqrt(MU * (2/r3 - 1/a_trans1));
        double v_trans2_perigee = FastMath.sqrt(MU * (2/r2 - 1/a_trans2));
        double v_trans2_apogee  = FastMath.sqrt(MU * (2/r3 - 1/a_trans2));
        double deltaV1 = v_trans1_perigee - v_circ1;
        double deltaV2 = v_trans2_apogee - v_trans1_apogee;
        double deltaV3 = v_circ2 - v_trans2_perigee;
        System.out.println("ΔV1: " + deltaV1 + " m/s, ΔV2: " + deltaV2 + " m/s, ΔV3: " + deltaV3 + " m/s");
        double mAfter1 = calculateFinalMass(initialMass, deltaV1, ISP, g0);
        double mAfter2 = calculateFinalMass(mAfter1, deltaV2, ISP, g0);
        double mAfter3 = calculateFinalMass(mAfter2, deltaV3, ISP, g0);
        System.out.println("Fuel for burn 1: " + (initialMass - mAfter1) + " kg");
        System.out.println("Fuel for burn 2: " + (mAfter1 - mAfter2) + " kg");
        System.out.println("Fuel for burn 3: " + (mAfter2 - mAfter3) + " kg");
        System.out.println("Total fuel consumed: " + (initialMass - mAfter3) + " kg");
    }

    @Override
    public void processReachOrbitTime() throws IOException {
        printAttributes();
        AbsoluteDate initialDate = Utils.parseDateFromTimestamp(DATE);
        System.out.println("Initial date: " + initialDate);
        // Compute time to first apogee (for second burn scheduling)
        double a_trans1 = (SMA + SMA_3) / 2.0;
        double timeToApogee = Math.PI * Math.sqrt(Math.pow(a_trans1, 3) / MU);
        AbsoluteDate apogeeDate = initialDate.shiftedBy(timeToApogee);
        Utils.logApsideDate(apsideFile, apogeeDate);
        System.out.println("Logged high-apogee reach time: " + apogeeDate);
    }
}
