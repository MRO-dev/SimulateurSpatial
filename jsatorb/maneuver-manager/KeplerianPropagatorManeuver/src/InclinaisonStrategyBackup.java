/*
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

*/
/**
 * This class implements an orbital inclination-change maneuver strategy.
 * It adjusts the orbit's inclination to a target value without altering the semi-major axis.
 * It uses vector-based rotation for the burn and includes implementations for calculating ergol consumption
 * and determining the "reach orbit time" (when the new orbit is fully established).
 *//*

public class InclinaisonStrategy extends AbstractManeuverStrategy {

    // Target inclination in radians (set from user input INC_2).
    private double INC_2;

    // Store the final spacecraft state after maneuver for later use.
    private SpacecraftState finalState;

    // ==============================
    //          Constructor
    // ==============================
    public InclinaisonStrategy(
            double sma,
            double ecc,
            double incDeg,
            double raanDeg,
            double aopDeg,
            double meanAnomDeg,
            double dryMass,
            double ergol,
            double isp,
            double inc_2,
            String dateStr,
            String modeParameter,
            String commandParameter
    ) throws IOException {
        super(modeParameter, commandParameter);
        this.SMA     = sma;
        this.ECC     = ecc;
        this.INC     = FastMath.toRadians(incDeg);
        this.RAAN    = FastMath.toRadians(raanDeg);
        this.PA      = FastMath.toRadians(aopDeg);
        this.ANO     = FastMath.toRadians(meanAnomDeg);
        this.DRYMASS = dryMass;
        this.ERGOL   = ergol;
        this.ISP     = isp;
        this.INC_2   = FastMath.toRadians(inc_2);
        this.DATE    = dateStr;
    }

    */
/**
     * Calculate final mass after an impulse maneuver using the rocket equation.
     * m_final = m_initial * exp(-|Δv|/(ISP * g0))
     *//*

    public static double calculateFinalMass(double initialMass, double deltaV, double ISP, double g0) {
        return initialMass * FastMath.exp(-Math.abs(deltaV) / (ISP * g0));
    }

    */
/** Validate basic orbital and mass parameters. *//*

    private void validateOrbitParameters() {
        if (SMA <= 0) {
            throw new IllegalArgumentException("Semi-major axis must be positive");
        }
        if (DRYMASS <= 0 || ERGOL < 0) {
            throw new IllegalArgumentException("Invalid mass parameters");
        }
    }

    */
/**
     * Validate that the maneuver is feasible (time and fuel constraints).
     *//*

    private void validateManeuverParameters(AbsoluteDate initialDate, AbsoluteDate apsideDate, double finalMass) {
        if (!isEqualOrAfterWithTolerance(initialDate, apsideDate, TIME_TOLERANCE_SECONDS)) {
            throw new IllegalArgumentException("Initial date must be equal or later than the reference apside date within tolerance");
        }
        double fuelThreshold = 1e-3;
        if (finalMass < DRYMASS || ERGOL < fuelThreshold) {
            throw new IllegalArgumentException("Insufficient propellant for inclination change");
        }
        if (ECC > 0.05) {
            throw new IllegalArgumentException("Eccentricity too high for near-circular plane change assumptions");
        }
    }

    @Override
    public void loadMassData(Boolean isMassCalculation) throws IOException {
        super.loadMassData(isMassCalculation);
        if (isMassCalculation && cmdData.containsKey("INC_2")) {
            this.INC_2 = FastMath.toRadians(Double.parseDouble(cmdData.get("INC_2")));
        }
    }

    @Override
    public void loadTimeData(Boolean isTimeCalculation) throws IOException {
        super.loadTimeData(isTimeCalculation);
        if (tipData.containsKey("INC_2")) {
            this.INC_2 = FastMath.toRadians(Double.parseDouble(tipData.get("INC_2")));
        }
        System.out.println(" - timeIntermediateParametersFile => ManeuverType: " + MANEUV_TYPE
                + ", Target Inclination: " + FastMath.toDegrees(INC_2) + "°");
    }

    */
/**
     * Find the next node crossing (ascending or descending) after the current state.
     *//*

    private AbsoluteDate findNextNodeDate(SpacecraftState state) {
        KeplerianOrbit kep = new KeplerianOrbit(state.getOrbit());
        double trueAnomaly = kep.getTrueAnomaly();
        double raan = kep.getRightAscensionOfAscendingNode();
        double aop  = kep.getPerigeeArgument();

        // Ascending node: TA = -(RAAN + AoP) mod 2π
        double taAsc = -(raan + aop) % (2 * FastMath.PI);
        if (taAsc < 0) {
            taAsc += 2 * FastMath.PI;
        }
        // Descending node: TA = (taAsc + π) mod 2π
        double taDesc = (taAsc + Math.PI) % (2 * FastMath.PI);

        double period = kep.getKeplerianPeriod();
        double diffAsc = (taAsc - trueAnomaly + 2 * FastMath.PI) % (2 * FastMath.PI);
        double diffDesc = (taDesc - trueAnomaly + 2 * FastMath.PI) % (2 * FastMath.PI);
        double timeToAsc = diffAsc / (2 * FastMath.PI) * period;
        double timeToDesc = diffDesc / (2 * FastMath.PI) * period;

        if (timeToAsc <= timeToDesc) {
            System.out.println("Next node is ascending, in " + timeToAsc + " s");
            return state.getDate().shiftedBy(timeToAsc);
        } else {
            System.out.println("Next node is descending, in " + timeToDesc + " s");
            return state.getDate().shiftedBy(timeToDesc);
        }
    }

    */

import org.hipparchus.util.FastMath;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.time.AbsoluteDate;

import java.io.IOException;

/**
     * Compute the exact Δv vector (impulse) to change the orbital inclination from the current
     * value to the target (INC_2) while preserving the speed.
     *//*

    private Vector3D computeInclinationChangeImpulse(SpacecraftState state, double incTarget) {
        KeplerianOrbit ko = new KeplerianOrbit(state.getOrbit());
        double incCurrent = ko.getI();
        double deltaInc = incTarget - incCurrent;
        if (FastMath.abs(deltaInc) < 1e-9) {
            return Vector3D.ZERO;
        }
        Vector3D v0 = state.getPVCoordinates().getVelocity();
        double speed = v0.getNorm();
        if (speed < 1e-6) {
            return Vector3D.ZERO;
        }
        // Rotation axis: (orbitNormal x v0)
        Vector3D orbitNormal = state.getPVCoordinates().getMomentum().normalize();
        Vector3D rotationAxis = orbitNormal.crossProduct(v0).normalize();
        org.hipparchus.geometry.euclidean.threed.Rotation rot =
                new org.hipparchus.geometry.euclidean.threed.Rotation(rotationAxis, deltaInc,
                        org.hipparchus.geometry.euclidean.threed.RotationConvention.VECTOR_OPERATOR);
        Vector3D v1 = rot.applyTo(v0);
        double newSpeed = v1.getNorm();
        if (FastMath.abs(newSpeed - speed) > 1e-7) {
            v1 = v1.scalarMultiply(speed / newSpeed);
        }
        return v1.subtract(v0);
    }

    @Override
    public void computeAndExecute() throws IOException, ParseException {
        double startTime = System.currentTimeMillis();

        // Load Orekit data
        manager.addProvider(new DirectoryCrawler(orekitData));
        validateOrbitParameters();

        Frame eme2000 = FramesFactory.getEME2000();
        MqttService mqttService = new MqttService();

        AbsoluteDate dateTLE = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
        AbsoluteDate initialDate = dateTLE.shiftedBy(manoeuverRelativeDate);
        AbsoluteDate apsideDate = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
        AbsoluteDate endHorizonDate = new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());

        // Create initial orbit at reference apside date.
        KeplerianOrbit orbitAtApside = new KeplerianOrbit(
                SMA, ECC, INC, PA, RAAN, ANO, PositionAngle.MEAN,
                eme2000, apsideDate, MU
        );
        double initialMassGlobal = DRYMASS + ERGOL;
        KeplerianPropagator propagator = new KeplerianPropagator(orbitAtApside);

        // Propagate to the maneuver start date.
        SpacecraftState stateAtInitialDate = propagator.propagate(initialDate);
        KeplerianOrbit initialOrbit = new KeplerianOrbit(stateAtInitialDate.getOrbit());
        SpacecraftState initialState = new SpacecraftState(initialOrbit, initialMassGlobal);

        // Validate maneuver feasibility.
        double estFinalMass = initialMassGlobal;
        validateManeuverParameters(initialDate, apsideDate, estFinalMass);

        // Print inclination info.
        double currentInc = initialOrbit.getI();
        double deltaInc = INC_2 - currentInc;
        System.out.println("Current inclination: " + FastMath.toDegrees(currentInc) + "°");
        System.out.println("Target inclination:  " + FastMath.toDegrees(INC_2) + "°");
        System.out.println("Delta inclination:   " + FastMath.toDegrees(deltaInc) + "°");

        // Decide which node to use.
        boolean incIncreasing = deltaInc > 0;
        double trueAnomaly = initialOrbit.getTrueAnomaly();
        double raan = initialOrbit.getRightAscensionOfAscendingNode();
        double aop = initialOrbit.getPerigeeArgument();
        double taAscending = -(raan + aop) % (2 * Math.PI);
        if (taAscending < 0) { taAscending += 2 * Math.PI; }
        double taDescending = (taAscending + Math.PI) % (2 * Math.PI);
        double orbitPeriod = initialOrbit.getKeplerianPeriod();
        double diffAsc = (taAscending - trueAnomaly + 2 * Math.PI) % (2 * Math.PI);
        double diffDesc = (taDescending - trueAnomaly + 2 * Math.PI) % (2 * Math.PI);
        double timeToNode;
        if (incIncreasing) {
            timeToNode = diffAsc / (2 * Math.PI) * orbitPeriod;
            if (timeToNode < TIME_TOLERANCE_SECONDS) { timeToNode = orbitPeriod; }
            System.out.println("Selected ascending node for first burn, in " + timeToNode + " s");
        } else {
            timeToNode = diffDesc / (2 * Math.PI) * orbitPeriod;
            if (timeToNode < TIME_TOLERANCE_SECONDS) { timeToNode = orbitPeriod; }
            System.out.println("Selected descending node for first burn, in " + timeToNode + " s");
        }
        AbsoluteDate firstBurnDate = initialState.getDate().shiftedBy(timeToNode);
        double burnTimeShift = 0.001;
        AbsoluteDate burnTriggerDate = firstBurnDate.shiftedBy(burnTimeShift);

        // Propagate to node to get state at burn.
        propagator = new KeplerianPropagator(initialOrbit);
        SpacecraftState stateAtNode = propagator.propagate(firstBurnDate);
        stateAtNode = new SpacecraftState(stateAtNode.getOrbit(), initialMassGlobal);
        Vector3D positionNode = stateAtNode.getPVCoordinates().getPosition();
        Vector3D velocityNode = stateAtNode.getPVCoordinates().getVelocity();
        double speedNode = velocityNode.getNorm();

        // Quick estimates for single vs. two-burn strategies.
        double deltaIncAngle = FastMath.abs(deltaInc);
        double dv_single = 2.0 * speedNode * FastMath.sin(deltaIncAngle / 2.0);
        double dv_half = 2.0 * speedNode * FastMath.sin(deltaIncAngle / 4.0);
        double dv_twoTotal = 2.0 * dv_half;
        System.out.println("One-burn ΔV:      " + dv_single + " m/s");
        System.out.println("Two-burn total ΔV:" + dv_twoTotal + " m/s");

        boolean useTwoBurns = (dv_twoTotal < dv_single) && (deltaIncAngle > FastMath.toRadians(10.0));
        System.out.println("Selected strategy: " + (useTwoBurns ? "Two burns" : "Single burn"));

        if (!useTwoBurns) {
            // SINGLE BURN MANEUVER
            Vector3D radialUnit = positionNode.normalize();
            Vector3D northUnit = new Vector3D(0, 0, 1);
            Vector3D eastUnit = northUnit.crossProduct(radialUnit).normalize();
            double signNorth = FastMath.signum(velocityNode.getZ());
            double vHorizFinal = FastMath.cos(INC_2) * speedNode;
            double vVertFinal = FastMath.sin(INC_2) * speedNode * signNorth;
            Vector3D vDesired = eastUnit.scalarMultiply(vHorizFinal).add(northUnit.scalarMultiply(vVertFinal));
            Vector3D deltaVVector = vDesired.subtract(velocityNode);
            DateDetector burnTrigger = new DateDetector(burnTriggerDate);
            ImpulseManeuver singleBurn = new ImpulseManeuver(burnTrigger, deltaVVector, ISP);
            KeplerianPropagator singleProp = new KeplerianPropagator(initialOrbit);
            singleProp.addEventDetector(singleBurn);
            SpacecraftState stateAfterBurn = singleProp.propagate(burnTriggerDate);
            double deltaV_actual = deltaVVector.getNorm();
            double massAfterBurn = calculateFinalMass(stateAtNode.getMass(), deltaV_actual, ISP, g0);
            SpacecraftState finalStateLocal = new SpacecraftState(stateAfterBurn.getOrbit(), massAfterBurn);
            double fuelConsumed = stateAtNode.getMass() - massAfterBurn;
            System.out.println("Fuel consumed (single-burn): " + fuelConsumed + " kg");
            Utils.appendLineToFile(consommationErgolsFile, String.valueOf(fuelConsumed));
            Utils.logResults(resultFileName, finalStateLocal, finalStateLocal.getOrbit(), stateAtNode.getMass(), DRYMASS);
            writeManeuverTimestamp(lastManeuverDateFile, stateAfterBurn.getDate());
            mqttService.sendFileViaMQTT(resultFileName, publishTopic);
            double newErgol = finalStateLocal.getMass() - DRYMASS;
            if (newErgol < 0) { newErgol = 0; }
            this.ERGOL = newErgol;
            this.finalState = finalStateLocal;
        } else {
            // TWO-BURN MANEUVER
            double incIntermediate = (currentInc + INC_2) / 2.0;
            Vector3D radialUnit = positionNode.normalize();
            Vector3D northUnit = new Vector3D(0, 0, 1);
            Vector3D eastUnit = northUnit.crossProduct(radialUnit).normalize();
            double signNorth = FastMath.signum(velocityNode.getZ());
            double vHorizInter = FastMath.cos(incIntermediate) * speedNode;
            double vVertInter = FastMath.sin(incIntermediate) * speedNode * signNorth;
            Vector3D vDesired1 = eastUnit.scalarMultiply(vHorizInter).add(northUnit.scalarMultiply(vVertInter));
            Vector3D deltaV1 = vDesired1.subtract(velocityNode);
            DateDetector burnTrigger1 = new DateDetector(burnTriggerDate);
            ImpulseManeuver firstBurn = new ImpulseManeuver(burnTrigger1, deltaV1, ISP);
            KeplerianPropagator firstProp = new KeplerianPropagator(initialOrbit);
            firstProp.addEventDetector(firstBurn);
            SpacecraftState stateAfterFirstBurn = firstProp.propagate(burnTriggerDate);
            double deltaV1_actual = deltaV1.getNorm();
            double massAfter1 = calculateFinalMass(stateAtNode.getMass(), deltaV1_actual, ISP, g0);
            SpacecraftState state1 = new SpacecraftState(stateAfterFirstBurn.getOrbit(), massAfter1);
            double fuelConsumed1 = stateAtNode.getMass() - massAfter1;
            System.out.println("Fuel consumed (1st burn): " + fuelConsumed1 + " kg");
            Utils.appendLineToFile(consommationErgolsFile, String.valueOf(fuelConsumed1));

            // SECOND BURN: Find next node from state1.
            AbsoluteDate secondNodeDate = findNextNodeDate(state1);
            KeplerianPropagator secondProp = new KeplerianPropagator(state1.getOrbit());
            SpacecraftState stateAtSecondNode = secondProp.propagate(secondNodeDate);
            stateAtSecondNode = new SpacecraftState(stateAtSecondNode.getOrbit(), state1.getMass());
            Vector3D position2 = stateAtSecondNode.getPVCoordinates().getPosition();
            Vector3D velocity2 = stateAtSecondNode.getPVCoordinates().getVelocity();
            double speed2 = velocity2.getNorm();
            Vector3D radial2 = position2.normalize();
            Vector3D east2 = new Vector3D(0, 0, 1).crossProduct(radial2).normalize();
            double signNorth2 = FastMath.signum(velocity2.getZ());
            double vHorizFinal2 = FastMath.cos(INC_2) * speed2;
            double vVertFinal2 = FastMath.sin(INC_2) * speed2 * signNorth2;
            Vector3D vDesired2 = east2.scalarMultiply(vHorizFinal2).add(new Vector3D(0, 0, 1).scalarMultiply(vVertFinal2));
            Vector3D deltaV2 = vDesired2.subtract(velocity2);
            DateDetector burnTrigger2 = new DateDetector(secondNodeDate.shiftedBy(burnTimeShift));
            ImpulseManeuver secondBurn = new ImpulseManeuver(burnTrigger2, deltaV2, ISP);
            KeplerianPropagator secondBurnProp = new KeplerianPropagator(stateAtSecondNode.getOrbit());
            secondBurnProp.addEventDetector(secondBurn);
            SpacecraftState stateAfterSecondBurn = secondBurnProp.propagate(secondNodeDate.shiftedBy(burnTimeShift));
            double deltaV2_actual = deltaV2.getNorm();
            double massAfter2 = calculateFinalMass(stateAtSecondNode.getMass(), deltaV2_actual, ISP, g0);
            SpacecraftState finalStateLocal = new SpacecraftState(stateAfterSecondBurn.getOrbit(), massAfter2);
            double fuelConsumed2 = stateAtSecondNode.getMass() - massAfter2;
            double totalFuelConsumed = fuelConsumed1 + fuelConsumed2;
            System.out.println("Fuel consumed (2nd burn): " + fuelConsumed2 + " kg");
            System.out.println("Total fuel consumed: " + totalFuelConsumed + " kg");
            Utils.appendLineToFile(consommationErgolsFile, String.valueOf(totalFuelConsumed));
            Utils.logResults(resultFileName, finalStateLocal, finalStateLocal.getOrbit(), stateAtNode.getMass(), DRYMASS);
            writeManeuverTimestamp(lastManeuverDateFile, finalStateLocal.getDate());
            mqttService.sendFileViaMQTT(resultFileName, publishTopic);
            double newErgol = finalStateLocal.getMass() - DRYMASS;
            if (newErgol < 0) { newErgol = 0; }
            this.ERGOL = newErgol;
            this.finalState = finalStateLocal;
        }

        double endTime = System.currentTimeMillis();
        System.out.println("Inclination maneuver completed in " + (endTime - startTime) + " ms");
    }

    @Override
    public void calculateErgolConsumption() throws IOException {
        // Calculate ergol (fuel) consumption for the planned inclination change maneuver.
        double initialMass = DRYMASS + ERGOL;
        double deltaInc = FastMath.abs(INC_2 - INC);
        // For a near-circular orbit, orbital speed = sqrt(MU/SMA)
        double orbitalSpeed = FastMath.sqrt(MU / SMA);
        // Single-burn deltaV required
        double dv_single = 2 * orbitalSpeed * FastMath.sin(deltaInc / 2.0);
        // Two-burn: two equal burns
        double dv_half = 2 * orbitalSpeed * FastMath.sin(deltaInc / 4.0);
        double dv_twoTotal = 2 * dv_half;

        // Compute final mass using the rocket equation for both strategies
        double finalMass_single = initialMass * FastMath.exp(-dv_single / (ISP * g0));
        double fuelConsumed_single = initialMass - finalMass_single;

        double finalMass_firstBurn = initialMass * FastMath.exp(-dv_half / (ISP * g0));
        double finalMass_two = finalMass_firstBurn * FastMath.exp(-dv_half / (ISP * g0));
        double fuelConsumed_two = initialMass - finalMass_two;

        boolean useTwoBurns = (dv_twoTotal < dv_single) && (deltaInc > FastMath.toRadians(10.0));

        System.out.println("Calculated ergol consumption estimates:");
        System.out.println("Single-burn ΔV: " + dv_single + " m/s, fuel consumed: " + fuelConsumed_single + " kg");
        System.out.println("Two-burn total ΔV: " + dv_twoTotal + " m/s, fuel consumed: " + fuelConsumed_two + " kg");

        if (useTwoBurns) {
            System.out.println("Optimal strategy: Two burns");
            System.out.println("Ergol consumption: " + fuelConsumed_two + " kg");
            Utils.appendLineToFile(consommationErgolsFile, String.valueOf(fuelConsumed_two));
        } else {
            System.out.println("Optimal strategy: Single burn");
            System.out.println("Ergol consumption: " + fuelConsumed_single + " kg");
            Utils.appendLineToFile(consommationErgolsFile, String.valueOf(fuelConsumed_single));
        }
    }

    @Override
    public void processReachOrbitTime() throws IOException {
        // This method calculates the "reach orbit time" – the time when the new orbit is fully established.
        // For an impulsive maneuver on a near-circular orbit, a common heuristic is to propagate from the initial date
        // (parsed from the command data) for half an orbital period.
        try {
            printAttributes();
            double start = System.currentTimeMillis();
            // 1) Parse the initial date
            AbsoluteDate initialDate = Utils.parseDateFromTimestamp(DATE);
            System.out.println("Initial date: " + initialDate);
            // 2) Validate orbit parameters
            validateOrbitParameters();
            // 3) Reconstruct the final orbit using the target inclination (assuming SMA remains unchanged).
            Frame eme2000 = FramesFactory.getEME2000();
            KeplerianOrbit finalOrbit = new KeplerianOrbit(SMA, ECC, INC_2, PA, RAAN, ANO, PositionAngle.MEAN, eme2000, initialDate, MU);
            double period = finalOrbit.getKeplerianPeriod();
            // 4) Define reach orbit time as initialDate + half an orbital period.
            AbsoluteDate reachOrbitTime = initialDate.shiftedBy(period / 2.0);
            // 5) Log it.
            Utils.logApsideDate(apsideFile, reachOrbitTime);
            System.out.println("Reach Orbit Time: " + reachOrbitTime);
            double duration = (System.currentTimeMillis() - start) / 1000.0;
            System.out.println("Execution time: " + duration + " s");
        } catch (Exception e) {
            System.err.println("Error in reach orbit time calculation: " + e.getMessage());
            throw e;
        }
    }
}
*/

//
//
//@Override
//public void calculateErgolConsumption() throws IOException {
//    // Calculate ergol (fuel) consumption for the planned inclination change maneuver.
//    double initialMass = DRYMASS + ERGOL;
//    double deltaInc = FastMath.abs(INC_2 - INC);
//    // For a near-circular orbit, orbital speed = sqrt(MU/SMA)
//    double orbitalSpeed = FastMath.sqrt(MU / SMA);
//    // Single-burn deltaV required
//    double dv_single = 2 * orbitalSpeed * FastMath.sin(deltaInc / 2.0);
//    // Two-burn: two equal burns
//    double dv_half = 2 * orbitalSpeed * FastMath.sin(deltaInc / 4.0);
//    double dv_twoTotal = 2 * dv_half;
//
//    // Compute final mass using the rocket equation for both strategies
//    double finalMass_single = initialMass * FastMath.exp(-dv_single / (ISP * g0));
//    double fuelConsumed_single = initialMass - finalMass_single;
//
//    double finalMass_firstBurn = initialMass * FastMath.exp(-dv_half / (ISP * g0));
//    double finalMass_two = finalMass_firstBurn * FastMath.exp(-dv_half / (ISP * g0));
//    double fuelConsumed_two = initialMass - finalMass_two;
//
//    boolean useTwoBurns = (dv_twoTotal < dv_single) && (deltaInc > FastMath.toRadians(10.0));
//
//    System.out.println("Calculated ergol consumption estimates:");
//    System.out.println("Single-burn ΔV: " + dv_single + " m/s, fuel consumed: " + fuelConsumed_single + " kg");
//    System.out.println("Two-burn total ΔV: " + dv_twoTotal + " m/s, fuel consumed: " + fuelConsumed_two + " kg");
//
//    if (useTwoBurns) {
//        System.out.println("Optimal strategy: Two burns");
//        System.out.println("Ergol consumption: " + fuelConsumed_two + " kg");
//        Utils.appendLineToFile(consommationErgolsFile, String.valueOf(fuelConsumed_two));
//    } else {
//        System.out.println("Optimal strategy: Single burn");
//        System.out.println("Ergol consumption: " + fuelConsumed_single + " kg");
//        Utils.appendLineToFile(consommationErgolsFile, String.valueOf(fuelConsumed_single));
//    }
//}
//
//@Override
//public void processReachOrbitTime() throws IOException {
//    // This method now calculates the exact moment when the burn will be executed
//    // instead of a theoretical time based on half orbital period
//    try {
//        printAttributes();
//        double start = System.currentTimeMillis();
//
//        // 1) Parse the initial date from command
//        AbsoluteDate initialDate = Utils.parseDateFromTimestamp(DATE);
//        System.out.println("Initial date: " + initialDate);
//
//        // 2) Validate orbit parameters
//        validateOrbitParameters();
//
//        // 3) Setup environment and initial orbit
//        Frame eme2000 = FramesFactory.getEME2000();
//        KeplerianOrbit initialOrbit = new KeplerianOrbit(
//                SMA, ECC, INC, PA, RAAN, ANO, PositionAngle.MEAN, eme2000, initialDate, MU
//        );
//
//        // 4) Calculate time to node crossing (exactly like in computeAndExecute)
//        double currentInc = initialOrbit.getI();
//        double deltaInc = INC_2 - currentInc;
//        boolean incIncreasing = deltaInc > 0;
//
//        double trueAnomaly = initialOrbit.getTrueAnomaly();
//        double raan = initialOrbit.getRightAscensionOfAscendingNode();
//        double aop = initialOrbit.getPerigeeArgument();
//        double taAscending = -(raan + aop) % (2 * Math.PI);
//        if (taAscending < 0) { taAscending += 2 * Math.PI; }
//        double taDescending = (taAscending + Math.PI) % (2 * Math.PI);
//        double orbitPeriod = initialOrbit.getKeplerianPeriod();
//        double diffAsc = (taAscending - trueAnomaly + 2 * Math.PI) % (2 * Math.PI);
//        double diffDesc = (taDescending - trueAnomaly + 2 * Math.PI) % (2 * Math.PI);
//
//        // 5) Determine time to node based on whether inclination is increasing/decreasing
//        double timeToNode;
//        if (incIncreasing) {
//            timeToNode = diffAsc / (2 * Math.PI) * orbitPeriod;
//            if (timeToNode < TIME_TOLERANCE_SECONDS) { timeToNode = orbitPeriod; }
//            System.out.println("Selected ascending node for burn, in " + timeToNode + " s");
//        } else {
//            timeToNode = diffDesc / (2 * Math.PI) * orbitPeriod;
//            if (timeToNode < TIME_TOLERANCE_SECONDS) { timeToNode = orbitPeriod; }
//            System.out.println("Selected descending node for burn, in " + timeToNode + " s");
//        }
//
//        // 6) Calculate burn execution date
//        AbsoluteDate burnExecutionDate = initialDate.shiftedBy(timeToNode);
//        double burnTimeShift = 0.001; // Small shift for burn trigger
//        AbsoluteDate burnTriggerDate = burnExecutionDate.shiftedBy(burnTimeShift);
//
//        // 7) Check for two-burn strategy
//        double orbitalSpeed = FastMath.sqrt(MU / SMA);
//        double deltaIncAngle = FastMath.abs(deltaInc);
//        double dv_single = 2.0 * orbitalSpeed * FastMath.sin(deltaIncAngle / 2.0);
//        double dv_half = 2.0 * orbitalSpeed * FastMath.sin(deltaIncAngle / 4.0);
//        double dv_twoTotal = 2.0 * dv_half;
//        boolean useTwoBurns = (dv_twoTotal < dv_single) && (deltaIncAngle > FastMath.toRadians(10.0));
//
//        // 8) For two-burn strategy, calculate the second burn time
//        AbsoluteDate finalBurnDate;
//        if (useTwoBurns) {
//            // For two-burn scenario, we need to estimate when the second burn would happen
//            // This is more complex but we'll use a simplified approach based on orbital period
//            System.out.println("Two-burn strategy selected");
//            // Simple estimate: It takes about half an orbit to go from one node to the next
//            finalBurnDate = burnExecutionDate.shiftedBy(orbitPeriod / 2.0);
//        } else {
//            // For single-burn, the final burn date is the same as the first burn
//            System.out.println("Single-burn strategy selected");
//            finalBurnDate = burnExecutionDate;
//        }
//
//        // 9) Log the actual burn execution time
//        Utils.logApsideDate(apsideFile, finalBurnDate);
//        System.out.println("Actual Burn Execution Time: " + finalBurnDate);
//
//        double duration = (System.currentTimeMillis() - start) / 1000.0;
//        System.out.println("Execution time: " + duration + " s");
//    } catch (Exception e) {
//        System.err.println("Error in reach orbit time calculation: " + e.getMessage());
//        throw e;
//    }
//}
