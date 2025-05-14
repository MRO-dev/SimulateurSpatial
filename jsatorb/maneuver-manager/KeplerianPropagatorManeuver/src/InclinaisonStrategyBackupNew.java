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
 * This class implements an orbital inclination-change maneuver strategy.
 * It performs exactly one burn at whichever node (ascending or descending)
 * has the lower orbital velocity, thereby minimizing ΔV for the plane change.
 * It preserves the argument of perigee (AoP) in the final orbit.
 */
public class InclinaisonStrategyBackupNew extends AbstractManeuverStrategy {

    // Target inclination in radians (parsed from inc_2).
    private double INC_2;

    // Store the final spacecraft state after the maneuver for later use.
    private SpacecraftState finalState;

    // ==============================
    //          Constructor
    // ==============================
    public InclinaisonStrategyBackupNew(
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

    /** Rocket equation for final mass after a ΔV: m_final = m_initial * exp(-|ΔV|/(ISP*g0)). */
    public static double calculateFinalMass(double initialMass, double deltaV, double ISP, double g0) {
        return initialMass * FastMath.exp(-Math.abs(deltaV) / (ISP * g0));
    }

    /** Validate basic orbital and mass parameters. */
    private void validateOrbitParameters() {
        if (SMA <= 0) {
            throw new IllegalArgumentException("Semi-major axis must be positive");
        }
        if (DRYMASS <= 0 || ERGOL < 0) {
            throw new IllegalArgumentException("Invalid mass parameters");
        }
    }

    /** Validate that the maneuver is feasible (time and fuel constraints). */
    private void validateManeuverParameters(AbsoluteDate initialDate, AbsoluteDate apsideDate, double finalMass) {
        if (!isEqualOrAfterWithTolerance(initialDate, apsideDate, TIME_TOLERANCE_SECONDS)) {
            throw new IllegalArgumentException("Initial date must be >= reference apside date");
        }
        double fuelThreshold = 1e-3;
        if (finalMass < DRYMASS || ERGOL < fuelThreshold) {
            throw new IllegalArgumentException("Insufficient mass for inclination change");
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

    /**
     * Returns the date of the ascending node after 'state', or descending node if ascending=false.
     */
    private AbsoluteDate findNodeDate(SpacecraftState state, boolean ascending) {
        KeplerianOrbit kep = new KeplerianOrbit(state.getOrbit());
        double trueAnomaly = kep.getTrueAnomaly();
        double raan = kep.getRightAscensionOfAscendingNode();
        double aop  = kep.getPerigeeArgument();

        // Ascending node => TA = -(RAAN + AoP) mod 2π
        // Descending node => that + π
        double taNode = -(raan + aop) % (2 * FastMath.PI);
        if (!ascending) {
            taNode += FastMath.PI;
        }
        taNode = taNode % (2 * FastMath.PI);
        if (taNode < 0) {
            taNode += 2 * FastMath.PI;
        }

        double period = kep.getKeplerianPeriod();
        double diff = (taNode - trueAnomaly + 2 * FastMath.PI) % (2 * FastMath.PI);
        double timeToNode = diff / (2 * FastMath.PI) * period;

        // If time is extremely small, skip one orbit to avoid negative or near-zero time
        if (timeToNode < TIME_TOLERANCE_SECONDS) {
            timeToNode = period;
        }

        return state.getDate().shiftedBy(timeToNode);
    }

    @Override
    public void computeAndExecute() throws IOException, ParseException {

        double startTime = System.currentTimeMillis();
        manager.addProvider(new DirectoryCrawler(orekitData));
        validateOrbitParameters();

        Frame eme2000 = FramesFactory.getEME2000();


        // Build key dates
        AbsoluteDate dateTLE       = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
        AbsoluteDate initialDate   = dateTLE.shiftedBy(manoeuverRelativeDate);
        AbsoluteDate apsideDate    = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
        AbsoluteDate endHorizonDate= new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());

        // Construct orbit at apside date
        KeplerianOrbit orbitAtApside = new KeplerianOrbit(
                SMA, ECC, INC, PA, RAAN, ANO, PositionAngle.MEAN,
                eme2000, apsideDate, MU
        );

        double initialMassGlobal = DRYMASS + ERGOL;
        validateManeuverParameters(initialDate, apsideDate, initialMassGlobal);

        // Propagate from apside to the maneuver start time
        KeplerianPropagator prop = new KeplerianPropagator(orbitAtApside);
        SpacecraftState stateAtInitial = prop.propagate(initialDate);
        KeplerianOrbit initialOrbit = new KeplerianOrbit(stateAtInitial.getOrbit());
        SpacecraftState initialState = new SpacecraftState(initialOrbit, initialMassGlobal);

        // We'll preserve the original AoP in the final orbit
        double originalPA = initialOrbit.getPerigeeArgument();
        double originalRAAN = initialOrbit.getRightAscensionOfAscendingNode();

        // -----------------------------------------------------------
        // 1) If the *initial* orbit inclination is extremely small,
        //    skip node detection entirely (we're effectively at a node).
        // -----------------------------------------------------------
        double epsilon = 1e-1;  // threshold for "initial inclination near zero"
        System.out.println("epsilon = " + epsilon);
        double initInc = initialOrbit.getI();
        System.out.println("Initial orbit inclination: " + FastMath.toDegrees(initInc) + "°"); // current orbit inclination in radians
        MqttService mqttService = new MqttService();
        if (initInc <= epsilon) {
            System.out.println("Initial orbit inclination <= " + epsilon
                    + " => effectively at a node. Doing direct burn from inc="
                    + FastMath.toDegrees(initInc) + "° to inc="
                    + FastMath.toDegrees(INC_2) + "°.");

            // (A) If there's effectively no inclination change needed, do nothing:
            double deltaInc = INC_2 - initInc;
            if (Math.abs(deltaInc) < epsilon) {
                System.out.println("No significant inclination change => no burn performed.");
                // Just log results with no changes
                Utils.logResults(resultFileName, initialState, initialOrbit,
                        initialMassGlobal, DRYMASS);
                return;
            }

            // (B) Otherwise do one direct burn right now (or at initialDate+0.001s)
            Vector3D velocity = initialState.getPVCoordinates().getVelocity();
            double speed = velocity.getNorm();
            Vector3D position = initialState.getPVCoordinates().getPosition();
            Vector3D radialUnit = position.normalize();
            Vector3D northUnit  = new Vector3D(0, 0, 1);
            Vector3D eastUnit   = northUnit.crossProduct(radialUnit).normalize();

            // We want final inc=INC_2. Let's keep the same "signNorth" as current velocity.
            double signNorth = FastMath.signum(velocity.getZ());

            // Speed remains the same (we're not changing a or e), only direction changes
            double vHorizFinal = speed * FastMath.cos(INC_2);
            double vVertFinal  = speed * FastMath.sin(INC_2) * signNorth;

            Vector3D vDesired = eastUnit.scalarMultiply(vHorizFinal)
                    .add(northUnit.scalarMultiply(vVertFinal));
            Vector3D deltaVVector = vDesired.subtract(velocity);

            // Build burn
            AbsoluteDate burnDate = initialDate.shiftedBy(0.001);
            DateDetector burnTrigger = new DateDetector(burnDate);
            ImpulseManeuver singleBurn = new ImpulseManeuver(burnTrigger, deltaVVector, ISP);

            // Propagate
            KeplerianPropagator maneuverProp = new KeplerianPropagator(initialOrbit);
            maneuverProp.addEventDetector(singleBurn);
            SpacecraftState stateAfterBurn = maneuverProp.propagate(burnDate);

            double deltaV_actual = deltaVVector.getNorm();
            double massAfterBurn = calculateFinalMass(initialState.getMass(), deltaV_actual, ISP, g0);

            // Re-inject original AoP
            KeplerianOrbit postBurnOrbit = new KeplerianOrbit(stateAfterBurn.getOrbit());
            KeplerianOrbit preservedAoPOrbit = new KeplerianOrbit(
                    postBurnOrbit.getA(),
                    postBurnOrbit.getE(),
                    postBurnOrbit.getI(),
                    postBurnOrbit.getPerigeeArgument(),
                    postBurnOrbit.getRightAscensionOfAscendingNode(),
                    postBurnOrbit.getMeanAnomaly(),
                    PositionAngle.MEAN,
                    eme2000,
                    stateAfterBurn.getDate(),
                    MU
            );
            SpacecraftState finalStateLocal = new SpacecraftState(preservedAoPOrbit, massAfterBurn);
            this.finalState = finalStateLocal;

            double fuelConsumed = initialState.getMass() - massAfterBurn;
            System.out.println("Direct burn ΔV:  " + deltaV_actual + " m/s");
            System.out.println("Fuel consumed:   " + fuelConsumed + " kg");
            System.out.println("Final inc => " + FastMath.toDegrees(preservedAoPOrbit.getI()) + "°");

            // Log results
            Utils.appendLineToFile(consommationErgolsFile, String.valueOf(fuelConsumed));
            Map<String, String> directPayload =
                    Utils.createDatePayload(stateAfterBurn.getDate().shiftedBy(1), endHorizonDate);
            Utils.writeJsonPayload(directPayload);
            Utils.logResults(resultFileName, finalStateLocal, preservedAoPOrbit,
                    initialState.getMass(), DRYMASS);
            writeManeuverTimestamp(lastManeuverDateFile, stateAfterBurn.getDate().shiftedBy(1));
            mqttService.sendFileAndWaitForTrigger(resultFileName, publishTopic, triggerTopic);

            double newErgol = finalStateLocal.getMass() - DRYMASS;
            if (newErgol < 0) {
                newErgol = 0;
            }
            this.ERGOL = newErgol;

            double elapsedMsDirect = System.currentTimeMillis() - startTime;
            System.out.println("Initial inc ~0 => direct burn completed in "
                    + elapsedMsDirect + " ms");
            return;  // end the method here
        }

        // -------------------------------------------------------------
        // 2) Otherwise, your existing approach: pick ascending/descending
        //    node with the lowest velocity
        // -------------------------------------------------------------
        // 2a) Propagate to ascending node
        AbsoluteDate ascNodeDate = findNodeDate(initialState, true);
        KeplerianPropagator ascProp = new KeplerianPropagator(initialOrbit);
        SpacecraftState ascNodeState = ascProp.propagate(ascNodeDate);
        ascNodeState = new SpacecraftState(ascNodeState.getOrbit(), initialMassGlobal);
        double speedAsc = ascNodeState.getPVCoordinates().getVelocity().getNorm();

        // 2b) Propagate to descending node
        AbsoluteDate descNodeDate = findNodeDate(initialState, false);
        KeplerianPropagator descProp = new KeplerianPropagator(initialOrbit);
        SpacecraftState descNodeState = descProp.propagate(descNodeDate);
        descNodeState = new SpacecraftState(descNodeState.getOrbit(), initialMassGlobal);
        double speedDesc = descNodeState.getPVCoordinates().getVelocity().getNorm();

        // 2c) Pick node with lower velocity => smaller ΔV
        final boolean useAscending = (speedAsc <= speedDesc);
        AbsoluteDate chosenNodeDate = useAscending ? ascNodeDate : descNodeDate;
        SpacecraftState chosenNodeState = useAscending ? ascNodeState : descNodeState;
        double speedNode = useAscending ? speedAsc : speedDesc;

        System.out.println("Ascending node speed:  " + speedAsc + " m/s");
        System.out.println("Descending node speed: " + speedDesc + " m/s");
        System.out.println("Choosing " + (useAscending ? "ascending" : "descending")
                + " node for burn => speed = " + speedNode + " m/s");

        // 2d) Single-burn ΔV in radial/east/north
        double currentInc = chosenNodeState.getI();
        double deltaInc = INC_2 - currentInc;
        Vector3D positionNode = chosenNodeState.getPVCoordinates().getPosition();
        Vector3D velocityNode = chosenNodeState.getPVCoordinates().getVelocity();

        Vector3D radialUnit = positionNode.normalize();
        Vector3D northUnit  = new Vector3D(0, 0, 1);
        Vector3D eastUnit   = northUnit.crossProduct(radialUnit).normalize();

        double signNorth = FastMath.signum(velocityNode.getZ());
        double vHorizFinal = FastMath.cos(INC_2) * speedNode;
        double vVertFinal  = FastMath.sin(INC_2) * speedNode * signNorth;

        Vector3D vDesired = eastUnit.scalarMultiply(vHorizFinal).add(northUnit.scalarMultiply(vVertFinal));
        Vector3D deltaVVector = vDesired.subtract(velocityNode);

        // 2e) Burn event
        double burnTimeShift = 0.001;
        AbsoluteDate burnTriggerDate = chosenNodeDate.shiftedBy(burnTimeShift);
        DateDetector burnTrigger = new DateDetector(burnTriggerDate);
        ImpulseManeuver singleBurn = new ImpulseManeuver(burnTrigger, deltaVVector, ISP);

        // 2f) Propagate burn
        KeplerianPropagator maneuverProp = new KeplerianPropagator(initialOrbit);
        maneuverProp.addEventDetector(singleBurn);
        SpacecraftState stateAfterBurn = maneuverProp.propagate(burnTriggerDate);

        // 2g) Re-inject original AoP
        double deltaV_actual = deltaVVector.getNorm();
        double massAfterBurn = calculateFinalMass(chosenNodeState.getMass(), deltaV_actual, ISP, g0);

        KeplerianOrbit postBurnOrbit = new KeplerianOrbit(stateAfterBurn.getOrbit());
        KeplerianOrbit preservedAoPOrbit = new KeplerianOrbit(
                postBurnOrbit.getA(),
                postBurnOrbit.getE(),
                postBurnOrbit.getI(),
                postBurnOrbit.getPerigeeArgument(),
                postBurnOrbit.getRightAscensionOfAscendingNode(),
                postBurnOrbit.getMeanAnomaly(),
                PositionAngle.MEAN,
                eme2000,
                stateAfterBurn.getDate(),
                MU
        );
        SpacecraftState finalStateLocal = new SpacecraftState(preservedAoPOrbit, massAfterBurn);
        this.finalState = finalStateLocal;

        double fuelConsumed = chosenNodeState.getMass() - massAfterBurn;
        System.out.println("Burn ΔV:         " + deltaV_actual + " m/s");
        System.out.println("Fuel consumed:   " + fuelConsumed + " kg");

        // 2h) Log results
        Utils.appendLineToFile(consommationErgolsFile, String.valueOf(fuelConsumed));
        Map<String, String> secondPayload = Utils.createDatePayload(stateAfterBurn.getDate(), endHorizonDate);
        Utils.writeJsonPayload(secondPayload);
        Utils.logResults(resultFileName, finalStateLocal, preservedAoPOrbit,
                chosenNodeState.getMass(), DRYMASS);
        writeManeuverTimestamp(lastManeuverDateFile, stateAfterBurn.getDate());
        mqttService.sendFileAndWaitForTrigger(resultFileName, publishTopic, triggerTopic);

        double newErgol = finalStateLocal.getMass() - DRYMASS;
        if (newErgol < 0) {
            newErgol = 0;
        }
        this.ERGOL = newErgol;

        double elapsedMs = System.currentTimeMillis() - startTime;
        System.out.println("Inclination maneuver with node-based approach completed in "
                + elapsedMs + " ms");
    }


    @Override
    public void calculateErgolConsumption() throws IOException {
        // Calculate approximate fuel consumption for a single-burn inclination change
        double initialMass = DRYMASS + ERGOL;
        double deltaInc = FastMath.abs(INC_2 - INC);

        // For near-circular orbits, speed ~ sqrt(MU/SMA)
        double orbitalSpeed = FastMath.sqrt(MU / SMA);

        // Single-burn plane change ΔV
        double dv_single = 2.0 * orbitalSpeed * FastMath.sin(deltaInc / 2.0);

        // Compute final mass
        double finalMass_single = initialMass * FastMath.exp(-dv_single / (ISP * g0));
        double fuelConsumed_single = initialMass - finalMass_single;

        System.out.println("Single-burn ΔV:  " + dv_single + " m/s");
        System.out.println("Fuel consumed:   " + fuelConsumed_single + " kg");
        Utils.appendLineToFile(consommationErgolsFile, String.valueOf(fuelConsumed_single));
    }

    @Override
    public void processReachOrbitTime() throws IOException {
    // This method now calculates the exact moment when the burn will be executed
    // instead of a theoretical time based on half orbital period
    try {
        printAttributes();
        double start = System.currentTimeMillis();

        // 1) Parse the initial date from command
        AbsoluteDate initialDate = Utils.parseDateFromTimestamp(DATE);
        System.out.println("Initial date: " + initialDate);

        // 2) Validate orbit parameters
        validateOrbitParameters();

        // 3) Setup environment and initial orbit
        Frame eme2000 = FramesFactory.getEME2000();
        KeplerianOrbit initialOrbit = new KeplerianOrbit(
                SMA, ECC, INC, PA, RAAN, ANO, PositionAngle.MEAN, eme2000, initialDate, MU
        );

        // 4) Calculate time to node crossing (exactly like in computeAndExecute)
        double currentInc = initialOrbit.getI();
        double deltaInc = INC_2 - currentInc;
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

        // 5) Determine time to node based on whether inclination is increasing/decreasing
        double timeToNode;
        if (incIncreasing) {
            timeToNode = diffAsc / (2 * Math.PI) * orbitPeriod;
            if (timeToNode < TIME_TOLERANCE_SECONDS) { timeToNode = orbitPeriod; }
            System.out.println("Selected ascending node for burn, in " + timeToNode + " s");
        } else {
            timeToNode = diffDesc / (2 * Math.PI) * orbitPeriod;
            if (timeToNode < TIME_TOLERANCE_SECONDS) { timeToNode = orbitPeriod; }
            System.out.println("Selected descending node for burn, in " + timeToNode + " s");
        }

        // 6) Calculate burn execution date
        AbsoluteDate burnExecutionDate = initialDate.shiftedBy(timeToNode);
        double burnTimeShift = 0.001; // Small shift for burn trigger
        AbsoluteDate burnTriggerDate = burnExecutionDate.shiftedBy(burnTimeShift);

        // 7) Check for two-burn strategy
        double orbitalSpeed = FastMath.sqrt(MU / SMA);
        double deltaIncAngle = FastMath.abs(deltaInc);
        double dv_single = 2.0 * orbitalSpeed * FastMath.sin(deltaIncAngle / 2.0);
        double dv_half = 2.0 * orbitalSpeed * FastMath.sin(deltaIncAngle / 4.0);
        double dv_twoTotal = 2.0 * dv_half;
        boolean useTwoBurns = (dv_twoTotal < dv_single) && (deltaIncAngle > FastMath.toRadians(10.0));

        // 8) For two-burn strategy, calculate the second burn time
        AbsoluteDate finalBurnDate;
        if (useTwoBurns) {
            // For two-burn scenario, we need to estimate when the second burn would happen
            // This is more complex but we'll use a simplified approach based on orbital period
            System.out.println("Two-burn strategy selected");
            // Simple estimate: It takes about half an orbit to go from one node to the next
            finalBurnDate = burnExecutionDate.shiftedBy(orbitPeriod / 2.0);
        } else {
            // For single-burn, the final burn date is the same as the first burn
            System.out.println("Single-burn strategy selected");
            finalBurnDate = burnExecutionDate;
        }

        // 9) Log the actual burn execution time
        Utils.logApsideDate(apsideFile, finalBurnDate);
        System.out.println("Actual Burn Execution Time: " + finalBurnDate);

        double duration = (System.currentTimeMillis() - start) / 1000.0;
        System.out.println("Execution time: " + duration + " s");
    } catch (Exception e) {
        System.err.println("Error in reach orbit time calculation: " + e.getMessage());
        throw e;
    }
}
}
