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
import java.util.Map;

*/
/**
 * This class implements an orbital inclination-change maneuver strategy.
 * It performs exactly one burn at whichever node (ascending or descending)
 * has the lower orbital velocity, thereby minimizing ΔV for the plane change.
 * It preserves the argument of perigee (AoP) in the final orbit.
 *//*

public class InclinaisonStrategy extends AbstractManeuverStrategy {

    // Target inclination in radians (parsed from inc_2).
    private double INC_2;

    // Store the final spacecraft state after the maneuver for later use.
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
            double inc_2,//parametre final inclinaison 2 (degrés)
            String dateStr,
            String modeParameter,
            String commandParameter
    ) throws IOException {
        super(modeParameter, commandParameter);
        this.SMA = sma;
        this.ECC = ecc;
        this.INC = FastMath.toRadians(incDeg);
        this.RAAN = FastMath.toRadians(raanDeg);
        this.PA = FastMath.toRadians(aopDeg);
        this.ANO = FastMath.toRadians(meanAnomDeg);
        this.DRYMASS = dryMass;
        this.ERGOL = ergol;
        this.ISP = isp;
        this.INC_2 = FastMath.toRadians(inc_2);
        this.DATE = dateStr;
    }

    */
/**
     * Rocket equation for final mass after a ΔV: m_final = m_initial * exp(-|ΔV|/(ISP*g0)).
     *//*

    public static double calculateFinalMass(double initialMass, double deltaV, double ISP, double g0) {
        return initialMass * FastMath.exp(-Math.abs(deltaV) / (ISP * g0));
    }

    */
/**
     * Validate basic orbital and mass parameters.
     *//*

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

    */
/**
     * Returns the date of the ascending node after 'state', or descending node if ascending=false.
     *//*

    private AbsoluteDate findNodeDate(SpacecraftState state, boolean ascending) {
        KeplerianOrbit kep = new KeplerianOrbit(state.getOrbit());
        double trueAnomaly = kep.getTrueAnomaly();
        double raan = kep.getRightAscensionOfAscendingNode();
        double aop = kep.getPerigeeArgument();

        // Ascending node => TA = -(RAAN + AoP) mod 2π
        // Descending node => that + π
        double taNode = -(raan + aop) % (2 * FastMath.PI);
        if (!ascending) {
            taNode += FastMath.PI;
        }
        taNode = taNode % (2 * FastMath.PI);
        // Si taNode < 0 : on ajoute 2π
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

        // System.out.println("timeToNode = " + timeToNode);
        return state.getDate().shiftedBy(timeToNode);
    }

    */
/**
     * Propagate from a known origin orbit to a target date => get state.
     *//*

    private SpacecraftState propagateToDate(KeplerianOrbit originOrbit, AbsoluteDate date) {
        KeplerianPropagator kp = new KeplerianPropagator(originOrbit);
        return kp.propagate(date);
    }

    // Execute la maneuvre de changement d'inclinaison
    @Override
    public void computeAndExecute() throws IOException, ParseException {

        double startTime = System.currentTimeMillis();
        manager.addProvider(new DirectoryCrawler(orekitData));
        validateOrbitParameters();

        Frame eme2000 = FramesFactory.getEME2000();
        MqttService mqttService = new MqttService();

        // Build times
        AbsoluteDate dateTLE = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
        AbsoluteDate initialDate = dateTLE.shiftedBy(manoeuverRelativeDate);
        AbsoluteDate apsideDate = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
        AbsoluteDate endHorizonDate = new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());

        // Construct orbit at apside date
        KeplerianOrbit orbitAtApside = new KeplerianOrbit(
                SMA, ECC, INC, PA, RAAN, ANO, PositionAngle.MEAN,
                eme2000, apsideDate, MU
        );

        double initialMass = DRYMASS + ERGOL;
        validateManeuverParameters(initialDate, apsideDate, initialMass);

        // Propagate to initial date
        KeplerianPropagator prop = new KeplerianPropagator(orbitAtApside);
        SpacecraftState stateAtInitial = prop.propagate(initialDate);
        KeplerianOrbit initialOrbit = new KeplerianOrbit(stateAtInitial.getOrbit());
        SpacecraftState initialState = new SpacecraftState(initialOrbit, initialMass);

        // -------------- Main Logic --------------
        double incEps = 1e-3;   // threshold for near-0 inc (radians)
        double eccThreshold = 0.05;   // "high e" threshold
        double currentInc = initialOrbit.getI();

        System.out.println("Current inc= " + FastMath.toDegrees(currentInc)
                + "°, ECC= " + ECC);

        // Case A: near eq + high e => wait for apogee
        if (currentInc <= incEps && ECC >= eccThreshold) {
            System.out.println("Case A: near-eq + high e => wait for apogee.");

            // 1) find next apogee
            AbsoluteDate apogeeDate = findNextApogee(initialState);
            double smallShift = 1.0;
            AbsoluteDate burnDate = apogeeDate.shiftedBy(smallShift);

            // 2) get velocity at apogee
            SpacecraftState stateApogee = propagateToDate(initialOrbit, apogeeDate);
            Vector3D velocityApg = stateApogee.getPVCoordinates().getVelocity();
            double speedApg = velocityApg.getNorm();

            // Build local radial/east/north
            Vector3D positionApg = stateApogee.getPVCoordinates().getPosition();
            Vector3D radialUnit = positionApg.normalize();
            Vector3D northUnit = new Vector3D(0, 0, 1);
            Vector3D eastUnit = northUnit.crossProduct(radialUnit).normalize();
            double signNorth = FastMath.signum(velocityApg.getZ());

            double vHorizFinal = speedApg * FastMath.cos(INC_2);
            double vVertFinal = speedApg * FastMath.sin(INC_2) * signNorth;

            Vector3D vDesired = eastUnit.scalarMultiply(vHorizFinal)
                    .add(northUnit.scalarMultiply(vVertFinal));
            Vector3D deltaVVector = vDesired.subtract(velocityApg);

            // build impulse
            DateDetector burnTrigger = new DateDetector(burnDate);
            ImpulseManeuver singleBurn = new ImpulseManeuver(burnTrigger, deltaVVector, ISP);

            // propagate
            KeplerianPropagator mp = new KeplerianPropagator(initialOrbit);
            mp.addEventDetector(singleBurn);
            SpacecraftState stateAfterBurn = mp.propagate(burnDate);

            double dvApogee = deltaVVector.getNorm();
            double massAfterBurn = calculateFinalMass(initialMass, dvApogee, ISP, g0);

            KeplerianOrbit orbitAfter = new KeplerianOrbit(stateAfterBurn.getOrbit());
            KeplerianOrbit finalOrbit = new KeplerianOrbit(
                    orbitAfter.getA(),
                    orbitAfter.getE(),
                    orbitAfter.getI(),
                    orbitAfter.getPerigeeArgument(),
                    orbitAfter.getRightAscensionOfAscendingNode(),
                    orbitAfter.getMeanAnomaly(),
                    PositionAngle.MEAN,
                    eme2000,
                    stateAfterBurn.getDate(),
                    MU
            );
            this.finalState = new SpacecraftState(finalOrbit, massAfterBurn);

            double fuelUsed = initialMass - massAfterBurn;
            System.out.println("Apogee burn => dv= " + dvApogee
                    + " m/s, fuel= " + fuelUsed + " kg, final inc= "
                    + FastMath.toDegrees(finalOrbit.getI()) + "°");

            Map<String, String> secondPayload = Utils.createDatePayload(stateAfterBurn.getDate(), endHorizonDate);
            Utils.writeJsonPayload(secondPayload);
            Utils.logResults(resultFileName, finalState, finalOrbit,
                    massAfterBurn, DRYMASS);
            writeManeuverTimestamp(lastManeuverDateFile, stateAfterBurn.getDate());
            mqttService.sendFileAndWaitForTrigger(resultFileName, publishTopic, triggerTopic);

            return;
        }

        // Case B: near eq + low e => immediate single-burn
        else if (currentInc <= 0.1) {
            System.out.println("Case B: near-eq + low e => immediate single-burn.");

            // If the difference is negligible, do nothing
            double deltaInc = INC_2 - currentInc;
            if (Math.abs(deltaInc) < 1e-4) {
                System.out.println("No significant inc change => no burn.");
                // Just log results & return
                Utils.logResults(resultFileName, initialState, initialOrbit, initialMass, DRYMASS);
                return;
            }

            // Build local frame for immediate burn
            Vector3D velocity = initialState.getPVCoordinates().getVelocity();
            double speed = velocity.getNorm();

            Vector3D position = initialState.getPVCoordinates().getPosition();
            Vector3D radialUnit = position.normalize();
            Vector3D northUnit = new Vector3D(0, 0, 1);
            Vector3D eastUnit = northUnit.crossProduct(radialUnit).normalize();

            double signNorth = FastMath.signum(velocity.getZ());
            double vHorizFinal = speed * FastMath.cos(INC_2);
            double vVertFinal = speed * FastMath.sin(INC_2) * signNorth;
            Vector3D vDesired = eastUnit.scalarMultiply(vHorizFinal)
                    .add(northUnit.scalarMultiply(vVertFinal));
            Vector3D deltaVVector = vDesired.subtract(velocity);

            AbsoluteDate burnDate = initialDate.shiftedBy(0.001);
            DateDetector burnTrigger = new DateDetector(burnDate);
            ImpulseManeuver singleBurn = new ImpulseManeuver(burnTrigger, deltaVVector, ISP);

            KeplerianPropagator manProp = new KeplerianPropagator(initialOrbit);
            manProp.addEventDetector(singleBurn);
            SpacecraftState stateAfterBurn = manProp.propagate(burnDate);

            double dvImmediate = deltaVVector.getNorm();
            double massAfter = calculateFinalMass(initialMass, dvImmediate, ISP, g0);

            KeplerianOrbit orbitPost = new KeplerianOrbit(stateAfterBurn.getOrbit());
            KeplerianOrbit finalOrbit = new KeplerianOrbit(
                    orbitPost.getA(),
                    orbitPost.getE(),
                    orbitPost.getI(),
                    orbitPost.getPerigeeArgument(),
                    orbitPost.getRightAscensionOfAscendingNode(),
                    orbitPost.getMeanAnomaly(),
                    PositionAngle.MEAN,
                    eme2000,
                    stateAfterBurn.getDate(),
                    MU
            );
            this.finalState = new SpacecraftState(finalOrbit, massAfter);

            double fuelUsed = initialMass - massAfter;
            System.out.println("Immediate burn => dv= " + dvImmediate
                    + " m/s, fuel= " + fuelUsed + " kg, final inc= "
                    + FastMath.toDegrees(finalOrbit.getI()) + "°");

            Map<String, String> secondPayload = Utils.createDatePayload(stateAfterBurn.getDate(), endHorizonDate);
            Utils.writeJsonPayload(secondPayload);
            Utils.logResults(resultFileName, finalState, finalOrbit,
                    massAfter, DRYMASS);
            writeManeuverTimestamp(lastManeuverDateFile, stateAfterBurn.getDate());
            mqttService.sendFileAndWaitForTrigger(resultFileName, publishTopic, triggerTopic);

            return;
        }

        // Case C: standard approach => pick node with lower velocity
        else {
            System.out.println("Case C: inc >= threshold => do node-based approach.");

            // 1) ascend node
            KeplerianPropagator ascendProp = new KeplerianPropagator(initialOrbit);
            AbsoluteDate ascNodeDate = findNodeDate(initialState, true);
            SpacecraftState ascNodeState = ascendProp.propagate(ascNodeDate);
            ascNodeState = new SpacecraftState(ascNodeState.getOrbit(), initialMass);
            double speedAsc = ascNodeState.getPVCoordinates().getVelocity().getNorm();

            // 2) descend node
            KeplerianPropagator descendProp = new KeplerianPropagator(initialOrbit);
            AbsoluteDate descNodeDate = findNodeDate(initialState, false);
            SpacecraftState descNodeState = descendProp.propagate(descNodeDate);
            descNodeState = new SpacecraftState(descNodeState.getOrbit(), initialMass);
            double speedDesc = descNodeState.getPVCoordinates().getVelocity().getNorm();

            // choose slower node
            final boolean useAscending = (speedAsc <= speedDesc);
            AbsoluteDate chosenNodeDate = useAscending ? ascNodeDate : descNodeDate;
            SpacecraftState chosenNodeState = useAscending ? ascNodeState : descNodeState;
            double speedNode = useAscending ? speedAsc : speedDesc;

            System.out.println("Ascending node speed= " + speedAsc
                    + " m/s, descending= " + speedDesc + " m/s => picking "
                    + (useAscending ? "ascending" : "descending") + " node, speed= " + speedNode);

            // 3) build out-of-plane burn
            Vector3D nodePos = chosenNodeState.getPVCoordinates().getPosition();
            Vector3D nodeVel = chosenNodeState.getPVCoordinates().getVelocity();
            double signNorth = FastMath.signum(nodeVel.getZ());

            Vector3D radialUnit = nodePos.normalize();
            Vector3D northUnit = new Vector3D(0, 0, 1);
            Vector3D eastUnit = northUnit.crossProduct(radialUnit).normalize();

            double vHorizFinal = speedNode * FastMath.cos(INC_2);
            double vVertFinal = speedNode * FastMath.sin(INC_2) * signNorth;
            Vector3D vDesired = eastUnit.scalarMultiply(vHorizFinal)
                    .add(northUnit.scalarMultiply(vVertFinal));
            Vector3D deltaVVector = vDesired.subtract(nodeVel);

            double burnShift = 0.001;
            AbsoluteDate burnTriggerDate = chosenNodeDate.shiftedBy(burnShift);
            DateDetector burnTrigger = new DateDetector(burnTriggerDate);
            ImpulseManeuver singleBurn = new ImpulseManeuver(burnTrigger, deltaVVector, ISP);

            // 4) propagate
            KeplerianPropagator manPropNode = new KeplerianPropagator(initialOrbit);
            manPropNode.addEventDetector(singleBurn);
            SpacecraftState stateAfterBurn = manPropNode.propagate(burnTriggerDate);

            double dvNode = deltaVVector.getNorm();
            double massAfterNode = calculateFinalMass(initialMass, dvNode, ISP, g0);

            KeplerianOrbit orbitPost = new KeplerianOrbit(stateAfterBurn.getOrbit());
            KeplerianOrbit finalOrbit = new KeplerianOrbit(
                    orbitPost.getA(),
                    orbitPost.getE(),
                    orbitPost.getI(),
                    orbitPost.getPerigeeArgument(),
                    orbitPost.getRightAscensionOfAscendingNode(),
                    orbitPost.getMeanAnomaly(),
                    PositionAngle.MEAN,
                    eme2000,
                    stateAfterBurn.getDate(),
                    MU
            );
            this.finalState = new SpacecraftState(finalOrbit, massAfterNode);

            double fuelUsed = initialMass - massAfterNode;
            System.out.println("Node-based => dv= " + dvNode
                    + " m/s, used= " + fuelUsed
                    + " kg, final inc= " + FastMath.toDegrees(finalOrbit.getI()) + "°");

            Map<String, String> secondPayload = Utils.createDatePayload(stateAfterBurn.getDate(), endHorizonDate);
            Utils.writeJsonPayload(secondPayload);
            Utils.logResults(resultFileName, finalState, finalOrbit,
                    massAfterNode, DRYMASS);
            writeManeuverTimestamp(lastManeuverDateFile, stateAfterBurn.getDate());
            mqttService.sendFileAndWaitForTrigger(resultFileName, publishTopic, triggerTopic);

        }

        double elapsedSec = (System.currentTimeMillis() - startTime) * 1e-3;
        System.out.println("Inclination maneuver completed in " + elapsedSec + " s");
    }

    // =========================
    //  ADAPTED ERGOL CONSUMPTION
    // =========================
    @Override
    public void calculateErgolConsumption() throws IOException {
        // We'll replicate the "Case A/B/C" logic in an approximate form.

        double initialMass = DRYMASS + ERGOL;
        double incEps = 1e-3;     // near-0 inc threshold
        double eccThresh = 0.05;     // "high e" threshold
        double currentInc = this.INC; // from your field
        double deltaInc = FastMath.abs(INC_2 - currentInc);

        double dv_approx;
        System.out.println("Calculating Ergol for inclination change from "
                + FastMath.toDegrees(currentInc) + "° to "
                + FastMath.toDegrees(INC_2) + "°, e=" + ECC);

        if (currentInc <= incEps && ECC >= eccThresh) {
            System.out.println("Ergol: Case A => near eq + high e => approximate single-burn at apogee.");

            // For approximate "apogee" approach, the local velocity at apogee
            // ~ sqrt(MU*(1-e)/(a*(1+e)))
            // or simpler: we do a rough: v_apogee = sqrt(MU/a*(2/(rA) - 1/a)) if we want.
            // But let's do a simpler approach if we haven't a specialized formula:
            // We'll approximate an elliptical apogee speed as
            //   vA ~ sqrt(MU * (1-e)/(a*(1+e)))
            //   => a = SMA, e = ECC
            double a = SMA;
            double e = ECC;
            double vApogee = FastMath.sqrt(MU * (1.0 - e) / (a * (1.0 + e)));

            // single-burn plane change => dv ~ 2*vApogee * sin(deltaInc/2)
            dv_approx = 2.0 * vApogee * FastMath.sin(deltaInc / 2.0);

        } else if (currentInc <= 0.1) {
            System.out.println("Ergol: Case B => near eq + low e => immediate approach => dv=2v sin(...) for speed~ sqrt(MU/a).");

            // approximate speed ~ sqrt(MU / a)
            double speed = FastMath.sqrt(MU / SMA);
            dv_approx = 2.0 * speed * FastMath.sin(deltaInc / 2.0);

        } else {
            System.out.println("Ergol: Case C => node-based approach => pick ascending/descending => approx using orbital speed at that node.");

            // We'll do a simple approach:
            // If e is not extreme, the difference ascending vs. descending might be small.
            // We'll pick "lowest" => effectively an "average" speed ~ sqrt(MU / a) if near circular.
            // Or we can do a minimal approach if e>0 => assume apogee speed if ascending node near apogee
            // or perigee speed if descending node near perigee.
            // For an approximate approach, let's just do
            // dv ~ 2 * v_node * sin(deltaInc/2).
            // v_node ~ sqrt(MU / a*(1 ± e)) etc.
            // We'll pick the "lowest" => v=(some function)...

            double e = ECC;
            double a = SMA;
            // pretend ascending node near apogee => vAsc= ...
            // pretend descending node near perigee => vDesc= ...
            // We'll pick whichever is smaller.
            double vAsc = FastMath.sqrt(MU * (1.0 - e) / (a * (1.0 + e))); // approximate apogee
            double vDesc = FastMath.sqrt(MU * (1.0 + e) / (a * (1.0 - e))); // approximate perigee

            double vNodeUsed = FastMath.min(vAsc, vDesc);
            dv_approx = 2.0 * vNodeUsed * FastMath.sin(deltaInc / 2.0);
        }

        // once we have dv_approx, compute final mass & fuel used
        double finalMass = calculateFinalMass(initialMass, dv_approx, ISP, g0);
        double fuelUsed = initialMass - finalMass;

        System.out.println("Calculated dv_approx= " + dv_approx + " m/s, fuel= " + fuelUsed + " kg");
        Utils.appendLineToFile(consommationErgolsFile, String.valueOf(fuelUsed));
    }

    */
/**
     * Find next apogee from 'state'. Apogee => TA=π.
     *//*

    private AbsoluteDate findNextApogee(SpacecraftState state) {
        KeplerianOrbit kep = new KeplerianOrbit(state.getOrbit());
        double currentTA = kep.getTrueAnomaly();
        double period = kep.getKeplerianPeriod();

        double diff = (FastMath.PI - currentTA + 2 * FastMath.PI) % (2 * FastMath.PI);
        if (diff < 1e-10) {
            diff += 2 * FastMath.PI;
        }
        double fraction = diff / (2 * FastMath.PI);
        return state.getDate().shiftedBy(fraction * period);
    }

    // =========================
    //  ADAPTED REACH ORBIT TIME
    // =========================
    @Override
    public void processReachOrbitTime() throws IOException {
        double startTime = System.currentTimeMillis();
        manager.addProvider(new DirectoryCrawler(orekitData));
        validateOrbitParameters();

        Frame eme2000 = FramesFactory.getEME2000();
        MqttService mqttService = new MqttService();

        // Build times
        AbsoluteDate initialDate = Utils.parseDateFromTimestamp(DATE);
        AbsoluteDate apsideDate = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
        AbsoluteDate endHorizonDate = new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());

        // Construct orbit at apside date
        KeplerianOrbit orbitAtApside = new KeplerianOrbit(
                SMA, ECC, INC, PA, RAAN, ANO, PositionAngle.MEAN,
                eme2000, apsideDate, MU
        );

        double initialMass = DRYMASS + ERGOL;
        validateManeuverParameters(initialDate, apsideDate, initialMass);

        // Propagate to initial date
        KeplerianPropagator prop = new KeplerianPropagator(orbitAtApside);
        SpacecraftState stateAtInitial = prop.propagate(initialDate);
        KeplerianOrbit initialOrbit = new KeplerianOrbit(stateAtInitial.getOrbit());
        SpacecraftState initialState = new SpacecraftState(initialOrbit, initialMass);

        // For reference or partial re-injection
        double originalPA = initialOrbit.getPerigeeArgument();
        double originalRAAN = initialOrbit.getRightAscensionOfAscendingNode();

        // -------------- Main Logic --------------
        double incEps = 1e-3;   // threshold for near-0 inc (radians)
        double eccThreshold = 0.05;   // "high e" threshold
        double currentInc = initialOrbit.getI();

        System.out.println("Current inc= " + FastMath.toDegrees(currentInc)
                + "°, ECC= " + ECC);

        // Case A: near eq + high e => wait for apogee
        if (currentInc <= incEps && ECC >= eccThreshold) {
            System.out.println("Case A: near-eq + high e => wait for apogee.");

            // 1) find next apogee
            AbsoluteDate apogeeDate = findNextApogee(initialState);
            double smallShift = 1.0;
            AbsoluteDate burnDate = apogeeDate.shiftedBy(smallShift);

            // 2) get velocity at apogee
            SpacecraftState stateApogee = propagateToDate(initialOrbit, apogeeDate);
            Vector3D velocityApg = stateApogee.getPVCoordinates().getVelocity();
            double speedApg = velocityApg.getNorm();

            // Build local radial/east/north
            Vector3D positionApg = stateApogee.getPVCoordinates().getPosition();
            Vector3D radialUnit = positionApg.normalize();
            Vector3D northUnit = new Vector3D(0, 0, 1);
            Vector3D eastUnit = northUnit.crossProduct(radialUnit).normalize();
            double signNorth = FastMath.signum(velocityApg.getZ());

            double vHorizFinal = speedApg * FastMath.cos(INC_2);
            double vVertFinal = speedApg * FastMath.sin(INC_2) * signNorth;

            Vector3D vDesired = eastUnit.scalarMultiply(vHorizFinal)
                    .add(northUnit.scalarMultiply(vVertFinal));
            Vector3D deltaVVector = vDesired.subtract(velocityApg);

            // build impulse
            DateDetector burnTrigger = new DateDetector(burnDate);
            ImpulseManeuver singleBurn = new ImpulseManeuver(burnTrigger, deltaVVector, ISP);

            // propagate
            KeplerianPropagator mp = new KeplerianPropagator(initialOrbit);
            mp.addEventDetector(singleBurn);
            SpacecraftState stateAfterBurn = mp.propagate(burnDate);

            double dvApogee = deltaVVector.getNorm();
            double massAfterBurn = calculateFinalMass(initialMass, dvApogee, ISP, g0);

            KeplerianOrbit orbitAfter = new KeplerianOrbit(stateAfterBurn.getOrbit());
            KeplerianOrbit finalOrbit = new KeplerianOrbit(
                    orbitAfter.getA(),
                    orbitAfter.getE(),
                    orbitAfter.getI(),
                    orbitAfter.getPerigeeArgument(),
                    orbitAfter.getRightAscensionOfAscendingNode(),
                    orbitAfter.getMeanAnomaly(),
                    PositionAngle.MEAN,
                    eme2000,
                    stateAfterBurn.getDate(),
                    MU
            );
            this.finalState = new SpacecraftState(finalOrbit, massAfterBurn);

            double fuelUsed = initialMass - massAfterBurn;
            System.out.println("Apogee burn => dv= " + dvApogee
                    + " m/s, fuel= " + fuelUsed + " kg, final inc= "
                    + FastMath.toDegrees(finalOrbit.getI()) + "°");


            Utils.logApsideDate(apsideFile, stateAfterBurn.getDate());

            return;
        }

        // Case B: near eq + low e => immediate single-burn
        else if (currentInc <= 0.1) {
            System.out.println("Case B: near-eq + low e => immediate single-burn.");

            // If the difference is negligible, do nothing
            double deltaInc = INC_2 - currentInc;
            if (Math.abs(deltaInc) < 1e-4) {
                System.out.println("No significant inc change => no burn.");
                // Just log results & return
                Utils.logResults(resultFileName, initialState, initialOrbit, initialMass, DRYMASS);
                return;
            }

            // Build local frame for immediate burn
            Vector3D velocity = initialState.getPVCoordinates().getVelocity();
            double speed = velocity.getNorm();

            Vector3D position = initialState.getPVCoordinates().getPosition();
            Vector3D radialUnit = position.normalize();
            Vector3D northUnit = new Vector3D(0, 0, 1);
            Vector3D eastUnit = northUnit.crossProduct(radialUnit).normalize();

            double signNorth = FastMath.signum(velocity.getZ());
            double vHorizFinal = speed * FastMath.cos(INC_2);
            double vVertFinal = speed * FastMath.sin(INC_2) * signNorth;
            Vector3D vDesired = eastUnit.scalarMultiply(vHorizFinal)
                    .add(northUnit.scalarMultiply(vVertFinal));
            Vector3D deltaVVector = vDesired.subtract(velocity);

            AbsoluteDate burnDate = initialDate.shiftedBy(0.001);
            DateDetector burnTrigger = new DateDetector(burnDate);
            ImpulseManeuver singleBurn = new ImpulseManeuver(burnTrigger, deltaVVector, ISP);

            KeplerianPropagator manProp = new KeplerianPropagator(initialOrbit);
            manProp.addEventDetector(singleBurn);
            SpacecraftState stateAfterBurn = manProp.propagate(burnDate);

            double dvImmediate = deltaVVector.getNorm();
            double massAfter = calculateFinalMass(initialMass, dvImmediate, ISP, g0);

            KeplerianOrbit orbitPost = new KeplerianOrbit(stateAfterBurn.getOrbit());
            KeplerianOrbit finalOrbit = new KeplerianOrbit(
                    orbitPost.getA(),
                    orbitPost.getE(),
                    orbitPost.getI(),
                    orbitPost.getPerigeeArgument(),
                    orbitPost.getRightAscensionOfAscendingNode(),
                    orbitPost.getMeanAnomaly(),
                    PositionAngle.MEAN,
                    eme2000,
                    stateAfterBurn.getDate(),
                    MU
            );
            this.finalState = new SpacecraftState(finalOrbit, massAfter);

            double fuelUsed = initialMass - massAfter;
            System.out.println("Immediate burn => dv= " + dvImmediate
                    + " m/s, fuel= " + fuelUsed + " kg, final inc= "
                    + FastMath.toDegrees(finalOrbit.getI()) + "°");

            Utils.logApsideDate(apsideFile, stateAfterBurn.getDate());

            return;
        }

        // Case C: standard approach => pick node with lower velocity
        else {
            System.out.println("Case C: inc >= threshold => do node-based approach.");

            // 1) ascend node
            KeplerianPropagator ascendProp = new KeplerianPropagator(initialOrbit);
            AbsoluteDate ascNodeDate = findNodeDate(initialState, true);
            SpacecraftState ascNodeState = ascendProp.propagate(ascNodeDate);
            ascNodeState = new SpacecraftState(ascNodeState.getOrbit(), initialMass);
            double speedAsc = ascNodeState.getPVCoordinates().getVelocity().getNorm();

            // 2) descend node
            KeplerianPropagator descendProp = new KeplerianPropagator(initialOrbit);
            AbsoluteDate descNodeDate = findNodeDate(initialState, false);
            SpacecraftState descNodeState = descendProp.propagate(descNodeDate);
            descNodeState = new SpacecraftState(descNodeState.getOrbit(), initialMass);
            double speedDesc = descNodeState.getPVCoordinates().getVelocity().getNorm();

            // choose slower node
            final boolean useAscending = (speedAsc <= speedDesc);
            AbsoluteDate chosenNodeDate = useAscending ? ascNodeDate : descNodeDate;
            SpacecraftState chosenNodeState = useAscending ? ascNodeState : descNodeState;
            double speedNode = useAscending ? speedAsc : speedDesc;

            System.out.println("Ascending node speed= " + speedAsc
                    + " m/s, descending= " + speedDesc + " m/s => picking "
                    + (useAscending ? "ascending" : "descending") + " node, speed= " + speedNode);

            // 3) build out-of-plane burn
            Vector3D nodePos = chosenNodeState.getPVCoordinates().getPosition();
            Vector3D nodeVel = chosenNodeState.getPVCoordinates().getVelocity();
            double signNorth = FastMath.signum(nodeVel.getZ());

            Vector3D radialUnit = nodePos.normalize();
            Vector3D northUnit = new Vector3D(0, 0, 1);
            Vector3D eastUnit = northUnit.crossProduct(radialUnit).normalize();

            double vHorizFinal = speedNode * FastMath.cos(INC_2);
            double vVertFinal = speedNode * FastMath.sin(INC_2) * signNorth;
            Vector3D vDesired = eastUnit.scalarMultiply(vHorizFinal)
                    .add(northUnit.scalarMultiply(vVertFinal));
            Vector3D deltaVVector = vDesired.subtract(nodeVel);

            double burnShift = 0.001;
            AbsoluteDate burnTriggerDate = chosenNodeDate.shiftedBy(burnShift);
            DateDetector burnTrigger = new DateDetector(burnTriggerDate);
            ImpulseManeuver singleBurn = new ImpulseManeuver(burnTrigger, deltaVVector, ISP);

            // 4) propagate
            KeplerianPropagator manPropNode = new KeplerianPropagator(initialOrbit);
            manPropNode.addEventDetector(singleBurn);
            SpacecraftState stateAfterBurn = manPropNode.propagate(burnTriggerDate);

            double dvNode = deltaVVector.getNorm();
            double massAfterNode = calculateFinalMass(initialMass, dvNode, ISP, g0);

            KeplerianOrbit orbitPost = new KeplerianOrbit(stateAfterBurn.getOrbit());
            KeplerianOrbit finalOrbit = new KeplerianOrbit(
                    orbitPost.getA(),
                    orbitPost.getE(),
                    orbitPost.getI(),
                    orbitPost.getPerigeeArgument(),
                    orbitPost.getRightAscensionOfAscendingNode(),
                    orbitPost.getMeanAnomaly(),
                    PositionAngle.MEAN,
                    eme2000,
                    stateAfterBurn.getDate(),
                    MU
            );
            this.finalState = new SpacecraftState(finalOrbit, massAfterNode);

            double fuelUsed = initialMass - massAfterNode;
            System.out.println("Node-based => dv= " + dvNode
                    + " m/s, used= " + fuelUsed
                    + " kg, final inc= " + FastMath.toDegrees(finalOrbit.getI()) + "°");

            Utils.logApsideDate(apsideFile, stateAfterBurn.getDate());

        }

        double elapsedSec = (System.currentTimeMillis() - startTime) * 1e-3;
        System.out.println("Inclination maneuver completed in " + elapsedSec + " s");
    }
}
*/
