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
 * This class implements the Hohmann transfer maneuver strategy.
 */
public class HohmannStrategy extends AbstractManeuverStrategy {
    private double SMA_2;
    // =============== Constructor ==================
    public HohmannStrategy(
            double sma,
            double ecc,
            double incDeg,
            double raanDeg,
            double aopDeg,
            double meanAnomDeg,
            double dryMass,
            double ergol,
            double isp,
            double sma2,
            String dateStr,
            String modeParameter,
            String commandParameter
    ) throws IOException {
        super(modeParameter, commandParameter);
        this.SMA             = sma;
        this.ECC             = ecc;
        this.INC = FastMath.toRadians(incDeg);  // Convert to radians
        this.RAAN = FastMath.toRadians(raanDeg);
        this.PA = FastMath.toRadians(aopDeg);
        this.ANO = FastMath.toRadians(meanAnomDeg);
        this.DRYMASS         = dryMass;
        this.ERGOL           = ergol;
        this.ISP             = isp;
        this.SMA_2            = sma2;
        this.DATE         = dateStr;
        this.modeParameter   = modeParameter;
        this.commandParameter= commandParameter;
    }

    public static double[] calculateDeltaVs(double SMA1, double SMA2) {
        double mu = MU;
        double r1 = SMA1;
        double r2 = SMA2;

        // Vitesses orbitales circulaires initiale et finale
        double v_circ1 = Math.sqrt(mu / r1);
        double v_circ2 = Math.sqrt(mu / r2);

        // Semi-grand axe de l'orbite de transfert
        double a_transfer = (r1 + r2) / 2.0;

        // Vitesses sur l'orbite de transfert au périgée et à l'apogée
        double v_transf_perigee = Math.sqrt(2 * mu * r2 / (r1 * (r1 + r2)));
        double v_transf_apogee = Math.sqrt(2 * mu * r1 / (r2 * (r1 + r2)));

        // Calcul des Delta-V avec les signes corrects
        double deltaV1 = v_transf_perigee - v_circ1;
        double deltaV2 = v_circ2 - v_transf_apogee;

        return new double[]{deltaV1, deltaV2};
    }

    // Méthode pour calculer la masse finale après une manœuvre impulsionnelle
    public static double calculateFinalMass(double initialMass, double deltaV, double ISP, double g0) {
        return initialMass * FastMath.exp(-Math.abs(deltaV)/ (ISP * g0));
    }

    // Method to calculate the time to apogee (or perigee) in a Hohmann transfer
    public static double calculateTimeToApogee(double SMA1_m, double SMA2_m) {
        double a_transfer = (SMA1_m + SMA2_m) / 2.0; // Semi-major axis of the transfer orbit in meters
        double period = 2 * Math.PI * Math.sqrt(Math.pow(a_transfer, 3) / MU); // Orbital period in seconds
        return period / 2.0; // Time to reach apogee (or perigee) in seconds
    }

    private void validateOrbitParameters() {
        if (SMA <= 0 || SMA_2 <= 6600) {
            throw new IllegalArgumentException("Semi-major axis must be positive or superior to 6600 km");
        }
        if (DRYMASS <= 0 || ERGOL < 0) {
            throw new IllegalArgumentException("Invalid mass parameters");
        }
    }

    /**
     * Validates the maneuver parameters for feasibility
     */
    private void validateManeuverParameters(AbsoluteDate initialDate, AbsoluteDate apsideDate, double finalMass) {
        // Check date tolerance
        if (!isEqualOrAfterWithTolerance(initialDate, apsideDate, TIME_TOLERANCE_SECONDS)) {
            throw new IllegalArgumentException(
                    "Initial date must be equal to or later than the apside date within the allowed tolerance.");
        }

        // Check propellant and mass constraints
        double carburant_limite = Math.pow(10, -3);
        System.out.println("Ergols : " + ERGOL);
        System.out.println("DRYLASS : " + DRYMASS);
        System.out.println("Final MASS"+finalMass);
        if (finalMass < DRYMASS || ERGOL < carburant_limite) {
            throw new IllegalArgumentException("La masse finale ne peut pas être inférieure à la masse à vide.");
        }

        if (ECC > 0.05) {
            throw new IllegalArgumentException("Eccentricité trop importante.");
        }
    }

    @Override
    public void loadMassData(Boolean isMassCalculation) throws IOException {
        super.loadMassData(isMassCalculation);
        if (isMassCalculation){
            this.SMA_2   = Double.parseDouble(cmdData.get("SMA_2")) * 1000.0;
        }
        System.out.println("SMA final : " + (SMA_2 / 1000.0) + " km");
    }

    @Override
    public void loadTimeData(Boolean isTimeCalculation) throws IOException {
        super.loadTimeData(isTimeCalculation);
        this.SMA_2 = Double.parseDouble(tipData.get("SMA_2")) * 1000.0;
        System.out.println(" - timeIntermediateParametersFile => ManeuverType: " + MANEUV_TYPE + ", SMA_2: " + SMA_2 + " m");
    }

//    @Override
//    public void computeAndExecute() throws IOException, ParseException {
//        double start = System.currentTimeMillis();
//        manager.addProvider(new DirectoryCrawler(orekitData));
//        // if eccentricity > 5E-3 :
//        // Hohmann --> Bielliptic
//        // Define frame and initial orbit
//        Frame eme2000 = FramesFactory.getEME2000();
//        MqttService mqttService = new MqttService();
//        AbsoluteDate dateTLE = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
//        AbsoluteDate initialDate = dateTLE.shiftedBy(manoeuverRelativeDate);
//
//        System.out.println("Apside Date: " + APSIDE_DATE);
//        AbsoluteDate apsideDate = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
//        AbsoluteDate endHorizonDate = new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());
//        // Print it out for demonstration
//        System.out.println("Parsed endDate: " + endHorizonDate);
//
//        // 1) Define an orbit at the apside date
//        KeplerianOrbit orbitAtApside = new KeplerianOrbit(
//                SMA, ECC, INC, PA, RAAN, ANO,
//                PositionAngle.MEAN, eme2000, apsideDate, MU
//        );
//        double initialMass = DRYMASS + ERGOL;
//        SpacecraftState stateAtApside = new SpacecraftState(orbitAtApside, initialMass);
//
//        // --- Print logs for orbit at apside date (SMA and Mean Anomaly) ---
//        System.out.println("Orbit at Apside Date:");
//        System.out.printf("  -> SMA = %.3f km%n", orbitAtApside.getA() / 1000.0);
//        System.out.printf("  -> Mean Anomaly = %.3f deg%n",
//                FastMath.toDegrees(orbitAtApside.getMeanAnomaly()));
//
//        // 2) Create a KeplerianPropagator with that orbit and propagate from apsideDate -> initialDate
//        KeplerianPropagator keplerProp = new KeplerianPropagator(orbitAtApside);
//        SpacecraftState stateAtInitialDate = keplerProp.propagate(initialDate);
//
//        // 3) Now define the "initial orbit" and state at that propagated date
//        KeplerianOrbit initialOrbit = new KeplerianOrbit(stateAtInitialDate.getOrbit());
//        SpacecraftState initialState = new SpacecraftState(initialOrbit, initialMass);
//
//        // --- Print logs for orbit at initial date (SMA and Mean Anomaly) ---
//        System.out.println("Orbit at Initial Date (after propagation from apside):");
//        System.out.printf("  -> SMA = %.3f km%n", initialOrbit.getA() / 1000.0);
//        System.out.printf("  -> Mean Anomaly = %.3f deg%n",
//                FastMath.toDegrees(initialOrbit.getMeanAnomaly()));
//
//        // Print initial state with logs
//        System.out.println("Initial orbit parameters after propagation from apsideDate:");
//        System.out.println("SMA (Initial semi-major axis): " + initialOrbit.getA() / 1000.0 + " km");
//        System.out.println("ECC (Initial eccentricity): " + initialOrbit.getE());
//        System.out.println("Mass (Initial mass): " + initialState.getMass() + " kg");
//        System.out.println("Initial Date: " + initialDate);
//        System.out.println("Apside Date: " + apsideDate);
//
//        // === The rest of your Hohmann code remains the same ===
//        double[] deltaVs = calculateDeltaVs(SMA, SMA_2);
//        double DV1 = deltaVs[0];
//        double DV2 = deltaVs[1];
//        System.out.println("Delta-V1 (First maneuver): " + DV1 + " m/s");
//        System.out.println("Delta-V2 (Second maneuver): " + DV2 + " m/s");
//
//        double massIntermediaire = calculateFinalMass(initialState.getMass(), DV1, ISP, g0);
//        double masseFinale = calculateFinalMass(massIntermediaire, DV2, ISP, g0);
//        System.out.println("Initial mass: " + initialMass);
//        System.out.println("Ergols: " + ERGOL);
//        System.out.println("ISP: " + ISP);
//        System.out.println("Final mass: " + masseFinale);
//
//        // Check date tolerance
//        if (!isEqualOrAfterWithTolerance(initialDate, apsideDate, TIME_TOLERANCE_SECONDS)) {
//            throw new IllegalArgumentException(
//                    "Initial date must be equal to or later than the apside date within the allowed tolerance.");
//        } else {
//            System.out.println("Initial date is equal to or after the apside date within tolerance.");
//        }
//
//        // Check propellant and mass constraints
//        double carburant_limite = Math.pow(10, -3);
//        if (masseFinale < DRYMASS || ERGOL < carburant_limite) {
//            throw new IllegalArgumentException("La masse finale ne peut pas être inférieure à la masse à vide.");
//        }
//        if (ECC > 0.05) {
//            throw new IllegalArgumentException("Eccentricité trop importante.");
//        }
//
//        // Calculate the time to apogee (or perigee) for the Hohmann transfer
//        double deltaT = calculateTimeToApogee(SMA, SMA_2);
//        System.out.println("Time to Apogee/Perigee (Delta T): " + deltaT + " s");
//
//        // --------------------------- First Maneuver (DV1) -------------------------
//        Vector3D velocityBeforeManeuver = initialState.getPVCoordinates().getVelocity();
//        Vector3D directionImpulse1 = velocityBeforeManeuver.normalize().scalarMultiply(DV1);
//
//        DateDetector maneuverTrigger = new DateDetector(initialDate.shiftedBy(0.001));
//        ImpulseManeuver firstManeuver = new ImpulseManeuver(maneuverTrigger, directionImpulse1, ISP);
//
//        // Propagate the orbit with first maneuver
//        KeplerianPropagator propagator = new KeplerianPropagator(initialOrbit);
//        propagator.addEventDetector(firstManeuver);
//        double finalMassAfterFirstManeuver = calculateFinalMass(initialState.getMass(), DV1, ISP, g0);
//        System.out.println("Final mass after first maneuver: " + finalMassAfterFirstManeuver + " kg");
//
//        SpacecraftState finalState = propagator.propagate(initialDate.shiftedBy(0.001));
//        // Date when the first maneuver finishes + half orbit to next apsis
//        AbsoluteDate manoeuverEndDate = finalState.getDate().shiftedBy(deltaT);
//        AbsoluteDate firstManeuverEndDate = manoeuverEndDate;
//        if (!isEqualOrAfterWithTolerance(endHorizonDate, firstManeuverEndDate, TIME_TOLERANCE_SECONDS)) {
//            throw new IllegalArgumentException(
//                    "Initial date must be before than the end of horizon end time !");
//        } else {
//            System.out.println("Initial date is equal to or after the apside date within tolerance.");
//        }
//
//        Map<String, String> firstPayload = Utils.createDatePayload(initialDate, firstManeuverEndDate);
//        System.out.println("First Maneuver JSON Payload: " + firstPayload);
//        // Write to file
//        Utils.writeJsonPayload(firstPayload);
//
//        Orbit finalOrbit = new KeplerianOrbit(finalState.getOrbit());
//        finalState = new SpacecraftState(finalState.getOrbit(), finalMassAfterFirstManeuver);
//
//        Utils.logResults(resultFileName, finalState, finalOrbit, m0, DRYMASS);
//
//        // Create a timestamp for the completed first maneuver
//        AbsoluteDate epoch = new AbsoluteDate(1970, 1, 1, 0, 0, 0.000, TimeScalesFactory.getUTC());
//        double durationInSeconds = manoeuverEndDate.durationFrom(epoch);
//        long timestampMillis = Math.round(durationInSeconds * 1000);
//        System.out.println("Time Stamp post-maneuver (first): " + manoeuverEndDate);
//        System.out.println("Manoeuver end date as timestamp in milliseconds (first): " + timestampMillis);
//
//        // Write that date to a text file
//        try (FileWriter writer = new FileWriter(postManeuverDateFileName, true);
//             BufferedWriter bufferedWriter = new BufferedWriter(writer)) {
//            bufferedWriter.newLine();
//            bufferedWriter.write("Date post - maneuver");
//            bufferedWriter.newLine();
//            bufferedWriter.write(String.valueOf(timestampMillis));
//            bufferedWriter.newLine();
//        }
//
//        // -------------------- WAIT HERE AFTER FIRST MANEUVER --------------------
//        // We wait for an MQTT trigger so the next step can start
//        mqttService.sendFileAndWaitForTrigger(resultFileName,publishTopic,triggerTopic);
//
//        // ----------------------- Propagate to the second maneuver ----------------
//        finalState = propagator.propagate(finalState.getDate().shiftedBy(deltaT));
//        manoeuverEndDate = finalState.getDate();
//        finalOrbit = new KeplerianOrbit(finalState.getOrbit());
//
//        // --------------------------- Second Maneuver (DV2) ------------------------
//        Vector3D velocityBeforeSecondManeuver = finalState.getPVCoordinates().getVelocity();
//        Vector3D directionImpulse2 = velocityBeforeSecondManeuver.normalize().scalarMultiply(DV2);
//
//        DateDetector secondManeuverTrigger = new DateDetector(finalState.getDate().shiftedBy(0.001));
//        ImpulseManeuver secondManeuver = new ImpulseManeuver(secondManeuverTrigger, directionImpulse2, ISP);
//
//        // Propagate the orbit with second maneuver
//        propagator = new KeplerianPropagator(finalState.getOrbit());
//        propagator.addEventDetector(secondManeuver);
//        double finalMassAfterSecondManeuver = calculateFinalMass(finalMassAfterFirstManeuver, DV2, ISP, g0);
//
//        System.out.println("Final mass after second maneuver: " + finalMassAfterSecondManeuver + " kg");
//        finalState = propagator.propagate(secondManeuverTrigger.getDate());
//        finalState = new SpacecraftState(finalState.getOrbit(), finalMassAfterSecondManeuver);
//        manoeuverEndDate = finalState.getDate();
//        finalOrbit = new KeplerianOrbit(finalState.getOrbit());
//        AbsoluteDate secondManeuverEndDate = manoeuverEndDate;
//
//        Map<String, String> secondPayload = Utils.createDatePayload(firstManeuverEndDate, endHorizonDate);
//        System.out.println("Second Maneuver JSON Payload: " + secondPayload);
//        Utils.writeJsonPayload(secondPayload);
//
//        // Log final orbit after second maneuver
//        Utils.logResults(resultFileName, finalState, finalOrbit, m0, DRYMASS);
//
//        epoch = new AbsoluteDate(1970, 1, 1, 0, 0, 0.000, TimeScalesFactory.getUTC());
//        durationInSeconds = manoeuverEndDate.durationFrom(epoch);
//        timestampMillis = Math.round(durationInSeconds * 1000);
//        System.out.println("Time Stamp post-maneuver (second): " + manoeuverEndDate);
//        System.out.println("Manoeuver end date as timestamp in milliseconds (second): " + timestampMillis);
//
//        try (FileWriter writer = new FileWriter(postManeuverDateFileName, true);
//             BufferedWriter bufferedWriter = new BufferedWriter(writer)) {
//            bufferedWriter.newLine();
//            bufferedWriter.write("Date post - maneuver");
//            bufferedWriter.newLine();
//            bufferedWriter.write(String.valueOf(timestampMillis));
//            bufferedWriter.newLine();
//        }
//
//        try (FileWriter writer = new FileWriter(lastManeuverDateFile, true);
//             BufferedWriter bufferedWriter = new BufferedWriter(writer)) {
//            bufferedWriter.newLine();
//            bufferedWriter.write("Date post - maneuver");
//            bufferedWriter.newLine();
//            bufferedWriter.write(String.valueOf(manoeuverEndDate));
//            bufferedWriter.newLine();
//        }
//
//        // -------------------- DO NOT WAIT HERE: simply send via MQTT -------------
//        mqttService.sendFileViaMQTT(resultFileName, publishTopic);
//
//        System.out.println("Hohmann transfer complete.");
//
//        double end = System.currentTimeMillis();
//        double duration = (end - start) / 1000.0;
//        System.out.println("Execution time:" + duration);
//    }

    @Override
    public void computeAndExecute() throws IOException, ParseException {
        double start = System.currentTimeMillis();
        manager.addProvider(new DirectoryCrawler(orekitData));

        // Setup dates and frames
        Frame eme2000 = FramesFactory.getEME2000();
        MqttService mqttService = new MqttService();
        AbsoluteDate dateTLE = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
        AbsoluteDate initialDate = dateTLE.shiftedBy(manoeuverRelativeDate);
        AbsoluteDate apsideDate = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
        AbsoluteDate endHorizonDate = new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());

        System.out.println("Apside Date: " + APSIDE_DATE);
        System.out.println("Parsed endDate: " + endHorizonDate);

        // Initialize orbit and propagation
        KeplerianOrbit orbitAtApside = new KeplerianOrbit(
                SMA, ECC, INC, PA, RAAN, ANO, PositionAngle.MEAN, eme2000, apsideDate, MU);
        double initialMass = DRYMASS + ERGOL;

        // Propagate from apside to initial date
        KeplerianPropagator keplerProp = new KeplerianPropagator(orbitAtApside);
        SpacecraftState stateAtInitialDate = keplerProp.propagate(initialDate);
        KeplerianOrbit initialOrbit = new KeplerianOrbit(stateAtInitialDate.getOrbit());
        SpacecraftState initialState = new SpacecraftState(initialOrbit, initialMass);

        // Calculate Hohmann transfer parameters
        double[] deltaVs = calculateDeltaVs(SMA, SMA_2);
        double DV1 = deltaVs[0], DV2 = deltaVs[1];
        double massAfterFirstManeuver = calculateFinalMass(initialState.getMass(), DV1, ISP, g0);
        double massAfterSecondManeuver = calculateFinalMass(massAfterFirstManeuver, DV2, ISP, g0);
        double deltaT = calculateTimeToApogee(SMA, SMA_2);

        // Validation checks
        validateManeuverParameters(initialDate, apsideDate, massAfterSecondManeuver);

        // ---------- First Maneuver (DV1) ----------
        Vector3D directionImpulse1 = initialState.getPVCoordinates().getVelocity().normalize().scalarMultiply(DV1);
        DateDetector firstManeuverTrigger = new DateDetector(initialDate.shiftedBy(0.001));
        ImpulseManeuver firstManeuver = new ImpulseManeuver(firstManeuverTrigger, directionImpulse1, ISP);

        KeplerianPropagator propagator = new KeplerianPropagator(initialOrbit);
        propagator.addEventDetector(firstManeuver);

        // Execute first maneuver and propagate
        SpacecraftState finalState = propagator.propagate(initialDate.shiftedBy(0.001));
        AbsoluteDate firstManeuverEndDate = finalState.getDate().shiftedBy(deltaT);

        // Validate end horizon date
        if (!isEqualOrAfterWithTolerance(endHorizonDate, firstManeuverEndDate, TIME_TOLERANCE_SECONDS)) {
            throw new IllegalArgumentException("Initial date must be before the end of horizon time!");
        }

        // Update state with new mass and log results
        Orbit finalOrbit = new KeplerianOrbit(finalState.getOrbit());
        finalState = new SpacecraftState(finalState.getOrbit(), massAfterFirstManeuver);

        // Write date payload and log results
        Map<String, String> firstPayload = Utils.createDatePayload(initialDate, firstManeuverEndDate);
        Utils.writeJsonPayload(firstPayload);
        Utils.logResults(resultFileName, finalState, finalOrbit, m0, DRYMASS);

        // Write post-maneuver timestamp to file
        writeManeuverTimestamp(postManeuverDateFileName, firstManeuverEndDate);

        // Wait for MQTT trigger before second maneuver
        mqttService.sendFileAndWaitForTrigger(resultFileName, publishTopic, triggerTopic);

        // ---------- Second Maneuver (DV2) ----------
        // Propagate to apogee/perigee for second maneuver
        finalState = propagator.propagate(finalState.getDate().shiftedBy(deltaT));
        AbsoluteDate secondManeuverTime = finalState.getDate();

        // Setup second maneuver
        Vector3D directionImpulse2 = finalState.getPVCoordinates().getVelocity().normalize().scalarMultiply(DV2);
        DateDetector secondManeuverTrigger = new DateDetector(secondManeuverTime.shiftedBy(0.001));
        ImpulseManeuver secondManeuver = new ImpulseManeuver(secondManeuverTrigger, directionImpulse2, ISP);

        // Execute second maneuver
        propagator = new KeplerianPropagator(finalState.getOrbit());
        propagator.addEventDetector(secondManeuver);
        finalState = propagator.propagate(secondManeuverTrigger.getDate());
        finalState = new SpacecraftState(finalState.getOrbit(), massAfterSecondManeuver);

        // Final maneuver time and orbit
        AbsoluteDate secondManeuverEndDate = finalState.getDate();
        finalOrbit = new KeplerianOrbit(finalState.getOrbit());

        // Write date payload and log results for second maneuver
        Map<String, String> secondPayload = Utils.createDatePayload(firstManeuverEndDate, endHorizonDate);
        Utils.writeJsonPayload(secondPayload);
        Utils.logResults(resultFileName, finalState, finalOrbit, m0, DRYMASS);

        // Write timestamps to files
        writeManeuverTimestamp(postManeuverDateFileName, secondManeuverEndDate);
        writeManeuverTimestamp(lastManeuverDateFile, secondManeuverEndDate);

        // Send MQTT notification
        mqttService.sendFileViaMQTT(resultFileName, publishTopic);

        System.out.println("Hohmann transfer complete.");
        System.out.println("Execution time:" + (System.currentTimeMillis() - start) / 1000.0);
    }

    @Override
    public void calculateErgolConsumption() throws IOException {
        try {
            // Calcul de la masse initiale
            double initialMass = DRYMASS + ERGOL;
            System.out.println("Initial mass: " + initialMass + " kg");

            // Calculer les Delta-Vs pour le transfert de Hohmann
            double[] deltaVs = calculateDeltaVs(SMA, SMA_2);
            double DV1 = deltaVs[0];
            double DV2 = deltaVs[1];
            // Affichage des Delta-V
            System.out.println("Delta-V1 : " + DV1 + " m/s");
            System.out.println("Delta-V2 : " + DV2 + " m/s");

            // Calcul de la masse après la première manœuvre
            double finalMassAfterFirstManeuver = calculateFinalMass(initialMass, DV1, ISP, g0);
            double ergolConsumedFirstManeuver = initialMass - finalMassAfterFirstManeuver;

            // Calcul de la masse après la deuxième manœuvre
            double finalMassAfterSecondManeuver = calculateFinalMass(finalMassAfterFirstManeuver, DV2, ISP, g0);
            double ergolConsumedSecondManeuver = finalMassAfterFirstManeuver - finalMassAfterSecondManeuver;

            // Calcul total de la consommation d'ergols
            double totalErgolConsumed = ergolConsumedFirstManeuver + ergolConsumedSecondManeuver;

            // Affichage de la consommation d'ergols
            System.out.println("Consommation d'ergols pour la première manœuvre : " + ergolConsumedFirstManeuver + " kg");
            System.out.println("Consommation d'ergols pour la deuxième manœuvre : " + ergolConsumedSecondManeuver + " kg");
            System.out.println("Consommation totale d'ergols : " + totalErgolConsumed + " kg");

            // Enregistrement de la consommation totale dans un fichier
            Utils.appendLineToFile(consommationErgolsFile,String.valueOf(totalErgolConsumed));
//            FileWriter writer = new FileWriter(consommationErgolsFile, true);
//            BufferedWriter bufferedWriter = new BufferedWriter(writer);
//            bufferedWriter.write(String.valueOf(totalErgolConsumed));
//            bufferedWriter.close();
        } catch (Exception e) {
            // Gestion des exceptions
            System.err.println("Une erreur s'est produite : " + e.getMessage());
            e.printStackTrace();
            Utils.appendLineToFile(consommationErgolsFile,String.valueOf(0.0));
        }

    }

    @Override
    public void processReachOrbitTime() throws IOException {
        printAttributes();
        double start = System.currentTimeMillis();
        // 1) Parse the initial date
        AbsoluteDate apsideDate = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
        AbsoluteDate initialDate = Utils.parseDateFromTimestamp(DATE);
        AbsoluteDate endHorizonDate = new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());
        System.out.println("Initial date: " + initialDate);
        try {
            // 2) (Optional) Validate. E.g., if you still want to ensure SMA>0, etc.
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

            // 3) Just compute the time to apogee for Hohmann
            double timeToApogee = calculateTimeToApogee(SMA, SMA_2);
            // 4) The 'second maneuver date' is initialDate + timeToApogee
            AbsoluteDate secondManeuverDate = initialDate.shiftedBy(timeToApogee);
            // 5) Log it
            Utils.logApsideDate(apsideFile, secondManeuverDate);
        } catch (Exception e) {
            System.err.println("Error in Hohmann transfer computation: " + e.getMessage());
            throw e;
        }
        double duration = (System.currentTimeMillis() - start) / 1000.0;
        System.out.println("Execution time: " + duration);
    }
}
