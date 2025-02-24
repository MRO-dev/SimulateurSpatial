import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.json.JSONObject;
import org.orekit.data.DataContext;
import org.orekit.data.DataProvidersManager;
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
import org.orekit.utils.PVCoordinates;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.ParseException;
import java.time.Instant;
import java.time.ZoneId;
import java.time.ZonedDateTime;
import java.util.List;
import java.util.Map;

/**
 * This class implements the Hohmann transfer maneuver strategy.
 */
public class HohmannStrategy implements ManeuverStrategy {

    // Default file names
    private String dataFile = "Data.txt";
    private String maneuvFile = "Maneuv.txt";
    private String lastManeuverDateFile = "LastManeuverDate.txt";
    private String timePersistenceFile = "time-persistence.json";
    private String publishTopic = "resultat/fichier";
    private String triggerTopic = "trigger/continuation";
    private String resultFileName = "Result.txt";
    private String postManeuverDateFileName = "PostManeuverDate.txt";
    private String ergolsFile = "Ergols.txt";
    private String consommationErgolsFile = "ConsommationErgols.txt";
    private String maneuverOrderFile = "ManeuverOrder.txt";
    private String commandDataFile = "CommandData.txt";
    private String timeIntermediateParametersFile = "TimeIntermediateParameters.txt";
    private String apsideFile = "ApsideDates.txt";
    // New variable for the extra parameter
    private String modeParameter = "";
    private String commandParameter = "";
    private static final double MU = 3.986004418E14;
    private static final double TIME_TOLERANCE_SECONDS = 1e-3; // 1 millisecond tolerance// Earth's gravitational parameter (m^3/s^2)
    // Load Orekit data
    static File orekitData = new File("orekit-data");
    static DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
    // Read input data from file
    private String DATE;
    private String APSIDE_DATE;
    private double SMA; // Initial semi-major axis (in meters)
    private double ECC;
    private double INC;
    private double RAAN;
    private double PA;
    private double ANO;
    private double DRYMASS;
    private double ISP;
    private double ERGOL;
    private String MANEUV_TYPE;
    private double SMA_2; // Final semi-major axis (in meters)
    private String endDateString;
    private double m0 = DRYMASS + ERGOL;
    private  double manoeuverRelativeDate;

    // Set constants
    private static double g0 = 9.80665;  // Earth's gravitational acceleration (m/s^2)

    private static int item = 1;


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
            double sma2,
            double isp,
            String dateStr,
            String modeParameter,
            String commandParameter
    ) {
        this.SMA             = sma;
        this.ECC             = ecc;
        this.INC          = incDeg;
        this.RAAN         = raanDeg;
        this.PA          = aopDeg;
        this.ANO     = meanAnomDeg;
        this.DRYMASS         = dryMass;
        this.ERGOL           = ergol;
        this.ISP             = isp;
        this.SMA_2            = sma2;
        this.DATE         = dateStr;
        this.modeParameter   = modeParameter;
        this.commandParameter= commandParameter;

        // If "blue2" => override filenames
        if ("blue2".equalsIgnoreCase(modeParameter)) {
            this.dataFile = "Data2.txt";
            this.maneuvFile = "Maneuv2.txt";
            this.lastManeuverDateFile="LastManeuverDate2.txt";
            this.publishTopic       = "resultat/fichier2";
            this.triggerTopic       = "trigger/continuation2";
            this.resultFileName     = "Result2.txt";
            this.postManeuverDateFileName   = "PostManeuverDate2.txt";
            this.ergolsFile = "Ergols2.txt";
            this.consommationErgolsFile= "ConsommationErgols2.txt";
            this.maneuverOrderFile="ManeuverOrder2.txt";
            this.commandDataFile="CommandData2.txt";
            this.timeIntermediateParametersFile="TimeIntermediateParameters2.txt";
            this.apsideFile="ApsideDates2.txt";
        }

        // If "red1" => override filenames
        if ("red1".equalsIgnoreCase(modeParameter)) {
            this.dataFile = "Data3.txt";
            this.maneuvFile = "Maneuv3.txt";
            this.lastManeuverDateFile="LastManeuverDate3.txt";
            this.publishTopic       = "resultat/fichier3";
            this.triggerTopic       = "trigger/continuation3";
            this.resultFileName     = "Result3.txt";
            this.postManeuverDateFileName   = "PostManeuverDate3.txt";
            this.ergolsFile= "Ergols3.txt";
            this.consommationErgolsFile= "ConsommationErgols3.txt";
            this.maneuverOrderFile="ManeuverOrder3.txt";
            this.commandDataFile="CommandData3.txt";
            this.timeIntermediateParametersFile="TimeIntermediateParameters3.txt";
            this.apsideFile="ApsideDates3.txt";
        }

        // If "red2" => override filenames
        if ("red2".equalsIgnoreCase(modeParameter)) {
            this.dataFile = "Data4.txt";
            this.maneuvFile = "Maneuv4.txt";
            this.lastManeuverDateFile="LastManeuverDate4.txt";
            this.publishTopic       = "resultat/fichier4";
            this.triggerTopic       = "trigger/continuation4";
            this.resultFileName     = "Result4.txt";
            this.postManeuverDateFileName   = "PostManeuverDate4.txt";
            this.ergolsFile= "Ergols4.txt";
            this.consommationErgolsFile= "ConsommationErgols4.txt";
            this.maneuverOrderFile="ManeuverOrder4.txt";
            this.commandDataFile="CommandData4.txt";
            this.timeIntermediateParametersFile="TimeIntermediateParameters4.txt";
            this.apsideFile="ApsideDates4.txt";
        }

        try {
            // Option 1: Override using system properties, if provided
            dataFile = System.getProperty("dataFile", dataFile);
            maneuvFile = System.getProperty("maneuvFile", maneuvFile);
            lastManeuverDateFile = System.getProperty("lastManeuverDateFile", lastManeuverDateFile);
            timePersistenceFile = System.getProperty("timePersistenceFile", timePersistenceFile);
            modeParameter = System.getProperty("modeParameter", "blue1");
            Map<String, String> data = DataParser.parseDataFile(dataFile);
            // Now use the data map instead of reading by index
            DATE = data.get("DATE");
            // Read the LastManeuverDate file
            List<String> allLines = Files.readAllLines(Paths.get(lastManeuverDateFile));
            String extractedApsideDate = null;
            for (int i = allLines.size() - 1; i >= 0; i--) {
                String line = allLines.get(i).trim();
                if (!line.isEmpty() && !line.startsWith("Date post - maneuvre")) {
                    extractedApsideDate = line;
                    break;
                }
            }
            if (extractedApsideDate == null) {
                throw new IOException("No valid apside date found in file.");
            }
            APSIDE_DATE = extractedApsideDate;
            // Read time-persistence JSON file
            String fileContent = new String(Files.readAllBytes(Paths.get(timePersistenceFile)));
            JSONObject jsonObject = new JSONObject(fileContent);
            endDateString = jsonObject.optString("enddate");
            if (endDateString == null || endDateString.isEmpty()) {
                throw new IllegalArgumentException("No valid 'enddate' found in " + timePersistenceFile);
            }
            System.out.println("End date string: " + endDateString);

            // Step 4: Convert the endDateString to an AbsoluteDate
            // (Assuming the date is in a standard ISO-8601 format recognized by Orekit)
            SMA = Double.parseDouble(data.get("SMA")) * 1000.0;
            ECC = Double.parseDouble(data.get("ECC"));
            INC = FastMath.toRadians(Double.parseDouble(data.get("INC")));
            RAAN = FastMath.toRadians(Double.parseDouble(data.get("RAAN")));
            PA = FastMath.toRadians(Double.parseDouble(data.get("AoP")));  // AoP is Argument of Perigee
            ANO = FastMath.toRadians(Double.parseDouble(data.get("MeanAnom")));
            DRYMASS = Double.parseDouble(data.get("Dry Mass"));
            ISP = Double.parseDouble(data.get("ISP"));
            ERGOL = Double.parseDouble(data.get("Ergol mass"));
            MANEUV_TYPE = data.get("ManeuverType");

            // Read the maneuv file for the maneuver relative date
            List<String> maneuvLines = Files.readAllLines(Paths.get(maneuvFile));
            manoeuverRelativeDate = Double.parseDouble(maneuvLines.get(item));
            if ("Hohmann".equalsIgnoreCase(MANEUV_TYPE)) {
                SMA_2 = Double.parseDouble(data.get("SMA_2")) * 1000.0;
                // Proceed with Hohmann transfer calculations
            } else if ("QLaw".equalsIgnoreCase(MANEUV_TYPE)) {
                // Retrieve additional parameters specific to QLaw
                String p1 = data.get("param1");
                String p2 = data.get("param2");
                String p3 = data.get("param3");
                String p4 = data.get("param4");
                // Proceed with QLaw calculations
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
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


    // Example method to compare two AbsoluteDate objects within a tolerance
    private static boolean isEqualOrAfterWithTolerance(AbsoluteDate dateToCheck, AbsoluteDate referenceDate, double toleranceSeconds) {
        // Calculate the difference: positive if dateToCheck is after referenceDate
        double difference = dateToCheck.durationFrom(referenceDate);

        // If difference is greater than or within the negative tolerance, consider it as valid
        return difference >= -toleranceSeconds;
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

    public void printAttributes() {
        System.out.println("========== CLASS ATTRIBUTES ==========");
        System.out.println("Physical Constants:");
        System.out.println(String.format("MU (Earth's gravitational parameter): %.6e m³/s²", MU));
        System.out.println(String.format("g0 (Earth's gravitational acceleration): %.6f m/s²", g0));

        System.out.println("\nFile Paths:");
        System.out.println("Orekit Data Path: " + orekitData.getAbsolutePath());

        System.out.println("\nOrbital Parameters:");
        System.out.println(String.format("Date: %s", DATE));
        System.out.println(String.format("Initial SMA (Semi-Major Axis): %.2f km", SMA/1000.0));
        System.out.println(String.format("Target SMA: %.2f km", SMA_2/1000.0));
        System.out.println(String.format("ECC (Eccentricity): %.6f", ECC));
        System.out.println(String.format("INC (Inclination): %.4f°", FastMath.toDegrees(INC)));
        System.out.println(String.format("RAAN (Right Ascension of Ascending Node): %.4f°", FastMath.toDegrees(RAAN)));
        System.out.println(String.format("PA (Perigee Argument): %.4f°", FastMath.toDegrees(PA)));
        System.out.println(String.format("ANO (Mean Anomaly): %.4f°", FastMath.toDegrees(ANO)));

        System.out.println("\nMass Properties:");
        System.out.println(String.format("Dry Mass: %.2f kg", DRYMASS));
        System.out.println(String.format("Propellant Mass (ERGOL): %.2f kg", ERGOL));
        System.out.println(String.format("Initial Total Mass (m0): %.2f kg", m0));

        System.out.println("\nPropulsion Parameters:");
        System.out.println(String.format("ISP (Specific Impulse): %.2f s", ISP));

        System.out.println("\nManeuver Configuration:");
        System.out.println("Maneuver Type: " + MANEUV_TYPE);
        System.out.println(String.format("Maneuver Relative Date: %.2f s", manoeuverRelativeDate));

        try {
            System.out.println("\nErgols File Data:");
            List<String> ergolsLines = Files.readAllLines(Paths.get(ergolsFile));
            System.out.println("Maneuver Type from Ergols: " + MANEUV_TYPE);
            //System.out.println(String.format("Delta-V from Ergols: %.2f m/s", DV));
        } catch (IOException e) {
            System.out.println("Error reading Ergols.txt: " + e.getMessage());
        }

        System.out.println("====================================");
    }

    /**
     * Logs the apside date to a specified file.
     *
     * @param filePath Path to the log file.
     * @param date     The date of the apside event.
     * @throws IOException If an I/O error occurs.
     */
    private static void logApsideDate(String filePath, AbsoluteDate date) throws IOException {
        try (BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(filePath, true))) {
            bufferedWriter.newLine();
            bufferedWriter.write("Apside reached at:");
            bufferedWriter.newLine();
            bufferedWriter.write(date.toString());
            bufferedWriter.newLine();
        }
    }


    @Override
    public void computeAndExecute() throws IOException, ParseException {
        double start = System.currentTimeMillis();
        manager.addProvider(new DirectoryCrawler(orekitData));
        // if eccentricity > 5E-3 :
        // Hohmann --> Bielliptic
        // Define frame and initial orbit
        Frame eme2000 = FramesFactory.getEME2000();
        MqttService mqttService = new MqttService();
        AbsoluteDate dateTLE = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
        AbsoluteDate initialDate = dateTLE.shiftedBy(manoeuverRelativeDate);

        System.out.println("Apside Date: " + APSIDE_DATE);
        AbsoluteDate apsideDate = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
        AbsoluteDate endHorizonDate = new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());
        // Print it out for demonstration
        System.out.println("Parsed endDate: " + endHorizonDate);

        // 1) Define an orbit at the apside date
        KeplerianOrbit orbitAtApside = new KeplerianOrbit(
                SMA, ECC, INC, PA, RAAN, ANO,
                PositionAngle.MEAN, eme2000, apsideDate, MU
        );
        double initialMass = DRYMASS + ERGOL;
        SpacecraftState stateAtApside = new SpacecraftState(orbitAtApside, initialMass);

        // --- Print logs for orbit at apside date (SMA and Mean Anomaly) ---
        System.out.println("Orbit at Apside Date:");
        System.out.printf("  -> SMA = %.3f km%n", orbitAtApside.getA() / 1000.0);
        System.out.printf("  -> Mean Anomaly = %.3f deg%n",
                FastMath.toDegrees(orbitAtApside.getMeanAnomaly()));

        // 2) Create a KeplerianPropagator with that orbit and propagate from apsideDate -> initialDate
        KeplerianPropagator keplerProp = new KeplerianPropagator(orbitAtApside);
        SpacecraftState stateAtInitialDate = keplerProp.propagate(initialDate);

        // 3) Now define the "initial orbit" and state at that propagated date
        KeplerianOrbit initialOrbit = new KeplerianOrbit(stateAtInitialDate.getOrbit());
        SpacecraftState initialState = new SpacecraftState(initialOrbit, initialMass);

        // --- Print logs for orbit at initial date (SMA and Mean Anomaly) ---
        System.out.println("Orbit at Initial Date (after propagation from apside):");
        System.out.printf("  -> SMA = %.3f km%n", initialOrbit.getA() / 1000.0);
        System.out.printf("  -> Mean Anomaly = %.3f deg%n",
                FastMath.toDegrees(initialOrbit.getMeanAnomaly()));

        // Print initial state with logs
        System.out.println("Initial orbit parameters after propagation from apsideDate:");
        System.out.println("SMA (Initial semi-major axis): " + initialOrbit.getA() / 1000.0 + " km");
        System.out.println("ECC (Initial eccentricity): " + initialOrbit.getE());
        System.out.println("Mass (Initial mass): " + initialState.getMass() + " kg");
        System.out.println("Initial Date: " + initialDate);
        System.out.println("Apside Date: " + apsideDate);

        // === The rest of your Hohmann code remains the same ===
        double[] deltaVs = calculateDeltaVs(SMA, SMA_2);
        double DV1 = deltaVs[0];
        double DV2 = deltaVs[1];
        System.out.println("Delta-V1 (First maneuver): " + DV1 + " m/s");
        System.out.println("Delta-V2 (Second maneuver): " + DV2 + " m/s");

        double massIntermediaire = calculateFinalMass(initialState.getMass(), DV1, ISP, g0);
        double masseFinale = calculateFinalMass(massIntermediaire, DV2, ISP, g0);
        System.out.println("Initial mass: " + initialMass);
        System.out.println("Ergols: " + ERGOL);
        System.out.println("ISP: " + ISP);
        System.out.println("Final mass: " + masseFinale);

        // Check date tolerance
        if (!isEqualOrAfterWithTolerance(initialDate, apsideDate, TIME_TOLERANCE_SECONDS)) {
            throw new IllegalArgumentException(
                    "Initial date must be equal to or later than the apside date within the allowed tolerance.");
        } else {
            System.out.println("Initial date is equal to or after the apside date within tolerance.");
        }

        // Check propellant and mass constraints
        double carburant_limite = Math.pow(10, -3);
        if (masseFinale < DRYMASS || ERGOL < carburant_limite) {
            throw new IllegalArgumentException("La masse finale ne peut pas être inférieure à la masse à vide.");
        }
        if (ECC > 0.05) {
            throw new IllegalArgumentException("Eccentricité trop importante.");
        }

        // Calculate the time to apogee (or perigee) for the Hohmann transfer
        double deltaT = calculateTimeToApogee(SMA, SMA_2);
        System.out.println("Time to Apogee/Perigee (Delta T): " + deltaT + " s");

        // --------------------------- First Maneuver (DV1) -------------------------
        Vector3D velocityBeforeManeuver = initialState.getPVCoordinates().getVelocity();
        Vector3D directionImpulse1 = velocityBeforeManeuver.normalize().scalarMultiply(DV1);

        DateDetector maneuverTrigger = new DateDetector(initialDate.shiftedBy(0.001));
        ImpulseManeuver firstManeuver = new ImpulseManeuver(maneuverTrigger, directionImpulse1, ISP);

        // Propagate the orbit with first maneuver
        KeplerianPropagator propagator = new KeplerianPropagator(initialOrbit);
        propagator.addEventDetector(firstManeuver);
        double finalMassAfterFirstManeuver = calculateFinalMass(initialState.getMass(), DV1, ISP, g0);
        System.out.println("Final mass after first maneuver: " + finalMassAfterFirstManeuver + " kg");

        SpacecraftState finalState = propagator.propagate(initialDate.shiftedBy(0.001));
        // Date when the first maneuver finishes + half orbit to next apsis
        AbsoluteDate manoeuverEndDate = finalState.getDate().shiftedBy(deltaT);
        AbsoluteDate firstManeuverEndDate = manoeuverEndDate;
        if (!isEqualOrAfterWithTolerance(endHorizonDate, firstManeuverEndDate, TIME_TOLERANCE_SECONDS)) {
            throw new IllegalArgumentException(
                    "Initial date must be before than the end of horizon end time !");
        } else {
            System.out.println("Initial date is equal to or after the apside date within tolerance.");
        }

        Map<String, String> firstPayload = Utils.createDatePayload(initialDate, firstManeuverEndDate);
        System.out.println("First Maneuver JSON Payload: " + firstPayload);
        // Write to file
        Utils.writeJsonPayload(firstPayload);

        Orbit finalOrbit = new KeplerianOrbit(finalState.getOrbit());
        finalState = new SpacecraftState(finalState.getOrbit(), finalMassAfterFirstManeuver);

        Utils.logResults(resultFileName, finalState, finalOrbit, m0, DRYMASS);

        // Create a timestamp for the completed first maneuver
        AbsoluteDate epoch = new AbsoluteDate(1970, 1, 1, 0, 0, 0.000, TimeScalesFactory.getUTC());
        double durationInSeconds = manoeuverEndDate.durationFrom(epoch);
        long timestampMillis = Math.round(durationInSeconds * 1000);
        System.out.println("Time Stamp post-maneuver (first): " + manoeuverEndDate);
        System.out.println("Manoeuver end date as timestamp in milliseconds (first): " + timestampMillis);

        // Write that date to a text file
        try (FileWriter writer = new FileWriter(postManeuverDateFileName, true);
             BufferedWriter bufferedWriter = new BufferedWriter(writer)) {
            bufferedWriter.newLine();
            bufferedWriter.write("Date post - maneuver");
            bufferedWriter.newLine();
            bufferedWriter.write(String.valueOf(timestampMillis));
            bufferedWriter.newLine();
        }

        // -------------------- WAIT HERE AFTER FIRST MANEUVER --------------------
        // We wait for an MQTT trigger so the next step can start
        mqttService.sendFileAndWaitForTrigger(resultFileName,publishTopic,triggerTopic);

        // ----------------------- Propagate to the second maneuver ----------------
        finalState = propagator.propagate(finalState.getDate().shiftedBy(deltaT));
        manoeuverEndDate = finalState.getDate();
        finalOrbit = new KeplerianOrbit(finalState.getOrbit());

        // --------------------------- Second Maneuver (DV2) ------------------------
        Vector3D velocityBeforeSecondManeuver = finalState.getPVCoordinates().getVelocity();
        Vector3D directionImpulse2 = velocityBeforeSecondManeuver.normalize().scalarMultiply(DV2);

        DateDetector secondManeuverTrigger = new DateDetector(finalState.getDate().shiftedBy(0.001));
        ImpulseManeuver secondManeuver = new ImpulseManeuver(secondManeuverTrigger, directionImpulse2, ISP);

        // Propagate the orbit with second maneuver
        propagator = new KeplerianPropagator(finalState.getOrbit());
        propagator.addEventDetector(secondManeuver);
        double finalMassAfterSecondManeuver = calculateFinalMass(finalMassAfterFirstManeuver, DV2, ISP, g0);

        System.out.println("Final mass after second maneuver: " + finalMassAfterSecondManeuver + " kg");
        finalState = propagator.propagate(secondManeuverTrigger.getDate());
        finalState = new SpacecraftState(finalState.getOrbit(), finalMassAfterSecondManeuver);
        manoeuverEndDate = finalState.getDate();
        finalOrbit = new KeplerianOrbit(finalState.getOrbit());
        AbsoluteDate secondManeuverEndDate = manoeuverEndDate;

        Map<String, String> secondPayload = Utils.createDatePayload(firstManeuverEndDate, endHorizonDate);
        System.out.println("Second Maneuver JSON Payload: " + secondPayload);
        Utils.writeJsonPayload(secondPayload);

        // Log final orbit after second maneuver
        Utils.logResults(resultFileName, finalState, finalOrbit, m0, DRYMASS);

        epoch = new AbsoluteDate(1970, 1, 1, 0, 0, 0.000, TimeScalesFactory.getUTC());
        durationInSeconds = manoeuverEndDate.durationFrom(epoch);
        timestampMillis = Math.round(durationInSeconds * 1000);
        System.out.println("Time Stamp post-maneuver (second): " + manoeuverEndDate);
        System.out.println("Manoeuver end date as timestamp in milliseconds (second): " + timestampMillis);

        try (FileWriter writer = new FileWriter(postManeuverDateFileName, true);
             BufferedWriter bufferedWriter = new BufferedWriter(writer)) {
            bufferedWriter.newLine();
            bufferedWriter.write("Date post - maneuver");
            bufferedWriter.newLine();
            bufferedWriter.write(String.valueOf(timestampMillis));
            bufferedWriter.newLine();
        }

        try (FileWriter writer = new FileWriter(lastManeuverDateFile, true);
             BufferedWriter bufferedWriter = new BufferedWriter(writer)) {
            bufferedWriter.newLine();
            bufferedWriter.write("Date post - maneuver");
            bufferedWriter.newLine();
            bufferedWriter.write(String.valueOf(manoeuverEndDate));
            bufferedWriter.newLine();
        }

        // -------------------- DO NOT WAIT HERE: simply send via MQTT -------------
        mqttService.sendFileViaMQTT(resultFileName, publishTopic);

        System.out.println("Hohmann transfer complete.");

        double end = System.currentTimeMillis();
        double duration = (end - start) / 1000.0;
        System.out.println("Execution time:" + duration);
    }
    @Override
    public void calculateErgolConsumption(){
        System.out.println("Using ergols input file: " + ergolsFile);
        System.out.println("Using ergols consumption output file: " + consommationErgolsFile);
        System.out.println("Mode parameter: " + modeParameter);
        // Chargement des données Orekit
        File orekitData = new File("orekit-data");
        DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
        manager.addProvider(new DirectoryCrawler(orekitData));
        try {
            DRYMASS = Double.parseDouble(Files.readAllLines(Paths.get(ergolsFile)).get(1));
            System.out.println("Masse à vide : " + DRYMASS + " kg");
            ERGOL = Double.parseDouble(Files.readAllLines(Paths.get(ergolsFile)).get(2));
            System.out.println("Masse d'ergols : " + ERGOL + " kg");
            ISP = Double.parseDouble(Files.readAllLines(Paths.get(ergolsFile)).get(3));
            System.out.println("Impulsion spécifique (ISP) : " + ISP + " s");
            SMA = Double.parseDouble(Files.readAllLines(Paths.get(ergolsFile)).get(5)) * 1000.0;
            System.out.println("SMA initial : " + SMA / 1000.0 + " km");
            SMA_2 = Double.parseDouble(Files.readAllLines(Paths.get(ergolsFile)).get(6)) * 1000.0; // SMA final en mètres
            System.out.println("SMA final : " + SMA_2 / 1000.0 + " km");

            // Calcul de la masse initiale
            double initialMass = DRYMASS + ERGOL;

            // Calculer les Delta-Vs pour le transfert de Hohmann
            double[] deltaVs = calculateDeltaVs(SMA, SMA_2);
            double DV1 = deltaVs[0];
            double DV2 = deltaVs[1];

            // Affichage des Delta-V
            System.out.println("Delta-V1 : " + DV1 + " m/s");
            System.out.println("Delta-V2 : " + DV2 + " m/s");

            // Accélération due à la gravité
            double g0 = 9.80665;

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
            FileWriter writer = new FileWriter(consommationErgolsFile, true);
            BufferedWriter bufferedWriter = new BufferedWriter(writer);
            bufferedWriter.write(String.valueOf(totalErgolConsumed));
            bufferedWriter.close();
        } catch (Exception e) {
            // Gestion des exceptions
            System.err.println("Une erreur s'est produite : " + e.getMessage());
            e.printStackTrace();
            try (BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(consommationErgolsFile, true))) {
                bufferedWriter.write(String.valueOf(0.0));
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }

    }

    @Override
    public void processReachOrbitTime() throws IOException {
        try {
            DATE = Files.readAllLines(Paths.get(maneuverOrderFile)).get(1);
            SMA = Double.parseDouble(Files.readAllLines(Paths.get(commandDataFile)).get(2)) * 1000.0;
            ECC = Double.parseDouble(Files.readAllLines(Paths.get(commandDataFile)).get(3));
            INC = FastMath.toRadians(Double.parseDouble(Files.readAllLines(Paths.get(commandDataFile)).get(4)));
            RAAN = FastMath.toRadians(Double.parseDouble(Files.readAllLines(Paths.get(commandDataFile)).get(5)));
            PA = FastMath.toRadians(Double.parseDouble(Files.readAllLines(Paths.get(commandDataFile)).get(6)));
            ANO = FastMath.toRadians(Double.parseDouble(Files.readAllLines(Paths.get(commandDataFile)).get(7)));
            DRYMASS = Double.parseDouble(Files.readAllLines(Paths.get(ergolsFile)).get(1));
            ISP = new Double(1.0);
            ERGOL = Double.parseDouble(Files.readAllLines(Paths.get(ergolsFile)).get(2));
            MANEUV_TYPE = Files.readAllLines(Paths.get(timeIntermediateParametersFile)).get(1);
            SMA_2 = Double.parseDouble(Files.readAllLines(Paths.get(timeIntermediateParametersFile)).get(2)) * 1000.0;
            int item = 1;
            m0 = DRYMASS + ERGOL;
            manoeuverRelativeDate = Double.parseDouble(Files.readAllLines(Paths.get(maneuvFile)).get(item));
            g0 = 9.80665;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        printAttributes();
        double start = (double) System.currentTimeMillis();
        manager.addProvider(new DirectoryCrawler(orekitData));
        // Define frame and initial date
        Frame eme2000 = FramesFactory.getEME2000();

        // Correctly parse the timestamp as UTC
        long timestamp = Long.parseLong(DATE);
        Instant instant = Instant.ofEpochMilli(timestamp);
        ZonedDateTime zdt = instant.atZone(ZoneId.of("UTC"));

        // Extract date and time components
        int year = zdt.getYear();
        int month = zdt.getMonthValue();
        int day = zdt.getDayOfMonth();
        int hour = zdt.getHour();
        int minute = zdt.getMinute();
        double second = zdt.getSecond() + zdt.getNano() / 1e9;

        // Create AbsoluteDate correctly using UTC
        AbsoluteDate initialDate = new AbsoluteDate(year, month, day, hour, minute, second, TimeScalesFactory.getUTC());

        System.out.println("Date: " + initialDate);

        // Shift by manoeuverRelativeDate to get the first maneuver date
        //AbsoluteDate initialDate = dateTLE.shiftedBy(manoeuverRelativeDate);

        try {
            List<String> lines;
            lines = Files.readAllLines(Paths.get(ergolsFile));
            // Define the initial Keplerian orbit with validation
            if (SMA <= 0 || SMA_2 <= 0) {
                throw new IllegalArgumentException("Semi-major axis must be positive");
            }

            KeplerianOrbit initialOrbit = new KeplerianOrbit(
                    SMA, ECC, INC, PA, RAAN, ANO,
                    PositionAngle.MEAN, eme2000, initialDate, MU
            );

            // Validate initial mass
            if (DRYMASS <= 0 || ERGOL < 0) {
                throw new IllegalArgumentException("Invalid mass parameters");
            }
            double initialMass = DRYMASS + ERGOL;
            SpacecraftState initialState = new SpacecraftState(initialOrbit, initialMass);

            // Calculate and validate Delta-Vs
            double[] deltaVs = calculateDeltaVs(SMA, SMA_2);
            double DV1 = deltaVs[0];
            double DV2 = deltaVs[1];

            // First maneuver calculation
            Vector3D velocityBeforeManeuver = initialState.getPVCoordinates().getVelocity();
            Vector3D directionImpulse1 = velocityBeforeManeuver.normalize().scalarMultiply(DV1);

            // Apply first maneuver
            Vector3D newVelocity = velocityBeforeManeuver.add(directionImpulse1);
            double massAfterDV1 = calculateFinalMass(initialMass, DV1, ISP, g0);

            SpacecraftState stateAfterFirstManeuver = new SpacecraftState(
                    new KeplerianOrbit(
                            new PVCoordinates(initialState.getPVCoordinates().getPosition(), newVelocity),
                            eme2000,
                            initialState.getDate(),
                            MU
                    ),
                    massAfterDV1
            );

            // Set up propagator
            KeplerianPropagator propagator = new KeplerianPropagator(stateAfterFirstManeuver.getOrbit());

            // Create ApsideDetector (not strictly necessary for circular final orbit)
            // Instead, schedule the second maneuver based on the transfer time

            // Calculate time to apogee in transfer orbit
            double timeToApogee = calculateTimeToApogee(SMA, SMA_2); // in seconds

            // Determine the date for the second maneuver
            AbsoluteDate secondManeuverDate = stateAfterFirstManeuver.getDate().shiftedBy(timeToApogee);

            // Apply second maneuver
            Vector3D velocityAtSecondManeuver = stateAfterFirstManeuver.getPVCoordinates().getVelocity();
            Vector3D directionImpulse2 = velocityAtSecondManeuver.normalize().scalarMultiply(DV2);
            Vector3D newVelocityAtSecondManeuver = velocityAtSecondManeuver.add(directionImpulse2);
            double massAfterDV2 = calculateFinalMass(massAfterDV1, DV2, ISP, g0);

            // Create final state after second maneuver
            SpacecraftState finalState = new SpacecraftState(
                    new KeplerianOrbit(
                            new PVCoordinates(stateAfterFirstManeuver.getPVCoordinates().getPosition(), newVelocityAtSecondManeuver),
                            eme2000,
                            secondManeuverDate,
                            MU
                    ),
                    massAfterDV2
            );

            // Log results with the date of the second maneuver
            //logResults("Result.txt", finalState, finalState.getOrbit(), initialMass, DRYMASS);
            logApsideDate(apsideFile, secondManeuverDate); // Log the exact maneuver date as final orbit apogee/perigee

        } catch (Exception e) {
            System.err.println("Error in Hohmann transfer computation: " + e.getMessage());
            throw e;
        }
        double end = (double) System.currentTimeMillis();
        double duration = (end - start) / 1000.0;
        System.out.println("Execution time:" + duration);

    }
}
