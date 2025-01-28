import org.eclipse.paho.client.mqttv3.MqttClient;
import org.eclipse.paho.client.mqttv3.MqttException;
import org.eclipse.paho.client.mqttv3.MqttMessage;
import org.eclipse.paho.client.mqttv3.persist.MemoryPersistence;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
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
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.text.ParseException;
import java.time.format.DateTimeFormatter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CountDownLatch;

public class HohmannTransfert {

    private static final double MU = 3.986004418E14;
    private static final double TIME_TOLERANCE_SECONDS = 1e-3; // 1 millisecond tolerance// Earth's gravitational parameter (m^3/s^2)
    // Load Orekit data
    static File orekitData = new File("orekit-data");
    static DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
    // Read input data from file
    static String DATE;
    static String APSIDE_DATE;
    static double SMA; // Initial semi-major axis (in meters)
    static double ECC;
    static double INC;
    static double RAAN;
    static double PA;
    static double ANO;
    static double DRYMASS;
    static double ISP;
    static double ERGOL;
    static String MANEUV_TYPE;
    static double SMA_int = 15000;
    static double SMA_2; // Final semi-major axis (in meters)
    static String endDateString;
    static {
        try {
            DATE = Files.readAllLines(Paths.get("Data.txt")).get(1);
            // Read all lines from ApsideDates.txt
            List<String> allLines = Files.readAllLines(Paths.get("LastManeuverDate.txt"));

            // Initialize APSIDE_DATE to null or empty
            String extractedApsideDate = null;

            // Iterate over the lines in reverse order
            for (int i = allLines.size() - 1; i >= 0; i--) {
                String line = allLines.get(i).trim();
                // Skip empty lines and lines starting with the prefix
                if (!line.isEmpty() && !line.startsWith("Date post - maneuvre")) {
                    extractedApsideDate = line;
                    break;
                }
            }

            // Validate that a date was found
            if (extractedApsideDate == null) {
                throw new IOException("No valid apside date found in file.");
            }

            // Assign the extracted date to APSIDE_DATE
            APSIDE_DATE = extractedApsideDate;
            String fileContent = new String(Files.readAllBytes(Paths.get("time-persistence.json")));

            // Step 2: Parse it using org.json
            JSONObject jsonObject = new JSONObject(fileContent);

            // Step 3: Retrieve the "enddate" field as a string
            endDateString = jsonObject.optString("enddate");
            if (endDateString == null || endDateString.isEmpty()) {
                throw new IllegalArgumentException("No valid 'enddate' found in time-persistence.json");
            }
            System.out.println("End date string: " + endDateString);
            // Step 4: Convert the endDateString to an AbsoluteDate
            // (Assuming the date is in a standard ISO-8601 format recognized by Orekit)
            SMA = Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(2)) * 1000.0;
            ECC = Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(3));
            INC = FastMath.toRadians(Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(4)));
            RAAN = FastMath.toRadians(Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(5)));
            PA = FastMath.toRadians(Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(6)));
            ANO = FastMath.toRadians(Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(7)));
            DRYMASS = Double.parseDouble(Files.readAllLines(Paths.get("Ergols.txt")).get(1));
            ISP = Double.parseDouble(Files.readAllLines(Paths.get("Ergols.txt")).get(3));
            ERGOL = Double.parseDouble(Files.readAllLines(Paths.get("Ergols.txt")).get(2));
            MANEUV_TYPE = Files.readAllLines(Paths.get("Data.txt")).get(13);
            SMA_2 = Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(14)) * 1000.0;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    static int item = 1;
    static double m0 = DRYMASS + ERGOL;
    static double manoeuverRelativeDate;

    static {
        try {
            manoeuverRelativeDate = Double.parseDouble(Files.readAllLines(Paths.get("Maneuv.txt")).get(item));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    // Set constants
    static double g0 = 9.80665;  // Earth's gravitational acceleration (m/s^2)

    public HohmannTransfert() throws IOException {
    }

    public static void main(String[] args) throws IOException, ParseException {
        double start = (double) System.currentTimeMillis();
        manager.addProvider(new DirectoryCrawler(orekitData));


        switch (MANEUV_TYPE) {
            case "Hohmann":
                System.out.println("HOHMANN");
                computeHohmann();
                break;
            default:
                System.out.println("PARAMETRES INVALIDES !");
                break;
        }
        double end = (double) System.currentTimeMillis();
        double duration = (end - start) / 1000.0;
        System.out.println("Execution time:" + duration);

    }

    // Method to log results in Result.txt
    private static void logResults(String filePath, SpacecraftState finalState, Orbit finalOrbit, double m0, double DRYMASS) throws IOException {
        FileWriter writer = new FileWriter(filePath, true);
        BufferedWriter bufferedWriter = new BufferedWriter(writer);
        bufferedWriter.newLine();
        bufferedWriter.write("Orbital parameters post-maneuver :");
        bufferedWriter.newLine();
        bufferedWriter.write(String.format("%.6f", DRYMASS));
        bufferedWriter.newLine();
        bufferedWriter.write(String.format("%.6f", m0 - finalState.getMass()));
        bufferedWriter.newLine();
        bufferedWriter.write(String.format("%.6f", finalState.getMass() - DRYMASS));
        bufferedWriter.newLine();
        bufferedWriter.write(finalOrbit.getDate().toString());
        bufferedWriter.newLine();
        bufferedWriter.write(String.format("%.6f", (finalState.getPVCoordinates().getPosition().getNorm() - 6378137.0) / 1000.0));
        bufferedWriter.newLine();
        bufferedWriter.write(String.format("%.12f", finalState.getA() / 1000.0));
        System.out.println("SMA (Semi-major axis): " + finalState.getA() / 1000.0 + " km");
        bufferedWriter.newLine();
        bufferedWriter.write(String.format("%.7f", finalState.getE()));
        bufferedWriter.newLine();
        bufferedWriter.write(String.format("%.4f", FastMath.toDegrees(finalState.getI())));
        bufferedWriter.newLine();
        bufferedWriter.write(String.format("%.4f", FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit) finalOrbit).getRightAscensionOfAscendingNode(), Math.PI))));
        bufferedWriter.newLine();
        bufferedWriter.write(String.format("%.8f", FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit) finalOrbit).getPerigeeArgument(), Math.PI))));
        bufferedWriter.newLine();
        bufferedWriter.write(String.format("%.8f", FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit) finalOrbit).getMeanAnomaly(), Math.PI))));
        bufferedWriter.newLine();
        bufferedWriter.write(String.format("%.14f", 86400.0 * finalState.getOrbit().getKeplerianMeanMotion() / 6.283185307179586));
        bufferedWriter.newLine();
        bufferedWriter.write(String.format("%.6f", finalState.getKeplerianPeriod()));
        bufferedWriter.newLine();
        bufferedWriter.close();

    }

    // ==========================================================
    // ==========   MQTT SEND + WAIT-FOR-TRIGGER   =============
    // ==========================================================
    /**
     * Publishes the file path (or its content) and then waits for a trigger message
     * from another MQTT publisher (e.g. Node-RED) on topic "trigger/continuation".
     */
    private static void sendFileAndWaitForTrigger(String filePath) {
        String broker = "tcp://mosquitto:1883";    // The MQTT broker URL
        String publishTopic = "resultat/fichier";  // Topic for sending results
        String triggerTopic = "trigger/continuation"; // Topic to receive the "continue" signal
        String clientId = "JavaMQTTSender_" + System.currentTimeMillis();

        // Latch to block until we receive the "continue" message
        CountDownLatch latch = new CountDownLatch(1);

        try {
            System.out.println("Preparing to send file '" + filePath + "' via MQTT...");
            MqttClient client = new MqttClient(broker, clientId, new MemoryPersistence());


            // Set callback to handle incoming messages
            client.setCallback(new org.eclipse.paho.client.mqttv3.MqttCallback() {
                @Override
                public void connectionLost(Throwable cause) {
                    System.err.println("Connection lost!");
                }

                @Override
                public void messageArrived(String topic, MqttMessage message) throws Exception {
                    System.out.println("Message arrived on topic: " + topic);
                    System.out.println("Payload: " + new String(message.getPayload()));
                    // If this is the trigger message, let's unblock
                    if (topic.equals(triggerTopic)) {
                        latch.countDown();
                    }
                }

                @Override
                public void deliveryComplete(org.eclipse.paho.client.mqttv3.IMqttDeliveryToken token) {
                    // Acknowledgement that the message was delivered
                }
            });

            // Connect and subscribe to the trigger topic before publishing
            client.connect();
            client.subscribe(triggerTopic);

            // Publish the file name or content
            // (You can read the actual file content if you want)
            MqttMessage message = new MqttMessage(("Sending info about file: " + filePath).getBytes());
            message.setQos(2);  // "Exactly once" QoS
            client.publish(publishTopic, message);
            System.out.println("File info published. Now waiting for the trigger message on '" + triggerTopic + "'...");

            // Block until we receive the trigger
            latch.await();

            System.out.println("Trigger message received! Resuming execution...");

            // Disconnect if you want to stop using the client
            client.disconnect();
            client.close();

        } catch (MqttException | InterruptedException e) {
            e.printStackTrace();
        }
    }

    // Method to send a file via MQTT
    private static void sendFileViaMQTT(String filePath) {
        String broker = "tcp://mosquitto:1883";  // Adresse du broker MQTT
        String topic = "resultat/fichier";  // Sujet sur lequel publier
        String clientId = "JavaMQTTSender";

        try {
            // Connexion au broker MQTT
            MqttClient client = new MqttClient(broker, clientId, new MemoryPersistence());


            client.connect();

            // Créer un message MQTT avec le contenu du fichier
            MqttMessage message = new MqttMessage("TEST".getBytes());
            message.setQos(2);  // Niveau de qualité de service, ibufferci "exactement une fois"

            // Publier le message
            client.publish(topic, message);

            System.out.println("Fichier '" + filePath + "' envoyé via MQTT !");
            client.disconnect();

        } catch (MqttException e) {
            e.printStackTrace();
        }
    }

    public static void computeHohmann() throws IOException {
        // if eccentricity > 5E-3 :
        // Hohmann --> Bielliptic
        // Define frame and initial orbit
        Frame eme2000 = FramesFactory.getEME2000();
        AbsoluteDate dateTLE = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
        AbsoluteDate initialDate = dateTLE.shiftedBy(manoeuverRelativeDate);
        System.out.println("Aosude Date: " + APSIDE_DATE);
        AbsoluteDate apsideDate = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
        AbsoluteDate endHorizonDate = new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());
        // Print it out for demonstration
        System.out.println("Parsed endDate: " + endHorizonDate);

        // Define the initial Keplerian orbit
        KeplerianOrbit initialOrbit = new KeplerianOrbit(SMA, ECC, INC, PA, RAAN, ANO,
                PositionAngle.MEAN, eme2000, initialDate, MU);
        double initialMass = DRYMASS + ERGOL;
        SpacecraftState initialState = new SpacecraftState(initialOrbit, initialMass);

        // Print initial state with logs
        System.out.println("Initial orbit parameters:");
        System.out.println("SMA (Initial semi-major axis): " + initialOrbit.getA() / 1000.0 + " km");
        System.out.println("ECC (Initial eccentricity): " + initialOrbit.getE());
        System.out.println("Mass (Initial mass): " + initialState.getMass() + " kg");
        System.out.println("Initial Date: " + initialDate);
        System.out.println("Apside Date: " + apsideDate);

        // Calculate the Delta-Vs for the Hohmann transfer
        double[] deltaVs = calculateDeltaVs(SMA, SMA_2);
        double DV1 = deltaVs[0];
        double DV2 = deltaVs[1];

        System.out.println("Delta-V1 (First maneuver): " + DV1 + " m/s");
        System.out.println("Delta-V2 (Second maneuver): " + DV2 + " m/s");
        double massIntermediaire = calculateFinalMass(initialState.getMass(), DV1, ISP, g0);
        double masseFinale = calculateFinalMass(massIntermediaire, DV2, ISP, g0);
        // Check final mass > DRYMASS
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

        Map<String, String> firstPayload = createDatePayload(initialDate.shiftedBy(180), firstManeuverEndDate.shiftedBy(-90));
        System.out.println("First Maneuver JSON Payload: " + firstPayload);
        // Write to file
        writeJsonPayload(firstPayload);

        Orbit finalOrbit = new KeplerianOrbit(finalState.getOrbit());
        finalState = new SpacecraftState(finalState.getOrbit(), finalMassAfterFirstManeuver);

        logResults("Result.txt", finalState, finalOrbit, m0, DRYMASS);

        // Create a timestamp for the completed first maneuver
        AbsoluteDate epoch = new AbsoluteDate(1970, 1, 1, 0, 0, 0.000, TimeScalesFactory.getUTC());
        double durationInSeconds = manoeuverEndDate.durationFrom(epoch);
        long timestampMillis = Math.round(durationInSeconds * 1000);
        System.out.println("Time Stamp post-maneuver (first): " + manoeuverEndDate);
        System.out.println("Manoeuver end date as timestamp in milliseconds (first): " + timestampMillis);

        // Write that date to a text file
        try (FileWriter writer = new FileWriter("PostManeuverDate.txt", true);
             BufferedWriter bufferedWriter = new BufferedWriter(writer)) {
            bufferedWriter.newLine();
            bufferedWriter.write("Date post - maneuver");
            bufferedWriter.newLine();
            bufferedWriter.write(String.valueOf(timestampMillis));
            bufferedWriter.newLine();
        }

        // -------------------- WAIT HERE AFTER FIRST MANEUVER --------------------
        // We wait for an MQTT trigger so the next step can start
        sendFileAndWaitForTrigger("Result.txt");

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

        Map<String, String> secondPayload = createDatePayload(firstManeuverEndDate.shiftedBy(90), endHorizonDate);
        System.out.println("Second Maneuver JSON Payload: " + secondPayload);
        writeJsonPayload(secondPayload);

        // Log final orbit after second maneuver
        logResults("Result.txt", finalState, finalOrbit, m0, DRYMASS);

        epoch = new AbsoluteDate(1970, 1, 1, 0, 0, 0.000, TimeScalesFactory.getUTC());
        durationInSeconds = manoeuverEndDate.durationFrom(epoch);
        timestampMillis = Math.round(durationInSeconds * 1000);
        System.out.println("Time Stamp post-maneuver (second): " + manoeuverEndDate);
        System.out.println("Manoeuver end date as timestamp in milliseconds (second): " + timestampMillis);

        try (FileWriter writer = new FileWriter("PostManeuverDate.txt", true);
             BufferedWriter bufferedWriter = new BufferedWriter(writer)) {
            bufferedWriter.newLine();
            bufferedWriter.write("Date post - maneuver");
            bufferedWriter.newLine();
            bufferedWriter.write(String.valueOf(timestampMillis));
            bufferedWriter.newLine();
        }

        try (FileWriter writer = new FileWriter("LastManeuverDate.txt", true);
             BufferedWriter bufferedWriter = new BufferedWriter(writer)) {
            bufferedWriter.newLine();
            bufferedWriter.write("Date post - maneuver");
            bufferedWriter.newLine();
            bufferedWriter.write(String.valueOf(manoeuverEndDate));
            bufferedWriter.newLine();
        }

        // -------------------- DO NOT WAIT HERE: simply send via MQTT -------------
        // Instead of waiting for an MQTT trigger, we only publish the result and end
        sendFileViaMQTT("Result.txt");

        System.out.println("Hohmann transfer complete.");
    }



    public static void computeBiElliptic() throws IOException {
        // Charger l'orbite initiale
        Frame eme2000 = FramesFactory.getEME2000();
        AbsoluteDate dateTLE = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
        AbsoluteDate initialDate = dateTLE.shiftedBy(manoeuverRelativeDate);

        // Définir l'orbite képlérienne initiale
        KeplerianOrbit initialOrbit = new KeplerianOrbit(SMA, ECC, INC, PA, RAAN, ANO, PositionAngle.MEAN, eme2000, initialDate, MU);
        double initialMass = DRYMASS + ERGOL;
        SpacecraftState initialState = new SpacecraftState(initialOrbit, initialMass);

        // Afficher les paramètres de l'orbite initiale
        System.out.println("Paramètres de l'orbite initiale :");
        System.out.println("SMA (Axe semi-majeur initial) : " + initialOrbit.getA() / 1000.0 + " km");
        System.out.println("ECC (Excentricité initiale) : " + initialOrbit.getE());
        System.out.println("Masse initiale : " + initialState.getMass() + " kg");
        System.out.println("Date initiale : " + initialDate);

        // Calculer les Delta-V pour le transfert bi-elliptique
        double SMA_intermediate = SMA_int; // Valeur du SMA intermédiaire (à définir ou à calculer)
        double[] deltaVs = calculateBiEllipticDeltaVs(SMA, SMA_2, SMA_intermediate);
        double DV1 = deltaVs[0];
        double DV2 = deltaVs[1];
        double DV3 = deltaVs[2];

        System.out.println("Delta-V1 (Premier manoeuvre) : " + DV1 + " m/s");
        System.out.println("Delta-V2 (Deuxième manoeuvre) : " + DV2 + " m/s");
        System.out.println("Delta-V3 (Troisième manoeuvre) : " + DV3 + " m/s");

        // Appliquer le premier manoeuvre (DV1)
        Vector3D velocityBeforeManeuver = initialState.getPVCoordinates().getVelocity();
        Vector3D directionImpulse1 = velocityBeforeManeuver.normalize().scalarMultiply(DV1);
        ImpulseManeuver firstManeuver = new ImpulseManeuver(new DateDetector(initialDate.shiftedBy(0.001)), directionImpulse1, ISP);

        // Propager l'orbite après le premier manoeuvre
        KeplerianPropagator propagator = new KeplerianPropagator(initialOrbit);
        propagator.addEventDetector(firstManeuver);
        double finalMassAfterFirstManeuver = calculateFinalMass(initialState.getMass() ,DV1 ,ISP ,g0);
        System.out.println("Masse après premier manoeuvre : " + finalMassAfterFirstManeuver + " kg");

        SpacecraftState stateAfterFirstManeuver = propagator.propagate(initialDate.shiftedBy(0.001));
        stateAfterFirstManeuver = new SpacecraftState(stateAfterFirstManeuver.getOrbit(), finalMassAfterFirstManeuver);

        // Appliquer le deuxième manoeuvre (DV2) à l'apogée de l'orbite intermédiaire
        double timeToApogee = calculatePropagationTime(SMA_intermediate);  // Temps pour atteindre l'apogée
        SpacecraftState stateAtApogee = propagator.propagate(stateAfterFirstManeuver.getDate().shiftedBy(timeToApogee));


        Vector3D velocityBeforeSecondManeuver = stateAtApogee.getPVCoordinates().getVelocity();
        Vector3D directionImpulse2 = velocityBeforeSecondManeuver.normalize().scalarMultiply(DV2);
        ImpulseManeuver secondManeuver = new ImpulseManeuver(new DateDetector(stateAtApogee.getDate().shiftedBy(0.001)), directionImpulse2, ISP);

        // Propager l'orbite après le deuxième manoeuvre
        propagator.addEventDetector(secondManeuver);
        double finalMassAfterSecondManeuver = calculateFinalMass(finalMassAfterFirstManeuver ,DV2 ,ISP ,g0);
        System.out.println("Masse après deuxième manoeuvre : " + finalMassAfterSecondManeuver + " kg");

        SpacecraftState stateAfterSecondManeuver = propagator.propagate(stateAtApogee.getDate().shiftedBy(0.001));
        stateAfterSecondManeuver = new SpacecraftState(stateAfterSecondManeuver.getOrbit(), finalMassAfterSecondManeuver);

        // Appliquer le troisième manoeuvre (DV3) pour circulariser l'orbite finale
        double timeToPerigee = calculatePropagationTime(SMA_2);  // Temps pour atteindre le périgée
        SpacecraftState stateAtPerigee = propagator.propagate(stateAfterSecondManeuver.getDate().shiftedBy(timeToPerigee));

        Vector3D velocityBeforeThirdManeuver = stateAtPerigee.getPVCoordinates().getVelocity();
        Vector3D directionImpulse3 = velocityBeforeThirdManeuver.normalize().scalarMultiply(DV3);
        ImpulseManeuver thirdManeuver = new ImpulseManeuver(new DateDetector(stateAtPerigee.getDate().shiftedBy(0.001)), directionImpulse3, ISP);

        // Propager l'orbite finale après le troisième manoeuvre
        propagator.addEventDetector(thirdManeuver);
        double finalMassAfterThirdManeuver = finalMassAfterSecondManeuver * FastMath.exp(-DV3 / (ISP * g0));
        System.out.println("Masse finale après troisième manoeuvre : " + finalMassAfterThirdManeuver + " kg");

        // Calcul de la masse après la première manœuvre
        finalMassAfterFirstManeuver = calculateFinalMass(initialMass, DV1, ISP, g0);
        double ergolConsumedFirstManeuver = initialMass - finalMassAfterFirstManeuver;

        // Calcul de la masse après la deuxième manœuvre
        finalMassAfterSecondManeuver = calculateFinalMass(finalMassAfterFirstManeuver, DV2, ISP, g0);
        double ergolConsumedSecondManeuver = finalMassAfterFirstManeuver - finalMassAfterSecondManeuver;

        // Calcul total de la consommation d'ergols
        double totalErgolConsumed = ergolConsumedFirstManeuver + ergolConsumedSecondManeuver;

        SpacecraftState finalState = propagator.propagate(stateAtPerigee.getDate().shiftedBy(0.001));
        finalState = new SpacecraftState(finalState.getOrbit(), finalMassAfterThirdManeuver);

        // Enregistrement des résultats finaux
        logResults("Result.txt", finalState, finalState.getOrbit(), m0, DRYMASS);

        // Envoyer les résultats via MQTT
        sendFileViaMQTT("Result.txt");

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


    // Méthode pour calculer le temps de propagation à partir du SMA (axe semi-majeur)
    public static double calculatePropagationTime(double SMA) {
        double orbitalPeriod = 2 * Math.PI * Math.sqrt(Math.pow(SMA, 3) / MU);  // Période orbitale en secondes
        return orbitalPeriod / 2.0;  // Temps pour atteindre l'apogée ou le périgée (la moitié de la période)
    }


    // Method to calculate the time to apogee (or perigee) in a Hohmann transfer
    public static double calculateTimeToApogee(double SMA1_m, double SMA2_m) {
        double a_transfer = (SMA1_m + SMA2_m) / 2.0; // Semi-major axis of the transfer orbit in meters
        double period = 2 * Math.PI * Math.sqrt(Math.pow(a_transfer, 3) / MU); // Orbital period in seconds
        return period / 2.0; // Time to reach apogee (or perigee) in seconds
    }


    // New method to check if near apogee or perigee
    public static boolean isNearApogeeOrPerigee(Vector3D position, KeplerianOrbit orbit) {
        double r = position.getNorm();  // Distance actuelle
        double apogee = orbit.getA() * (1 + orbit.getE());  // Distance à l'apogée
        double perigee = orbit.getA() * (1 - orbit.getE());  // Distance au périgée

        double tolerance = 1000.0;  // Tolérance ajustable
        if (Math.abs(r - apogee) < tolerance) {
            System.out.println("Near apogee.");
            return true;
        } else if (Math.abs(r - perigee) < tolerance) {
            System.out.println("Near perigee.");
            return true;
        } else {
            return false;
        }
    }

    public static double[] calculateBiEllipticDeltaVs(double SMA1, double SMA2, double SMA_intermediate) {
        // Premier Delta-V pour passer de SMA1 à l'orbite intermédiaire (SMA_intermediate)
        double DV1 = Math.sqrt(MU / SMA1) * (Math.sqrt(2 * SMA_intermediate / (SMA1 + SMA_intermediate)) - 1);

        // Deuxième Delta-V à l'apogée de l'orbite intermédiaire pour ajuster le périgée à l'orbite cible (SMA2)
        double DV2 = Math.sqrt(MU / SMA_intermediate) * (Math.sqrt(2 * SMA2 / (SMA2 + SMA_intermediate)) - Math.sqrt(2 * SMA1 / (SMA1 + SMA_intermediate)));

        // Troisième Delta-V pour circulariser l'orbite finale
        double DV3 = Math.sqrt(MU / SMA2) * (1 - Math.sqrt(2 * SMA_intermediate / (SMA2 + SMA_intermediate)));

        return new double[]{DV1, DV2, DV3};
    }

    // Méthode pour calculer la masse finale après une manœuvre impulsionnelle
    public static double calculateFinalMass(double initialMass, double deltaV, double ISP, double g0) {
        return initialMass * FastMath.exp(-Math.abs(deltaV)/ (ISP * g0));
    }

    // Example method to compare two AbsoluteDate objects within a tolerance
    private static boolean isEqualOrAfterWithTolerance(AbsoluteDate dateToCheck, AbsoluteDate referenceDate, double toleranceSeconds) {
        // Calculate the difference: positive if dateToCheck is after referenceDate
        double difference = dateToCheck.durationFrom(referenceDate);

        // If difference is greater than or within the negative tolerance, consider it as valid
        return difference >= -toleranceSeconds;
    }

    /**
     * Create a JSON-like payload with start and end dates (ISO format).
     */
    private static Map<String, String> createDatePayload(AbsoluteDate startDate, AbsoluteDate endDate) {
        DateTimeFormatter formatter = DateTimeFormatter.ISO_INSTANT;
        String startIso = formatter.format(startDate.toDate(TimeScalesFactory.getUTC()).toInstant());
        String endIso = formatter.format(endDate.toDate(TimeScalesFactory.getUTC()).toInstant());

        Map<String, String> payload = new HashMap<>();
        payload.put("startdate", startIso);
        payload.put("enddate", endIso);
        return payload;
    }

    private static void writeJsonPayload(Map<String, String> payload) throws IOException {
        // Construct a JSON string from the payload map
        StringBuilder jsonBuilder = new StringBuilder();
        jsonBuilder.append("{\n");
        for (Map.Entry<String, String> entry : payload.entrySet()) {
            jsonBuilder.append("  \"")
                    .append(entry.getKey())
                    .append("\": \"")
                    .append(entry.getValue())
                    .append("\",\n");
        }
        // Remove the last comma and newline, then close the JSON object
        if (!payload.isEmpty()) {
            jsonBuilder.setLength(jsonBuilder.length() - 2);
            jsonBuilder.append("\n");
        }
        jsonBuilder.append("}");

        // Convert to bytes using UTF-8 encoding
        byte[] jsonBytes = jsonBuilder.toString().getBytes(java.nio.charset.StandardCharsets.UTF_8);
        // Write to file, overwriting existing content
        Files.write(Paths.get("time-persistence-intermediate.json"), jsonBytes,
                StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);
    }



}
