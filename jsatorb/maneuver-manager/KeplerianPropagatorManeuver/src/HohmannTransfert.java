import org.eclipse.paho.client.mqttv3.MqttClient;
import org.eclipse.paho.client.mqttv3.MqttException;
import org.eclipse.paho.client.mqttv3.MqttMessage;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
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
import java.text.ParseException;
import java.util.List;

public class HohmannTransfert {

    private static final double MU = 3.986004418E14;  // Earth's gravitational parameter (m^3/s^2)
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
    static {
        try {
            DATE = Files.readAllLines(Paths.get("Data.txt")).get(1);
            try {
                List<String> lines = Files.readAllLines(Paths.get("LastManeuverDate.txt"));
                if (!lines.isEmpty()) {
                    APSIDE_DATE= lines.get(lines.size() - 1);
                    System.out.println("Last line: " + APSIDE_DATE);
                } else {
                    System.out.println("File is empty.");
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            SMA = Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(2)) * 1000.0;
            ECC = Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(3));
            INC = FastMath.toRadians(Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(4)));
            RAAN = FastMath.toRadians(Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(5)));
            PA = FastMath.toRadians(Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(6)));
            ANO = FastMath.toRadians(Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(7)));
            DRYMASS = Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(8));
            ISP = Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(11));
            ERGOL = Double.parseDouble(Files.readAllLines(Paths.get("Data.txt")).get(12));
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

    // Method to send a file via MQTT
    private static void sendFileViaMQTT(String filePath) {
        String broker = "tcp://mosquitto:1883";  // Adresse du broker MQTT
        String topic = "resultat/fichier";  // Sujet sur lequel publier
        String clientId = "JavaMQTTSender";

        try {
            // Connexion au broker MQTT
            MqttClient client = new MqttClient(broker, clientId);

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
        // if eccentriciy > 5E-3 :
        //Hohmann --> Bielliptic
        // Define frame and initial orbit
        Frame eme2000 = FramesFactory.getEME2000();
        AbsoluteDate dateTLE = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
        AbsoluteDate initialDate = dateTLE.shiftedBy(manoeuverRelativeDate);
        AbsoluteDate apsideDate = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());

        // Define the initial Keplerian orbit
        KeplerianOrbit initialOrbit = new KeplerianOrbit(SMA, ECC, INC, PA, RAAN, ANO, PositionAngle.MEAN, eme2000, initialDate, MU);
        double initialMass = DRYMASS + ERGOL;
        SpacecraftState initialState = new SpacecraftState(initialOrbit, initialMass);

        // Print initial state with logs
        System.out.println("Initial orbit parameters:");
        System.out.println("SMA (Initial semi-major axis): " + initialOrbit.getA() / 1000.0 + " km");
        System.out.println("ECC (Initial eccentricity): " + initialOrbit.getE());
        System.out.println("Mass (Initial mass): " + initialState.getMass() + " kg");
        System.out.println("Initial Date: " + initialDate);
        System.out.println("Apside Date: " + apsideDate);

        // Calculate the Delta-Vs for the Hohmann transfer and log them
        double[] deltaVs = calculateDeltaVs(SMA, SMA_2);
        double DV1 = deltaVs[0];
        double DV2 = deltaVs[1];

        System.out.println("Delta-V1 (First maneuver): " + DV1 + " m/s");
        System.out.println("Delta-V2 (Second maneuver): " + DV2 + " m/s");
        double massIntermediaire =  calculateFinalMass(initialState.getMass() ,DV1 ,ISP ,g0);
        double masseFinale = calculateFinalMass(massIntermediaire ,DV2 ,ISP ,g0);
        //Veerifier si masseFinale >= DRYMASS sinon erreur !
        System.out.println("Initial masse" + initialMass);
        System.out.println("Ergols " + ERGOL);
        System.out.println("ISP" + ISP);
        System.out.println("Masse finale" + masseFinale);
        if (initialDate.isBefore(apsideDate)) {
            System.out.println("initialDate is before apsideDate.");
        } else if (initialDate.isAfter(apsideDate)) {
            throw new IllegalArgumentException("La date initiale est posterieure à la date de l'apside.");
        } else {
            System.out.println("initialDate is equal to apsideDate.");
        }

        double carburant_limite = Math.pow(10, -3);
        if (masseFinale<DRYMASS || ERGOL<carburant_limite) {
            throw new IllegalArgumentException("La masse finale ne peut pas être inférieure à la masse à vide.");
        }
        if (ECC>0.05) {
            throw new IllegalArgumentException("Eccentricity trop importante");
        }
        // Calculate the time to apogee (or perigee) after the first maneuver
        double deltaT = calculateTimeToApogee(SMA, SMA_2);
        System.out.println("Time to Apogee/Perigee (Delta T): " + deltaT + " s");

        // First maneuver (DV1)
        // Appliquer le premier manoeuvre (DV1)
        Vector3D velocityBeforeManeuver = initialState.getPVCoordinates().getVelocity();
        Vector3D directionImpulse1 = velocityBeforeManeuver.normalize().scalarMultiply(DV1);
        Vector3D positionBeforeManeuver = initialState.getPVCoordinates().getPosition();


//        // Use new method to check if we are near apogee or perigee and decide the impulse direction
//        if (isNearApogeeOrPerigee(positionBeforeManeuver, initialOrbit)) {
//            System.out.println("Applying impulse without reversing direction.");
//            directionImpulse1 = velocityBeforeManeuver.normalize().scalarMultiply(DV1);
//        } else {
//            System.out.println("Reversing impulse direction.");
//            directionImpulse1 = velocityBeforeManeuver.normalize().scalarMultiply(DV1).negate();
//        }

        directionImpulse1 = velocityBeforeManeuver.normalize().scalarMultiply(DV1);


        DateDetector maneuverTrigger = new DateDetector(initialDate.shiftedBy(0.001));
        ImpulseManeuver firstManeuver = new ImpulseManeuver(maneuverTrigger, directionImpulse1, ISP);

        // Propagate the orbit with first maneuver
        KeplerianPropagator propagator = new KeplerianPropagator(initialOrbit);
        propagator.addEventDetector(firstManeuver);
        double finalMassAfterFirstManeuver =  calculateFinalMass(initialState.getMass() ,DV1 ,ISP ,g0);
        System.out.println("Final mass after first maneuver: " + finalMassAfterFirstManeuver + " kg");
//        if ((finalMassAfterFirstManeuver - DRYMASS) < 0) {
//            throw new IllegalArgumentException("La masse finale ne peut pas être inférieure à la masse à vide.");
//        }

        SpacecraftState finalState = propagator.propagate(initialDate.shiftedBy(0.001));

        // Second maneuver (DV2) at the correct moment (apogee)
        finalState = propagator.propagate(finalState.getDate().shiftedBy(deltaT));
        //Mettre Apside Detector et update la masse
        //finalState = new SpacecraftState(finalState.getOrbit(), finalMassAfterFirstManeuver);
        AbsoluteDate manoeuverEndDate = finalState.getDate();

        Orbit finalOrbit = new KeplerianOrbit(finalState.getOrbit());

        // Log first maneuver results to Result.txt
        logResults("Result.txt", finalState, finalOrbit, m0, DRYMASS);

        // Send file after first maneuver (intermédiaire)
        AbsoluteDate epoch = new AbsoluteDate(1970, 1, 1, 0, 0, 0.000, TimeScalesFactory.getUTC()); // Unix epoch in AbsoluteDate format

        double durationInSeconds = manoeuverEndDate.durationFrom(epoch);  // Duration from 1st Jan 1970
        long timestampMillis = Math.round(durationInSeconds * 1000);  // Convert seconds to milliseconds
        System.out.println("Time Stamp post-maneuvre :" + manoeuverEndDate);
        System.out.println("Manoeuver end date as timestamp in milliseconds: " + timestampMillis);

        FileWriter writer = new FileWriter("PostManeuverDate.txt", true);
        BufferedWriter bufferedWriter = new BufferedWriter(writer);
        bufferedWriter.newLine();
        bufferedWriter.write("Date post - maneuvre");
        bufferedWriter.newLine();
        bufferedWriter.write(String.valueOf(timestampMillis));
        bufferedWriter.newLine();
        bufferedWriter.close();

        sendFileViaMQTT("Result.txt");

        // Now calculate the second maneuver
        Vector3D velocityBeforeSecondManeuver = finalState.getPVCoordinates().getVelocity();
        Vector3D directionImpulse2 = velocityBeforeSecondManeuver.normalize().scalarMultiply(DV2);
        manoeuverEndDate = finalState.getDate().shiftedBy(0.001);

        DateDetector secondManeuverTrigger = new DateDetector(finalState.getDate().shiftedBy(0.001));
        ImpulseManeuver secondManeuver = new ImpulseManeuver(secondManeuverTrigger, directionImpulse2, ISP);

        // Propagate the orbit after the second maneuver
        propagator = new KeplerianPropagator(finalState.getOrbit());
        propagator.addEventDetector(secondManeuver);
        double finalMassAfterSecondManeuver = calculateFinalMass(finalMassAfterFirstManeuver ,DV2 ,ISP ,g0);

        System.out.println("Final mass after second maneuver: " + finalMassAfterSecondManeuver + " kg");
        finalState = propagator.propagate(secondManeuverTrigger.getDate());
        finalState = new SpacecraftState(finalState.getOrbit(), finalMassAfterSecondManeuver);
        manoeuverEndDate = finalState.getDate();
        finalOrbit = new KeplerianOrbit(finalState.getOrbit());
        logResults("Result.txt", finalState, finalOrbit, m0, DRYMASS);
        epoch = new AbsoluteDate(1970, 1, 1, 0, 0, 0.000, TimeScalesFactory.getUTC()); // Unix epoch in AbsoluteDate format
        durationInSeconds = manoeuverEndDate.durationFrom(epoch);  // Duration from 1st Jan 1970
        timestampMillis = Math.round(durationInSeconds * 1000);  // Convert seconds to milliseconds
        System.out.println("Time Stamp post-maneuvre :" + manoeuverEndDate);
        System.out.println("Manoeuver end date as timestamp in milliseconds: " + timestampMillis);
        writer = new FileWriter("PostManeuverDate.txt", true);
        bufferedWriter = new BufferedWriter(writer);
        bufferedWriter.newLine();
        bufferedWriter.write("Date post - maneuvre");
        bufferedWriter.newLine();
        bufferedWriter.write(String.valueOf(timestampMillis));
        bufferedWriter.newLine();
        bufferedWriter.close();
        writer = new FileWriter("LastManeuverDate.txt", true);
        bufferedWriter = new BufferedWriter(writer);
        bufferedWriter.newLine();
        bufferedWriter.write("Date post - maneuvre");
        bufferedWriter.newLine();
        bufferedWriter.write(String.valueOf(manoeuverEndDate));
        bufferedWriter.newLine();
        bufferedWriter.close();
        // Log final maneuver results to Result.txt


        // Send file after second maneuver (final)
        sendFileViaMQTT("Result.txt");

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

}
