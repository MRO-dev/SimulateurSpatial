//import jdk.jpackage.internal.Log;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.events.Action;
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
import org.orekit.propagation.events.ApsideDetector;
import org.orekit.propagation.events.DateDetector;
import org.orekit.propagation.events.handlers.EventHandler;
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
import java.time.Instant;
import java.time.LocalDateTime;
import java.time.ZoneId;
import java.time.ZonedDateTime;
import java.util.List;

public class TempsIntermediairePredefini {

    private static final double MU = 3.986004418E14;

    // Default file names for various input files
    private static String maneuverOrderFile = "ManeuverOrder.txt";
    private static String commandDataFile = "CommandData.txt";
    private static String ergolsFile = "Ergols.txt";
    private static String timeIntermediateParametersFile = "TimeIntermediateParameters.txt";
    private static String maneuvFile = "Maneuv.txt";
    private static String apsideFile = "ApsideDates.txt";

    // Mode parameter (default is "blue1")
    private static String modeParameter = "blue1";

    static {
        maneuverOrderFile = System.getProperty("maneuverOrderFile", maneuverOrderFile);
        commandDataFile = System.getProperty("commandDataFile", commandDataFile);
        ergolsFile = System.getProperty("ergolsFile", ergolsFile);
        timeIntermediateParametersFile = System.getProperty("timeIntermediateParametersFile", timeIntermediateParametersFile);
        maneuvFile = System.getProperty("maneuvFile", maneuvFile);
        modeParameter = System.getProperty("modeParameter", modeParameter);

        // If modeParameter is "blue2", override the file names with alternative ones.
        if ("blue2".equalsIgnoreCase(modeParameter)) {
            maneuverOrderFile = "ManeuverOrder2.txt";
            commandDataFile = "CommandData2.txt";
            ergolsFile = "Ergols2.txt";
            timeIntermediateParametersFile = "TimeIntermediateParameters2.txt";
            maneuvFile = "Maneuv2.txt";
            apsideFile = "ApsideDates2.txt";
        }

        // If modeParameter is "blue2", override the file names with alternative ones.
        if ("red1".equalsIgnoreCase(modeParameter)) {
            maneuverOrderFile = "ManeuverOrder3.txt";
            commandDataFile = "CommandData3.txt";
            ergolsFile = "Ergols3.txt";
            timeIntermediateParametersFile = "TimeIntermediateParameters3.txt";
            maneuvFile = "Maneuv3.txt";
            apsideFile = "ApsideDates3.txt";
        }

        // If modeParameter is "blue2", override the file names with alternative ones.
        if ("red2".equalsIgnoreCase(modeParameter)) {
            maneuverOrderFile = "ManeuverOrder4.txt";
            commandDataFile = "CommandData4.txt";
            ergolsFile = "Ergols4.txt";
            timeIntermediateParametersFile = "TimeIntermediateParameters4.txt";
            maneuvFile = "Maneuv4.txt";
            apsideFile = "ApsideDates4.txt";
        }
    }// Earth's gravitational parameter (m^3/s^2)
    // Load Orekit data
    static File orekitData = new File("orekit-data");
    static DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
    // Read input data from file
    static String DATE;
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
    static String ManeuverType;
    static double SMA_int = 15000;
    static double SMA_2;
    // Read input data from file
    static List<String> lines;
    static {
        try {
            lines = Files.readAllLines(Paths.get(ergolsFile));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    //static double DV = Double.parseDouble(lines.get(5));  // Delta-V unique
    static {
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
            ManeuverType = MANEUV_TYPE;
            SMA_2 = Double.parseDouble(Files.readAllLines(Paths.get(timeIntermediateParametersFile)).get(2)) * 1000.0;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    static int item = 1;
    static double m0 = DRYMASS + ERGOL;
    static double manoeuverRelativeDate;

    static {
        try {
            manoeuverRelativeDate = Double.parseDouble(Files.readAllLines(Paths.get(maneuvFile)).get(item));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    // Set constants
    static double g0 = 9.80665;  // Earth's gravitational acceleration (m/s^2)


    public TempsIntermediairePredefini() throws IOException {
    }

    public static void printAttributes() {
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
        System.out.println(String.format("Intermediate SMA: %.2f km", SMA_int));
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
            System.out.println("Maneuver Type from Ergols: " + ManeuverType);
            //System.out.println(String.format("Delta-V from Ergols: %.2f m/s", DV));
        } catch (IOException e) {
            System.out.println("Error reading Ergols.txt: " + e.getMessage());
        }

        System.out.println("====================================");
    }


    public static void main(String[] args) throws IOException, ParseException {
        printAttributes();
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



    public static void computeHohmann() throws IOException {
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
    }


    // Helper methods
    private static boolean isAtPerigee(SpacecraftState state) {
        KeplerianOrbit orbit = (KeplerianOrbit) state.getOrbit();
        double tolerance = 1e-6;
        double trueAnomaly = orbit.getTrueAnomaly();
        return Math.abs(trueAnomaly) < tolerance || Math.abs(trueAnomaly - 2 * Math.PI) < tolerance;
    }

    private static boolean isAtApogee(SpacecraftState state) {
        KeplerianOrbit orbit = (KeplerianOrbit) state.getOrbit();
        double tolerance = 1e-6;
        double trueAnomaly = orbit.getTrueAnomaly();
        return Math.abs(trueAnomaly - Math.PI) < tolerance;
    }

    private static void validateFinalOrbit(SpacecraftState finalState, double targetSMA) {
        double tolerance = 1e3; // 1 km tolerance
        if (Math.abs(finalState.getA() - targetSMA) > tolerance) {
            throw new IllegalStateException("Final orbit does not match target parameters");
        }
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

    /**
     * Determines if the current position is near apogee or perigee.
     *
     * @param position Current position vector.
     * @param orbit    Current orbit.
     * @return True if near apogee or perigee, false otherwise.
     */
    public static boolean isNearApogeeOrPerigeeTest(Vector3D position, KeplerianOrbit orbit) {
        double r = position.getNorm();  // Current distance
        double apogee = orbit.getA() * (1 + orbit.getE());  // Apogee distance
        double perigee = orbit.getA() * (1 - orbit.getE()); // Perigee distance

        double tolerance = 1000.0; // Tolerance in meters
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

        SpacecraftState finalState = propagator.propagate(stateAtPerigee.getDate().shiftedBy(0.001));
        finalState = new SpacecraftState(finalState.getOrbit(), finalMassAfterThirdManeuver);

        // Enregistrement des résultats finaux
        //logResults("Result.txt", finalState, finalState.getOrbit(), m0, DRYMASS);

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