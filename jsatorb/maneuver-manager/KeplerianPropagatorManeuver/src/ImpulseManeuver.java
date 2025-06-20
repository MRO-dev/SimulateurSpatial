import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.Locale;
import java.util.TimeZone;

import org.hipparchus.util.FastMath;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.MathUtils;
import org.json.JSONObject;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.attitudes.InertialProvider;
import org.orekit.data.DataContext;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.PVCoordinates;

/*
 * Class to handle impulse maneuvers.
 * Doesn't handle the continuous case.
 */
public class ImpulseManeuver {
    // Default file names (these can be overridden via command-line arguments)
    private static String dataFile               = "Data.txt";
    private static String maneuvFile             = "Maneuv.txt";
    private static String maneuverOrderManualFile = "ManeuverOrderManual.txt";
    private static String lastManeuverDateFile   = "LastManeuverDate.txt";
    private static String timePersistenceFile    = "time-persistence.json";
    private static String ergolsFile             = "Ergols.txt";
    private static String consommationErgolsFile = "ConsommationErgols.txt";
    private static String resultFileName = "Result.txt";
    // Additional mode parameter
    private static String modeParameter = "";
    private static String commandParameter = "";
    private static final double TIME_TOLERANCE_SECONDS = 0.999e-3; // 1 ms tolerance
    static String APSIDE_DATE;
    static String MANEUVER_START_DATE;
    static String endDateString;

    public ImpulseManeuver() {
    }

    /**
     * Reads the file names from command-line arguments (if provided)
     * in the order: dataFile, maneuvFile, lastManeuverDateFile, timePersistenceFile.
     */
    private static void initializeConfiguration(String[] args) {

        try {
            // --- (a) Possibly override default filenames using System properties
            dataFile               = System.getProperty("dataFile", dataFile);
            maneuvFile             = System.getProperty("maneuvFile", maneuvFile);
            lastManeuverDateFile   = System.getProperty("lastManeuverDateFile", lastManeuverDateFile);
            timePersistenceFile    = System.getProperty("timePersistenceFile", timePersistenceFile);
            ergolsFile             = System.getProperty("ergolsFile", ergolsFile);
            consommationErgolsFile = System.getProperty("consommationErgolsFile", consommationErgolsFile);
            resultFileName = System.getProperty("resultFileName", resultFileName);
            maneuverOrderManualFile=System.getProperty("maneuverOrderManualFile", maneuverOrderManualFile);
            modeParameter          = System.getProperty("modeParameter", "blue1");
            commandParameter = System.getProperty("commandParameter", commandParameter);


            // ------------------------------------------------------------------
            // (b) Switch filenames based on modeParameter
            // ------------------------------------------------------------------
            if ("blue2".equalsIgnoreCase(modeParameter)) {
                dataFile               = "Data2.txt";
                maneuvFile             = "Maneuv2.txt";
                lastManeuverDateFile   = "LastManeuverDate2.txt";
                ergolsFile             = "Ergols2.txt";
                consommationErgolsFile = "ConsommationErgols2.txt";
                resultFileName = "Result2.txt";
                maneuverOrderManualFile="ManeuverOrderManual2.txt";
            }
            else if ("red1".equalsIgnoreCase(modeParameter)) {
                dataFile               = "Data3.txt";
                maneuvFile             = "Maneuv3.txt";
                lastManeuverDateFile   = "LastManeuverDate3.txt";
                ergolsFile             = "Ergols3.txt";
                consommationErgolsFile = "ConsommationErgols3.txt";
                resultFileName = "Result3.txt";
                maneuverOrderManualFile="ManeuverOrderManual3.txt";
            }
            else if ("red2".equalsIgnoreCase(modeParameter)) {
                dataFile               = "Data4.txt";
                maneuvFile             = "Maneuv4.txt";
                lastManeuverDateFile   = "LastManeuverDate4.txt";
                ergolsFile             = "Ergols4.txt";
                consommationErgolsFile = "ConsommationErgols4.txt";
                resultFileName = "Result4.txt";
                maneuverOrderManualFile="ManeuverOrderManual4.txt";
            }
        } catch (SecurityException e) {
            System.out.println("SecurityException: " + e.getMessage());
        }
    }

    public static void main(String[] args) throws NullPointerException, ClassCastException, IOException, ParseException {
        initializeConfiguration(args);
        switch (commandParameter) {
            case "executeManeuver":
                computeImpulsManeuver();
                break;
            case "consommationErgols":
                calculateErgolConsumption();
                break;
            default:
                System.out.println("Command parameter not recognized.");
                break;
        }
    }

    public static void computeImpulsManeuver() throws IOException {

        double start = System.currentTimeMillis();

        // -------------------------------------------------------------------
        // 1) Read APSIDE_DATE from the LastManeuverDate file
        // -------------------------------------------------------------------
        List<String> allLines = Files.readAllLines(Paths.get(lastManeuverDateFile));
        String extractedApsideDate = null;
        // Iterate in reverse to find the last non-empty line that isn’t a header
        for (int i = allLines.size() - 1; i >= 0; i--) {
            String line = allLines.get(i).trim();
            if (!line.isEmpty() && !line.startsWith("Apside reached at:")) {
                extractedApsideDate = line;
                break;
            }
        }
        if (extractedApsideDate == null) {
            throw new IOException("No valid apside date found in " + lastManeuverDateFile);
        }
        APSIDE_DATE = extractedApsideDate;
        MANEUVER_START_DATE=Files.readAllLines(Paths.get(maneuverOrderManualFile)).get(0);
        System.out.println("Maneuver start date: " + MANEUVER_START_DATE);

        // -------------------------------------------------------------------
        // 2) Read endDate from the timePersistence file (JSON)
        // -------------------------------------------------------------------
        String fileContent = new String(Files.readAllBytes(Paths.get(timePersistenceFile)));
        JSONObject jsonObject = new JSONObject(fileContent);
        endDateString = jsonObject.optString("enddate");
        if (endDateString == null || endDateString.isEmpty()) {
            throw new IllegalArgumentException("No valid 'enddate' found in " + timePersistenceFile);
        }
        System.out.println("End date string: " + endDateString);

        try {
            List<String> dataLines = Files.readAllLines(Paths.get(dataFile));
            // -------------------------------------------------------------------
            // 3) Parse Data.txt
            // -------------------------------------------------------------------
            // Assuming your file structure (one value per line, starting at index 1)
            // Adjust the indices as needed based on your file’s format.
            int num = 1;
            String DATE = dataLines.get(num);
            num = 2;
            double SMA = Double.parseDouble(dataLines.get(num));
            num = 3;
            double ECC = Double.parseDouble(dataLines.get(num));
            num = 4;
            double INC = Double.parseDouble(dataLines.get(num));
            num = 5;
            double RAAN = Double.parseDouble(dataLines.get(num));
            num = 6;
            double PA = Double.parseDouble(dataLines.get(num));
            num = 7;
            double ANO = Double.parseDouble(dataLines.get(num));

            // Convert to SI and radians
            double sma = SMA * 1000.0;
            double ecc = ECC;
            double inc = FastMath.toRadians(INC);
            double raan = FastMath.toRadians(RAAN);
            double aop = FastMath.toRadians(PA);
            double ano = FastMath.toRadians(ANO);
            double MU = 3.986004418E14;

            // More Data.txt
            Frame GCRF = FramesFactory.getEME2000();
            num = 8;
            double DRYMASS = Double.parseDouble(dataLines.get(num));
            double dryMass = DRYMASS;
            num = 10;
            double THRUST = Double.parseDouble(dataLines.get(num));
            num = 11;
            double ISP = Double.parseDouble(dataLines.get(num));
            num = 12;
            double ERGOL = Double.parseDouble(dataLines.get(num));
            num = 13;
            double DURA = Double.parseDouble(dataLines.get(num));
            num = 14;
            String maneuverType = dataLines.get(num);
            num = 15;
            double DV = Double.parseDouble(dataLines.get(num));
            System.out.println("maneuverType: " + maneuverType);

            double m0 = dryMass + ERGOL;

            try {
                // -------------------------------------------------------------------
                // 4) Parse Maneuv.txt for your maneuver vectors & times
                // -------------------------------------------------------------------
                List<String> maneuvLines = Files.readAllLines(Paths.get(maneuvFile));
                int item = 1;
                double manoeuverRelativeDate = Double.parseDouble(maneuvLines.get(item));
                item = 2;
                double DVx = Double.parseDouble(maneuvLines.get(item));
                item = 3;
                double DVy = Double.parseDouble(maneuvLines.get(item));
                item = 4;
                double DVz = Double.parseDouble(maneuvLines.get(item));
                item = 5;
                double durationOfManoeuver = Double.parseDouble(maneuvLines.get(item));
//                Vector3D direction = new Vector3D(DVx, DVy, DVz);
//                System.out.println("Vecteur direction X= " + direction.getX()
//                        + " Y= " + direction.getY()
//                        + " Z= " + direction.getZ());
                // Normalize the direction vector
                Vector3D directionNormalized = new Vector3D(DVx, DVy, DVz).normalize();
                System.out.println("Normalized Vector direction X= " + directionNormalized.getX()
                        + " Y= " + directionNormalized.getY()
                        + " Z= " + directionNormalized.getZ());

                try {
                    // -------------------------------------------------------------------
                    // 5) Setup Orekit Data
                    // -------------------------------------------------------------------
                    File home = new File("/app/maneuver-manager/");
                    File orekitData = new File(Paths.get("orekit-data").toString());
                    if (!orekitData.exists()) {
                        System.err.format(Locale.US, "Failed to find %s folder%n", orekitData.getAbsolutePath());
                        System.err.format(Locale.US,
                                "You need to download %s from %s, unzip it in %s and rename it 'orekit-data' for this tutorial to work%n",
                                "orekit-data-master.zip",
                                "https://gitlab.orekit.org/orekit/orekit-data/-/archive/master/orekit-data-master.zip",
                                home.getAbsolutePath());
                        System.exit(1);
                    }

                    DataProvidersManager providersManager = DataContext.getDefault().getDataProvidersManager();
                    providersManager.addProvider(new DirectoryCrawler(orekitData));

                    // Build the initial absolute date from TLE
                    AbsoluteDate dateTLE = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
                    System.out.println("DATE :" + DATE);
                    System.out.println("dateTLE :" + dateTLE);

                    // -------------------------------------------------------------------
                    // 6) Build the orbit at APSIDE_DATE
                    // -------------------------------------------------------------------
                    AbsoluteDate apsideDateObj = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
                    Orbit orbitAtApside = new KeplerianOrbit(
                            sma, ecc, inc, aop, raan, ano,
                            PositionAngle.MEAN, GCRF, apsideDateObj, MU
                    );
                    System.out.println("\nOrbit at Apside Date:");
                    System.out.println("  SMA: " + (orbitAtApside.getA() / 1000.0) + " km");
                    System.out.println("  Mean Anomaly: "
                            + FastMath.toDegrees(((KeplerianOrbit) orbitAtApside).getMeanAnomaly())
                            + " deg");

                    // -------------------------------------------------------------------
                    // 7) Propagate from apsideDateObj to manoeuverStartDate
                    // -------------------------------------------------------------------
                    double initMass = m0;
                    double mass_limite = 1e-3; // 0.001

                    KeplerianPropagator keplerProp = new KeplerianPropagator(orbitAtApside);

                    // The actual start date
                    // If MANEUVER_START_DATE is a timestamp in milliseconds
                    long timestampMs = Long.parseLong(MANEUVER_START_DATE);
                    Date date = new Date(timestampMs);
                    SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss.SSS");
                    sdf.setTimeZone(TimeZone.getTimeZone("UTC"));
                    String formattedDate = sdf.format(date);
                    AbsoluteDate manoeuverStartDate = new AbsoluteDate(formattedDate, TimeScalesFactory.getUTC());
                    System.out.println("manoeuverStartDate :" + manoeuverStartDate);
                    System.out.println("apsideDateObj :" + apsideDateObj);

                    // Tolerance check
                    if (!isEqualOrAfterWithTolerance(manoeuverStartDate, apsideDateObj, TIME_TOLERANCE_SECONDS)) {
                        throw new IllegalArgumentException(
                                "Initial date must be >= APSIDE_DATE within tolerance."
                        );
                    } else {
                        System.out.println("Initial date is equal or after the apside date within tolerance.");
                    }

                    // Propagate to that date
                    SpacecraftState stateAtManoeuverStartDate = keplerProp.propagate(manoeuverStartDate);

                    // Rebuild the orbit at that date
                    Orbit iniOrbit = new KeplerianOrbit(stateAtManoeuverStartDate.getOrbit());
                    SpacecraftState initialState = new SpacecraftState(iniOrbit, initMass);

                    // Log
                    System.out.println("\nOrbit at manoeuverStartDate:");
                    System.out.println("  SMA: " + (iniOrbit.getA() / 1000.0) + " km");
                    System.out.println("  Mean Anomaly: "
                            + FastMath.toDegrees(((KeplerianOrbit) iniOrbit).getMeanAnomaly())
                            + " deg");

                    System.out.println("iniOrbit:" + iniOrbit);
                    System.out.println("initMass :" + initMass);

                    // End horizon
                    AbsoluteDate endHorizonDate = new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());
                    AbsoluteDate manoeuverEndDate = manoeuverStartDate.shiftedBy(durationOfManoeuver);

                    // Possibly use the orbitType if needed
                    OrbitType orbitType = OrbitType.CIRCULAR;
                    if (ecc > 0.01) {
                        orbitType = OrbitType.KEPLERIAN;
                    }

                    // Log initial orbit
                    System.out.println("Initial orbit:");
                    System.out.println("Date:" + iniOrbit.getDate());
                    System.out.println("SMA:" + iniOrbit.getA() / 1000.0 + " km");
                    System.out.println("MASS:" + initialState.getMass() + " kg");
                    if (initialState.getMass() < mass_limite) {
                        throw new IllegalArgumentException("La masse initiale ne peut pas être négligeable.");
                    }
                    System.out.println("Mean:" + FastMath.toDegrees(new KeplerianOrbit(iniOrbit).getMeanAnomaly()) + " °");
                    System.out.println("RAAN:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit) iniOrbit).getRightAscensionOfAscendingNode(), Math.PI)));
                    System.out.println("ECC:" + iniOrbit.getE());
                    System.out.println("AOP:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit) iniOrbit).getPerigeeArgument(), Math.PI)));
                    System.out.println("Mean Motion :" + (86400.0 * initialState.getOrbit().getKeplerianMeanMotion() / (2 * Math.PI)));
                    System.out.println("Incli :" + FastMath.toDegrees(initialState.getI()));
                    System.out.println("Fram :" + initialState.getFrame());
                    System.out.println("PV:" + initialState.getOrbit().getPVCoordinates());
                    System.out.println("Kep:" + initialState.getOrbit());

                    // finalState is the post-burn result
                    SpacecraftState finalState = null;

                    // -------------------------------------------------------------------
                    // 8) The Maneuver
                    // -------------------------------------------------------------------
                    SpacecraftState postBurnState = null;
                    if (maneuverType.equals("Impulse")) {
                        // 8A) Pre-burn position/velocity
                        Vector3D posBeforeBurn = initialState.getPVCoordinates().getPosition();
                        Vector3D velBeforeBurn = initialState.getPVCoordinates().getVelocity();

                        // 8B) Build the delta-V vector using the normalized direction
                        Vector3D directionImp = directionNormalized.scalarMultiply(DV);
                        Vector3D velAfterBurn = velBeforeBurn.add(directionImp);

                        // 8C) Validate that the resulting orbit will be elliptic
                        double vAfterBurnMag = velAfterBurn.getNorm();
                        double specificOrbitalEnergy = vAfterBurnMag * vAfterBurnMag / 2 - MU / posBeforeBurn.getNorm();
                        if (specificOrbitalEnergy >= 0) {
                            throw new IllegalArgumentException("Delta-V results in a parabolic or hyperbolic orbit.");
                        }

                        // 8D) Build the post-burn orbit
                        Orbit postBurnOrbit = new KeplerianOrbit(
                                new PVCoordinates(posBeforeBurn, velAfterBurn),
                                iniOrbit.getFrame(),
                                manoeuverStartDate,
                                MU
                        );

                        // 8D) Tsiolkovsky for new mass
                        double g0 = 9.80665;
                        double actualDV = directionImp.getNorm();
                        double finalMass = initialState.getMass() * Math.exp(-actualDV / (ISP * g0));

                        // Debug prints
                        System.out.println("\n--- IMPULSE DEBUG ---");
                        System.out.printf("actualDV: %.3f m/s\n", actualDV);
                        System.out.printf("finalMass (Tsiolkovsky): %.3f kg\n", finalMass);

                        // 8E) Build a new state
                        postBurnState = new SpacecraftState(postBurnOrbit, finalMass);

                        // 8F) Construct a new KeplerianPropagator with inertial attitude
                        AttitudeProvider inertialProvider = new InertialProvider(FramesFactory.getEME2000());
                        KeplerianPropagator postBurnProp = new KeplerianPropagator(
                                postBurnState.getOrbit(),
                                inertialProvider,
                                postBurnState.getMass()
                        );

                        // 8G) Propagate from start to end of the burn
                        finalState = postBurnProp.propagate(manoeuverEndDate);
                        System.out.println("postBurnSate Mass:" + postBurnState.getMass() + " kg");


// Alternatively, check after constructing the orbit
                        if (postBurnOrbit.getE() >= 1.0 && postBurnOrbit.getA() > 0.0) {
                            throw new IllegalArgumentException("Resulting orbit is invalid: e >= 1.0 with a > 0.0");
                        }


                    } else if (maneuverType.equals("Continuous")) {
                        // (unchanged from your original approach)
                        double thrustStep = 1.0;
                        double isp = ISP;
                        double g0 = 9.80665;
                        double thrust = THRUST;

                        double initialMass = initialState.getMass();
                        double currentMass = initialMass;
                        dryMass = DRYMASS;

                        double totalDuration = durationOfManoeuver;
                        SpacecraftState currentState = initialState;

                        while (currentState.getDate().compareTo(manoeuverEndDate) < 0) {
                            double DVxPerStep = (DVx * DV) / totalDuration;
                            double DVyPerStep = (DVy * DV) / totalDuration;
                            double DVzPerStep = (DVz * DV) / totalDuration;

                            Vector3D thrustPerStep = new Vector3D(DVxPerStep, DVyPerStep, DVzPerStep);
                            Vector3D updatedVelocity = currentState.getPVCoordinates().getVelocity().add(thrustPerStep);

                            double ergolConsumedByStep = (thrust * thrustStep) / (isp * g0);
                            currentMass -= ergolConsumedByStep;

                            if (currentMass < dryMass) {
                                currentMass = dryMass;
                                break;
                            }

                            PVCoordinates updatedPV = new PVCoordinates(currentState.getPVCoordinates().getPosition(), updatedVelocity);
                            Orbit updatedOrbit = new KeplerianOrbit(updatedPV, currentState.getFrame(), currentState.getDate(), MU);
                            currentState = new SpacecraftState(updatedOrbit, currentMass);

                            KeplerianPropagator updatedPropagator = new KeplerianPropagator(updatedOrbit);
                            currentState = updatedPropagator.propagate(currentState.getDate().shiftedBy(thrustStep));

                            System.out.println("Masse actuelle après mise à jour : " + currentMass + " kg");
                        }

                        double ergolRestant = currentMass - dryMass;
                        finalState = new SpacecraftState(currentState.getOrbit(), currentMass);

                    } else {
                        // Maneuver: NONE
                        System.out.println("Maneuver : NONE");
                    }

                    // If no maneuver took place at all, finalState = initial
                    if (finalState == null) {
                        finalState = initialState;
                    }

                    // -------------------------------------------------------------------
                    // 9) Analyze final orbit
                    // -------------------------------------------------------------------
                    Orbit finalOrbit = new KeplerianOrbit(finalState.getOrbit());
                    System.out.println("\nAfter the last maneuver / final orbit:");
                    System.out.println("Date: " + finalState.getDate());
                    System.out.println("SMA: " + finalState.getA() / 1000.0 + " km");
                    System.out.println("MeanAnomaly = "
                            + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit) finalOrbit).getMeanAnomaly(), Math.PI))
                            + " deg");
                    System.out.println("manoeuverStartDate:" + manoeuverStartDate + " .");
                    System.out.println("manoeuverEndDate:" + manoeuverEndDate + " .");
                    System.out.println("RAAN:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit) finalOrbit).getRightAscensionOfAscendingNode(), Math.PI)));
                    System.out.println("ECC:" + finalOrbit.getE());
                    System.out.println("AOP:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit) finalOrbit).getPerigeeArgument(), Math.PI)));
                    System.out.println("Mean Motion :" + 86400.0 * finalState.getOrbit().getKeplerianMeanMotion() / (2 * Math.PI));
                    System.out.println("Incli :" + FastMath.toDegrees(finalState.getI()));
                    System.out.println("PV:" + finalState.getOrbit().getPVCoordinates());
                    System.out.println("Kep:" + finalState.getOrbit());
                    System.out.println("Final mass (after maneuver): " + postBurnState.getMass());

                    // -------------------------------------------------------------------
                    // 10) Check that final mass is not below dry mass
                    //     (allow exactly = dryMass, but not below)
                    // -------------------------------------------------------------------
                    double carburant_limite = 1e-3;
                    if (postBurnState.getMass() < (dryMass + carburant_limite)) {
                        throw new IllegalArgumentException("La masse finale ne peut pas être inférieure à la masse à vide.");
                    }

                    // Also check we didn't start with near-zero prop
                    if (ERGOL < carburant_limite) {
                        throw new IllegalArgumentException("ERGOL initial < 0.001 kg? Not valid.");
                    }

                    if (finalState.getE() >= 1.0) {
                        throw new IllegalStateException("hyperbolic orbits cannot be handled");
                    }
                    double n = postBurnState.getOrbit().getKeplerianMeanMotion();
                    System.out.println("n:" + n);
                    if (n < 1e-9) {
                        throw new IllegalStateException("hyperbolic orbits cannot be handled");
                    }



                    // -------------------------------------------------------------------
                    // 11) Write results to file
                    // -------------------------------------------------------------------
                    FileWriter writer = new FileWriter(resultFileName, true);
                    BufferedWriter bufferedWriter = new BufferedWriter(writer);
                    bufferedWriter.newLine();
                    bufferedWriter.write("Orbital parameters post-maneuver :");
                    bufferedWriter.newLine();
                    // 1) Dry mass
                    bufferedWriter.write(String.format("%.6f", dryMass));
                    bufferedWriter.newLine();
                    // 2) (m0 - finalState.mass)
                    bufferedWriter.write(String.format("%.6f", m0 - postBurnState.getMass()));
                    bufferedWriter.newLine();
                    // 3) (final mass - dry mass)
                    bufferedWriter.write(String.format("%.6f", postBurnState.getMass() - dryMass));
                    bufferedWriter.newLine();
                    // 4) finalOrbit date
                    bufferedWriter.write(postBurnState.getOrbit().getDate().toString());
                    bufferedWriter.newLine();
                    // 5) altitude
                    bufferedWriter.write(String.format("%.6f",
                            (postBurnState.getPVCoordinates().getPosition().getNorm() - 6378137.0) / 1000.0));
                    bufferedWriter.newLine();
                    // 6) SMA
                    bufferedWriter.write(String.format("%.12f", postBurnState.getA() / 1000.0));
                    bufferedWriter.newLine();
                    // 7) E
                    bufferedWriter.write(String.format("%.7f", postBurnState.getE()));
                    bufferedWriter.newLine();
                    // 8) i (deg)
                    bufferedWriter.write(String.format("%.4f", FastMath.toDegrees(postBurnState.getI())));
                    bufferedWriter.newLine();
                    // 9) RAAN
                    bufferedWriter.write(String.format("%.4f",
                            FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit) postBurnState.getOrbit()).getRightAscensionOfAscendingNode(), Math.PI))));
                    bufferedWriter.newLine();
                    // 10) Perigee argument
                    bufferedWriter.write(String.format("%.8f",
                            FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit) postBurnState.getOrbit()).getPerigeeArgument(), Math.PI))));
                    bufferedWriter.newLine();
                    // 11) Mean anomaly
                    bufferedWriter.write(String.format("%.8f",
                            FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit) postBurnState.getOrbit()).getMeanAnomaly(), Math.PI))));
                    bufferedWriter.newLine();
                    // 12) Mean motion (rev/day)
                    bufferedWriter.write(String.format("%.14f",
                            86400.0 * postBurnState.getOrbit().getKeplerianMeanMotion() / (2 * Math.PI)));
                    bufferedWriter.newLine();
                    // 13) Orbital period (s)
                    bufferedWriter.write(String.format("%.6f", postBurnState.getKeplerianPeriod()));
                    bufferedWriter.newLine();
                    bufferedWriter.close();

                    // Also record the post-maneuver date
                    writer = new FileWriter(lastManeuverDateFile, true);
                    bufferedWriter = new BufferedWriter(writer);
                    bufferedWriter.newLine();
                    bufferedWriter.write("Date post - maneuvre");
                    bufferedWriter.newLine();
                    bufferedWriter.write(String.valueOf(manoeuverEndDate));
                    bufferedWriter.newLine();
                    bufferedWriter.close();

                    double end = System.currentTimeMillis();
                    double duration = (end - start) / 1000.0;
                    System.out.println("Execution time:" + duration + " s");

                } catch (OrekitException var92) {
                    var92.printStackTrace();
                    System.exit(1);
                }
            } catch (IOException var93) {
                var93.printStackTrace();
            }
        } catch (IOException var94) {
            var94.printStackTrace();
        }
    }

    public static void calculateErgolConsumption() {
        double start = (double) System.currentTimeMillis();

        try {

            System.out.println("Using ergols input file: " + ergolsFile);
            System.out.println("Using ergols consumption output file: " + consommationErgolsFile);
            System.out.println("Mode parameter: " + modeParameter);
            // Load Orekit data
            File orekitData = new File("orekit-data");
            DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
            manager.addProvider(new DirectoryCrawler(orekitData));

            // Read input data from file
            List<String> lines = Files.readAllLines(Paths.get(ergolsFile));
            double DRYMASS = Double.parseDouble(lines.get(1));
            System.out.println("Dry Mass: " + DRYMASS);
            double ERGOL = Double.parseDouble(lines.get(2));
            System.out.println("Ergol: " + ERGOL);
            double ISP = Double.parseDouble(lines.get(5));
            System.out.println("ISP: " + ISP);
            String ManeuverType = lines.get(4);
            System.out.println("ManeuverType: " + ManeuverType);
            double DV = Double.parseDouble(lines.get(3));  // Delta-V unique
            System.out.println("Delta-V: " + DV);

            // Calcul de la masse initiale
            double initialMass = DRYMASS + ERGOL;

            // Calcul de la gravité à la surface de la Terre
            double g0 = 9.80665;

            // Calcul de la masse après la manœuvre
            double finalMassAfterManeuver = calculateFinalMass(initialMass, DV, ISP, g0);
            System.out.println(finalMassAfterManeuver);
            double ergolConsumed = initialMass - finalMassAfterManeuver;

            // Affichage de la consommation d'ergols
            System.out.println("Consommation d'ergols pour la manœuvre : " + ergolConsumed + " kg");

            // Écriture des résultats dans un fichier
            try (BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(consommationErgolsFile, true))) {
                bufferedWriter.write(String.valueOf(ergolConsumed));

            }

        } catch (Exception e) {
            // Gérer toutes les exceptions ici
            System.err.println("Une erreur s'est produite : " + e.getMessage());
            // Ajoutez vos instructions globales ici (log, affichage, etc.)
            e.printStackTrace();  // Affiche la pile d'appel pour déboguer (facultatif)
            try (BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(consommationErgolsFile, true))) {
                bufferedWriter.write(String.valueOf(0.0));
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }
    }
    /**
     * Returns true if dateToCheck >= referenceDate (within toleranceSeconds).
     */
    private static boolean isEqualOrAfterWithTolerance(AbsoluteDate dateToCheck, AbsoluteDate referenceDate, double toleranceSeconds) {
        double difference = dateToCheck.durationFrom(referenceDate);
        System.out.println("difference = " + difference);
        return difference >= -toleranceSeconds;
    }

    // Méthode pour calculer la masse finale après une manœuvre impulsionnelle
    public static double calculateFinalMass(double initialMass, double deltaV, double ISP, double g0) {
        // Tsiolkovsky rocket equation
        return initialMass * FastMath.exp(-(deltaV / (ISP * g0)));
    }


}
