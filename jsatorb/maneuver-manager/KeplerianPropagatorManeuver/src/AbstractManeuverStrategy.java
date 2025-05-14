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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Abstract parent class for maneuver strategies.
 * It factors out the common fields and initialization steps:
 *  - Handling modeParameter overrides
 *  - Defining typical file names
 *  - Possibly reading certain system properties or data
 *  - Setting up Orekit data, etc.
 */
public abstract class AbstractManeuverStrategy implements ManeuverStrategy {

    // -------------------------------------------------
    // Common files and parameters
    // -------------------------------------------------
    protected String dataFile = "Data.txt";
    protected String maneuvFile = "Maneuv.txt";
    protected String lastManeuverDateFile = "LastManeuverDate.txt";
    protected String timePersistenceFile = "time-persistence.json";
    protected String publishTopic = "resultat/fichier";
    protected String triggerTopic = "trigger/continuation";
    protected String resultFileName = "Result.txt";
    protected String postManeuverDateFileName = "PostManeuverDate.txt";
    protected String ergolsFile = "Ergols.txt";
    protected String consommationErgolsFile = "ConsommationErgols.txt";
    protected String maneuverOrderFile = "ManeuverOrder.txt";
    protected String commandDataFile = "CommandData.txt";
    protected String timeIntermediateParametersFile = "TimeIntermediateParameters.txt";
    protected String apsideFile = "ApsideDates.txt";

    // Typically used constants or fields
    protected static final double MU = 3.986004418E14;
    protected static double g0 = 9.80665;
    protected static final double TIME_TOLERANCE_SECONDS = 1e-3;
    protected File orekitData = new File("orekit-data");
    protected DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();

    // Additional parameters
    protected String modeParameter;
    protected String commandParameter;

    // Example orbit parameters to be loaded
    protected double SMA;
    protected double ECC;
    protected double INC;
    protected double RAAN;
    protected double PA;
    protected double ANO;
    protected double DRYMASS;
    protected double ISP;
    protected double ERGOL;
    protected String DATE;
    protected String APSIDE_DATE;
    protected String MANEUV_TYPE;
    protected String endDateString;
    protected double m0 = DRYMASS + ERGOL;
    protected  double manoeuverRelativeDate;
    protected static int item = 1;

    protected Map<String, String> ergolsMap = DataParser.parseDataFile(ergolsFile);
    protected Map<String, String> tipData = DataParser.parseDataFile(timeIntermediateParametersFile);
    protected Map<String, String> cmdData = DataParser.parseDataFile(commandDataFile);


    // -------------------------------------------------
    // Constructor
    // -------------------------------------------------
    public AbstractManeuverStrategy(String modeParameter, String commandParameter) throws IOException {
        this.modeParameter = modeParameter;
        this.commandParameter = commandParameter;

        // Optionally read some system properties
        // to override defaults if you like:
        this.dataFile = System.getProperty("dataFile", this.dataFile);
        this.maneuvFile = System.getProperty("maneuvFile", this.maneuvFile);
        this.lastManeuverDateFile = System.getProperty("lastManeuverDateFile", this.lastManeuverDateFile);
        this.timePersistenceFile = System.getProperty("timePersistenceFile", this.timePersistenceFile);

        // 2) Now override filenames based on modeParameter ("blue2", "red1", etc.)
        overrideFilenamesForMode(this.modeParameter);

        // Read the LastManeuverDate file
        List<String> allLines = null;
        try {
            allLines = Files.readAllLines(Paths.get(lastManeuverDateFile));
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
            // Read the maneuv file for the maneuver relative date
            List<String> maneuvLines = Files.readAllLines(Paths.get(maneuvFile));
            manoeuverRelativeDate = Double.parseDouble(maneuvLines.get(item));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        // Now handle the different "color+number" modes:
        overrideFilenamesForMode(this.modeParameter);

        // Load the Orekit data manager:
        manager.addProvider(new DirectoryCrawler(orekitData));
    }
    /**
     * Child classes can call this if they also need to read from
     * dataFile or other text files, but this method centralizes
     * reading the usual "Data.txt" (or "DataX.txt").
     */
    public void loadMassData(Boolean isMassCalculation) throws IOException{
        manager.addProvider(new DirectoryCrawler(orekitData));
        System.out.println("Using ergols input file: " + ergolsFile);
        System.out.println("Using ergols consumption output file: " + consommationErgolsFile);
        System.out.println("Mode parameter: " + modeParameter);
        // Chargement des données Orekit
        File orekitData = new File("orekit-data");
        DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
        manager.addProvider(new DirectoryCrawler(orekitData));
        // Parse the ergols file
        ergolsMap = DataParser.parseDataFile(ergolsFile);
        cmdData = DataParser.parseDataFile(commandDataFile);

        // Assign fields
        this.DRYMASS = Double.parseDouble(ergolsMap.get("DryMass"));
        this.ERGOL   = Double.parseDouble(ergolsMap.get("ErgolMass"));
        this.ISP     = Double.parseDouble(ergolsMap.get("ISP"));
        System.out.println("isMassCalculation: " + isMassCalculation);
        System.out.println(commandDataFile);
        System.out.println(cmdData);
        // If you want the maneuver type from this file, you can do:
        // MANEUV_TYPE = ergolsMap.get("ManeuverType");
       if (isMassCalculation) {
           this.SMA     = Double.parseDouble(cmdData.get("SMA")) * 1000.0;
       }
        // Convert SMA_1, SMA_2 from km to meters

        //this.SMA_2   = Double.parseDouble(ergolsMap.get("SMA_2")) * 1000.0;

        // Debug prints
        System.out.println("Using ergols input file: " + ergolsFile);
        System.out.println("Masse à vide : " + DRYMASS + " kg");
        System.out.println("Masse d'ergols : " + ERGOL + " kg");
        System.out.println("Impulsion spécifique (ISP) : " + ISP + " s");
        System.out.println("SMA initial : " + (SMA / 1000.0) + " km");
        //System.out.println("SMA final : " + (SMA_2 / 1000.0) + " km");
    };

    public void loadTimeData(Boolean isTimeCalculation) throws IOException{
        // 1) Load Orekit data provider
        manager.addProvider(new DirectoryCrawler(orekitData));
        System.out.println("TEST");
        if (isTimeCalculation) {
            this.DATE = Files.readAllLines(Paths.get(maneuverOrderFile)).get(1);
            System.out.println("Date from maneuverOrderFile: " + DATE);

        }
        cmdData = DataParser.parseDataFile(commandDataFile);
        // the second line is your date/timestamp
        // 3) Parse the commandDataFile via DataParser (CSV approach)
        //    If the first column is “DATE” but it’s blank (and you don’t want to use it),
        //    just skip it and parse the others as needed:

        // Now read the fields we care about (SMA, ECC, INC, RAAN, AoP, MeanAnom, etc.)
        // Each key must match the header line in commandDataFile
        // e.g. "SMA", "ECC", "INC", "RAAN", "AoP", "MeanAnom", "Dry Mass", ...
        this.SMA  = Double.parseDouble(cmdData.get("SMA"))  * 1000.0;
        this.ECC  = Double.parseDouble(cmdData.get("ECC"));
        // For INC, RAAN, AoP, MeanAnom, we do Math.toRadians(...) per your code
        this.INC  = FastMath.toRadians(Double.parseDouble(cmdData.get("INC")));
        this.RAAN = FastMath.toRadians(Double.parseDouble(cmdData.get("RAAN")));
        this.PA   = FastMath.toRadians(Double.parseDouble(cmdData.get("AoP")));
        this.ANO  = FastMath.toRadians(Double.parseDouble(cmdData.get("MeanAnom")));

        // 4) Parse the ergolsFile
        Map<String, String> ergolsData = DataParser.parseDataFile(ergolsFile);
        this.DRYMASS = Double.parseDouble(ergolsData.get("DryMass"));
        this.ERGOL   = Double.parseDouble(ergolsData.get("ErgolMass"));
        // You said you hard-coded ISP to 1.0 in the old code, but if you do want the real
        // value from ergolsFile:
        //  this.ISP = Double.parseDouble(ergolsData.get("ISP"));
        // or if you want to keep it forced to 1.0:
        this.ISP = 1.0;

        // 5) Parse the timeIntermediateParametersFile
        tipData = DataParser.parseDataFile(timeIntermediateParametersFile);
        this.MANEUV_TYPE = tipData.get("ManeuverType");  // e.g. "Hohmann"
        // Convert to meters


        int item = 1;
        this.m0 = DRYMASS + ERGOL;
        this.manoeuverRelativeDate = Double.parseDouble(Files.readAllLines(Paths.get(maneuvFile)).get(item));
        this.g0 = 9.80665;

        // Print or log to confirm:
        System.out.println("Loaded time/orbit data from:");
        System.out.println(" - maneuverOrderFile => DATE: " + DATE);
        System.out.println(" - commandDataFile => SMA: " + SMA + " m, ECC: " + ECC + ", etc.");
        System.out.println(" - ergolsFile => DryMass: " + DRYMASS + " kg, Ergol: " + ERGOL);
        System.out.println(" - maneuvFile => manoeuverRelativeDate: " + manoeuverRelativeDate);
    };


    /**
     * Example method that sets filenames, topics, etc.,
     * depending on the given modeParameter.
     */
    protected void overrideFilenamesForMode(String modeParameter) {
        if ("blue2".equalsIgnoreCase(modeParameter)) {
            this.dataFile                 = "Data2.txt";
            this.maneuvFile               = "Maneuv2.txt";
            this.lastManeuverDateFile     = "LastManeuverDate2.txt";
            this.publishTopic             = "resultat/fichier2";
            this.triggerTopic             = "trigger/continuation2";
            this.resultFileName           = "Result2.txt";
            this.postManeuverDateFileName = "PostManeuverDate2.txt";
            this.ergolsFile               = "Ergols2.txt";
            this.consommationErgolsFile   = "ConsommationErgols2.txt";
            this.maneuverOrderFile        = "ManeuverOrder2.txt";
            this.commandDataFile          = "CommandData2.txt";
            this.timeIntermediateParametersFile = "TimeIntermediateParameters2.txt";
            this.apsideFile               = "ApsideDates2.txt";
        }
        else if ("red1".equalsIgnoreCase(modeParameter)) {
            this.dataFile                 = "Data3.txt";
            this.maneuvFile               = "Maneuv3.txt";
            this.lastManeuverDateFile     = "LastManeuverDate3.txt";
            this.publishTopic             = "resultat/fichier3";
            this.triggerTopic             = "trigger/continuation3";
            this.resultFileName           = "Result3.txt";
            this.postManeuverDateFileName = "PostManeuverDate3.txt";
            this.ergolsFile               = "Ergols3.txt";
            this.consommationErgolsFile   = "ConsommationErgols3.txt";
            this.maneuverOrderFile        = "ManeuverOrder3.txt";
            this.commandDataFile          = "CommandData3.txt";
            this.timeIntermediateParametersFile = "TimeIntermediateParameters3.txt";
            this.apsideFile               = "ApsideDates3.txt";
        }
        else if ("red2".equalsIgnoreCase(modeParameter)) {
            this.dataFile                 = "Data4.txt";
            this.maneuvFile               = "Maneuv4.txt";
            this.lastManeuverDateFile     = "LastManeuverDate4.txt";
            this.publishTopic             = "resultat/fichier4";
            this.triggerTopic             = "trigger/continuation4";
            this.resultFileName           = "Result4.txt";
            this.postManeuverDateFileName = "PostManeuverDate4.txt";
            this.ergolsFile               = "Ergols4.txt";
            this.consommationErgolsFile   = "ConsommationErgols4.txt";
            this.maneuverOrderFile        = "ManeuverOrder4.txt";
            this.commandDataFile          = "CommandData4.txt";
            this.timeIntermediateParametersFile = "TimeIntermediateParameters4.txt";
            this.apsideFile               = "ApsideDates4.txt";
        }
        // else remain the defaults for "blue1" or unknown modes
    }

    // Example method to compare two AbsoluteDate objects within a tolerance
    protected static boolean isEqualOrAfterWithTolerance(AbsoluteDate dateToCheck, AbsoluteDate referenceDate, double toleranceSeconds) {
        // Calculate the difference: positive if dateToCheck is after referenceDate
        double difference = dateToCheck.durationFrom(referenceDate);

        // If difference is greater than or within the negative tolerance, consider it as valid
        return difference >= -toleranceSeconds;
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
//        System.out.println(String.format("Target SMA: %.2f km", SMA_2/1000.0));
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

    // -------------------------------------------------
    // The methods from ManeuverStrategy interface:
    //    computeAndExecute(), calculateErgolConsumption(),
    //    processReachOrbitTime().
    // Each child can override them differently.
    // -------------------------------------------------


    /**
     * Writes a timestamp to the specified file in the required format
     */
    public void writeManeuverTimestamp(String fileName, AbsoluteDate maneuverDate) throws IOException {
        // Calculate timestamp
        AbsoluteDate epoch = new AbsoluteDate(1970, 1, 1, 0, 0, 0.000, TimeScalesFactory.getUTC());
        double durationInSeconds = maneuverDate.durationFrom(epoch);
        long timestampMillis = Math.round(durationInSeconds * 1000);

        // Write to file
        try (FileWriter writer = new FileWriter(fileName, true);
             BufferedWriter bufferedWriter = new BufferedWriter(writer)) {
            bufferedWriter.newLine();
            bufferedWriter.write("Date post - maneuver");
            bufferedWriter.newLine();

            // For lastManeuverDateFile, write the AbsoluteDate string
            if (fileName.equals(lastManeuverDateFile)) {
                bufferedWriter.write(String.valueOf(maneuverDate));
            } else {
                // For other files, write the timestamp in milliseconds
                bufferedWriter.write(String.valueOf(timestampMillis));
            }

            bufferedWriter.newLine();
        }
    }

    @Override
    public abstract void computeAndExecute() throws IOException, java.text.ParseException;

    @Override
    public abstract void calculateErgolConsumption() throws IOException;

    @Override
    public abstract void processReachOrbitTime() throws IOException;

}
