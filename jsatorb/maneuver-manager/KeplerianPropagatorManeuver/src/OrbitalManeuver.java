import java.io.IOException;
import java.text.ParseException;
import java.util.Map;

/**
 * Main class responsible for:
 *  - Parsing arguments
 *  - Building the chosen ManeuverStrategy
 *  - Executing the strategy
 */
public class OrbitalManeuver {

    public static void main(String[] args) throws IOException, ParseException {

        // ------------------------------------------------------------------
        // 1) Read modeParameter, commandParameter from system properties
        // ------------------------------------------------------------------
        String modeParameter    = System.getProperty("modeParameter", "blue1");
        String commandParameter = System.getProperty("commandParameter", "compute");
        String dataFileToParse = getString(modeParameter);
        String commandDataFileToParse = getCommandDataString(modeParameter);

// Then parse:
        Map<String, String> userData = DataParser.parseDataFile(dataFileToParse);
        Map<String, String> commandData = DataParser.parseDataFile(commandDataFileToParse);

        String maneuverType      = commandData.get("ManeuverType");  // "Hohmann" or "QLaw", etc.

        // Suppose we also read some relevant orbit parameters from userData:
        // (Those were previously in your static fields.)
        double sma       = Double.parseDouble(commandData.get("SMA")) * 1000.0;
        double ecc       = Double.parseDouble(commandData.get("ECC"));
        double incDeg    = Double.parseDouble(commandData.get("INC"));    // degrees
        double raanDeg   = Double.parseDouble(commandData.get("RAAN"));   // degrees
        double aopDeg    = Double.parseDouble(commandData.get("AoP"));    // argument of perigee, deg
        double meanAnom  = Double.parseDouble(commandData.get("MeanAnom"));
        double dryMass   = Double.parseDouble(commandData.get("Dry Mass"));
        double isp       = Double.parseDouble(commandData.get("ISP"));
        double ergol     = Double.parseDouble(commandData.get("Ergol mass"));
        String dateStr   = userData.get("DATE");
        // etc.

        // ------------------------------------------------------------------
        // 2) Instantiate the appropriate ManeuverStrategy
        // ------------------------------------------------------------------
        ManeuverStrategy strategy = null;

        if ("Hohmann".equalsIgnoreCase(maneuverType)) {
            double sma2      = Double.parseDouble(commandData.get("SMA_2")) * 1000.0;
            //System.out.println("Inc "+incDeg);// for Hohmann target
            strategy = new HohmannStrategy(
                    sma,
                    ecc,
                    incDeg,
                    raanDeg,
                    aopDeg,
                    meanAnom,
                    dryMass,
                    ergol,
                    isp,
                    sma2,
                    dateStr,
                    modeParameter,
                    commandParameter
            );
        }else if ("Inclinaison".equalsIgnoreCase(maneuverType)) {
            // parse raw value
            double rawIncli2 = Double.parseDouble(commandData.get("INC_2"));

// define your smallest acceptable swing (in degrees, here)
            final double MIN_DELTA = 0.001;

// if the requested Δinclination is too small (including 0.0 or -0.0),
// bump it up to ±MIN_DELTA, preserving the original sign
            double incli2 = (Math.abs(rawIncli2) < MIN_DELTA)
                    ? Math.copySign(MIN_DELTA, rawIncli2)
                    : rawIncli2;

            System.out.println("Using target inclination (clamped): " + incli2);



            System.out.println("Inc " + incDeg);// for Incli target
            strategy = new InclinaisonStrategy(
                    sma,
                    ecc,
                    incDeg,
                    raanDeg,
                    aopDeg,
                    meanAnom,
                    dryMass,
                    ergol,
                    isp,
                    incli2,
                    dateStr,
                    modeParameter,
                    commandParameter
            );
        }else if ("Phasage".equalsIgnoreCase(maneuverType)) {
            double delta_theta = Double.parseDouble(commandData.get("DELTA_THETA"));
            System.out.println("Inc "+incDeg);// for Incli target
            strategy = new PhasageStrategy(
                    sma,
                    ecc,
                    incDeg,
                    raanDeg,
                    aopDeg,
                    meanAnom,
                    dryMass,
                    ergol,
                    isp,
                    delta_theta,
                    dateStr,
                    modeParameter,
                    commandParameter
            );
        }

//        }else if ("Bi-elliptic".equalsIgnoreCase(maneuverType)) {
//            double incli2 = Double.parseDouble(commandData.get("Inclination_2"));
//            System.out.println("Inc "+incDeg);// for Incli target
//            strategy = new InclinaisonStrategy(
//                    sma,
//                    ecc,
//                    incDeg,
//                    raanDeg,
//                    aopDeg,
//                    meanAnom,
//                    dryMass,
//                    ergol,
//                    isp,
//                    sma2,
//                    dateStr,
//                    modeParameter,
//                    commandParameter
//            );}
//
        // else if ("QLaw".equalsIgnoreCase(maneuverType)) {
        //    strategy = new QLawStrategy(...);
        // }
        else {
            // fallback or error
            System.err.println("Unsupported maneuverType: " + maneuverType);
            System.exit(1);
        }

        // ------------------------------------------------------------------
        // 3) Execute the strategy
        //    The strategy can handle:
        //      - reading/writing the relevant output files
        //      - calling the MQTT service
        //      - doing the actual orbit calculations
        // ------------------------------------------------------------------
        if (strategy != null && "compute".equalsIgnoreCase(commandParameter)) {
            strategy.loadTimeData(Boolean.FALSE);
            strategy.loadMassData(Boolean.FALSE);
            strategy.computeAndExecute();
        } else if (strategy != null && "calculateMass".equalsIgnoreCase(commandParameter)) {
            strategy.loadMassData(Boolean.TRUE);
            strategy.calculateErgolConsumption();
        } else if (strategy != null && "determineOrbitInstant".equalsIgnoreCase(commandParameter)) {
            System.out.println("Load Time Data");
            strategy.loadTimeData(Boolean.TRUE);
            strategy.processReachOrbitTime();
        }else {
            System.out.printf("No valid strategy or commandParameter '%s' => nothing to do.%n", commandParameter);
        }
    }

    private static String getString(String modeParameter) {
        String dataFileToParse = "Data.txt"; // default is "blue1"

// If modeParameter is "blue2" => Data2.txt
        if ("blue2".equalsIgnoreCase(modeParameter)) {
            dataFileToParse = "Data2.txt";
        }
// If modeParameter is "red1" => Data3.txt
        else if ("red1".equalsIgnoreCase(modeParameter)) {
            dataFileToParse = "Data3.txt";
        }
// If modeParameter is "red2" => Data4.txt
        else if ("red2".equalsIgnoreCase(modeParameter)) {
            dataFileToParse = "Data4.txt";
        }
        return dataFileToParse;
    }

    private static String getCommandDataString(String modeParameter) {
        String dataFileToParse = "CommandData.txt"; // default is "blue1"

// If modeParameter is "blue2" => Data2.txt
        if ("blue2".equalsIgnoreCase(modeParameter)) {
            dataFileToParse = "CommandData2.txt";
        }
// If modeParameter is "red1" => Data3.txt
        else if ("red1".equalsIgnoreCase(modeParameter)) {
            dataFileToParse = "CommandData3.txt";
        }
// If modeParameter is "red2" => Data4.txt
        else if ("red2".equalsIgnoreCase(modeParameter)) {
            dataFileToParse = "CommandData4.txt";
        }
        return dataFileToParse;
    }

}
