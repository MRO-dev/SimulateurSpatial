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
        // 1) Parse arguments or system properties
        // ------------------------------------------------------------------
        // For brevity, we'll assume you have a small helper that
        // reads your data file, maneuvFile, etc.
        // (You can reuse your existing DataParser, etc.)
        Map<String, String> userData = DataParser.parseDataFile("Data.txt");  // Just an example

        // or read them from system properties or from the args:
        String modeParameter     = System.getProperty("modeParameter", "blue1");
        String commandParameter  = System.getProperty("commandParameter", "compute");
        String maneuverType      = userData.get("ManeuverType");  // "Hohmann" or "QLaw", etc.

        // Suppose we also read some relevant orbit parameters from userData:
        // (Those were previously in your static fields.)
        double sma       = Double.parseDouble(userData.get("SMA")) * 1000.0;
        double ecc       = Double.parseDouble(userData.get("ECC"));
        double incDeg    = Double.parseDouble(userData.get("INC"));    // degrees
        double raanDeg   = Double.parseDouble(userData.get("RAAN"));   // degrees
        double aopDeg    = Double.parseDouble(userData.get("AoP"));    // argument of perigee, deg
        double meanAnom  = Double.parseDouble(userData.get("MeanAnom"));
        double dryMass   = Double.parseDouble(userData.get("Dry Mass"));
        double isp       = Double.parseDouble(userData.get("ISP"));
        double ergol     = Double.parseDouble(userData.get("Ergol mass"));
        String dateStr   = userData.get("DATE");
        // etc.

        // ------------------------------------------------------------------
        // 2) Instantiate the appropriate ManeuverStrategy
        // ------------------------------------------------------------------
        ManeuverStrategy strategy = null;

        if ("Hohmann".equalsIgnoreCase(maneuverType)) {
            double sma2      = Double.parseDouble(userData.get("SMA_2")) * 1000.0; // for Hohmann target
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
        }
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
            strategy.computeAndExecute();
        } else if (strategy != null && "calculateMass".equalsIgnoreCase(commandParameter)) {
            strategy.calculateErgolConsumption();
        } else if (strategy != null && "determineOrbitInstant".equalsIgnoreCase(commandParameter)) {
            strategy.processReachOrbitTime();
        }else {
            System.out.printf("No valid strategy or commandParameter '%s' => nothing to do.%n", commandParameter);
        }
    }
}
