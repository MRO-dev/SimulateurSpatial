import java.io.IOException;
import java.text.ParseException;

/**
 * Strategy interface for different orbital maneuvers.
 */
public interface ManeuverStrategy {
    /**
     * Perform the necessary computations and produce the final results,
     * including any file outputs or MQTT steps.
     */
    void computeAndExecute() throws IOException, ParseException;
    void calculateErgolConsumption() throws IOException;
    void processReachOrbitTime() throws IOException;

    void loadMassData() throws IOException;
    void loadTimeData() throws IOException;
}
