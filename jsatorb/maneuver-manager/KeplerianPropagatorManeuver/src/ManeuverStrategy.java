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
    void computeAndExecute() throws IOException, ParseException, InterruptedException;
    void calculateErgolConsumption() throws IOException;
    void processReachOrbitTime() throws IOException;

    void loadMassData(Boolean isMassCalculation) throws IOException;
    void loadTimeData(Boolean isTimeCalculation) throws IOException;
}
