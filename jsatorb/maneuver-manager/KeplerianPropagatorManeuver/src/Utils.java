import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.jetbrains.annotations.Nullable;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.propagation.SpacecraftState;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.time.Instant;
import java.time.ZoneId;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;
import java.util.HashMap;
import java.util.Map;
/*
 * This class contains utility methods for logging results and creating JSON-like payloads.
 * It is used in the HohmannTransfert and HohmannStrategy classes.
 */
public class Utils {

    // Method to log results in Result.txt
    static void logResults(String filePath, SpacecraftState finalState, Orbit finalOrbit, double m0, double DRYMASS) throws IOException {
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
        // convert to degrees
        double degI      = FastMath.toDegrees(finalState.getI());
// define your minimum swing (in degrees)
        final double MIN_DELTA = 0.001;

// if it’s too small, bump to ±MIN_DELTA preserving the original sign
        double outI = (FastMath.abs(degI) < MIN_DELTA)
                ? FastMath.copySign(MIN_DELTA, degI)
                : degI;

// write it out
        bufferedWriter.newLine();
        bufferedWriter.write(String.format("%.4f", outI));

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

    /**
     * Logs the apside date to a specified file.
     *
     * @param filePath Path to the log file.
     * @param date     The date of the apside event.
     * @throws IOException If an I/O error occurs.
     */
    static void logApsideDate(String filePath, @Nullable AbsoluteDate date) throws IOException {
        try (BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(filePath, true))) {
            if (date == null) {
                bufferedWriter.newLine();
                bufferedWriter.write("");
                bufferedWriter.newLine();
                bufferedWriter.write("");
                bufferedWriter.newLine();
            }else{
                bufferedWriter.newLine();
                bufferedWriter.write("Apside reached at:");
                bufferedWriter.newLine();
                bufferedWriter.write(date.toString());
                bufferedWriter.newLine();
            }
        }
    }

    static void appendLineToFile(String filePath, String content) throws IOException {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(filePath, true))) {
            bw.write(content);
            bw.close();
        }
    }

    /**
     * Create a JSON-like payload with start and end dates (ISO format).
     */
    static Map<String, String> createDatePayload(AbsoluteDate startDate, AbsoluteDate endDate) {
        DateTimeFormatter formatter = DateTimeFormatter.ISO_INSTANT;
        String startIso = formatter.format(startDate.toDate(TimeScalesFactory.getUTC()).toInstant());
        String endIso = formatter.format(endDate.toDate(TimeScalesFactory.getUTC()).toInstant());

        Map<String, String> payload = new HashMap<>();
        payload.put("startdate", startIso);
        payload.put("enddate", endIso);
        return payload;
    }

    static void writeJsonPayload(Map<String, String> payload) throws IOException {
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

    static AbsoluteDate parseDateFromTimestamp(String timeStampStr) {
        long ts = Long.parseLong(timeStampStr);
        Instant instant = Instant.ofEpochMilli(ts);
        ZonedDateTime zdt = instant.atZone(ZoneId.of("UTC"));

        int year = zdt.getYear();
        int month = zdt.getMonthValue();
        int day = zdt.getDayOfMonth();
        int hour = zdt.getHour();
        int minute = zdt.getMinute();
        double second = zdt.getSecond() + zdt.getNano() / 1e9;

        return new AbsoluteDate(year, month, day, hour, minute, second, TimeScalesFactory.getUTC());
    }
}
