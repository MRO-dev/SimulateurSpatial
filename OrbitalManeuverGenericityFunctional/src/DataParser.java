import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class DataParser {

    /**
     * Parses the given file into a Map where keys are header names and values are the corresponding values.
     * Each file is expected to have exactly one header line, followed by exactly N lines of data
     * (where N = number of header fields). Trailing blank lines are ignored.
     *
     * @param filePath The path to the data file.
     * @return A map of parameter names to their values.
     * @throws IOException if the file is empty, or if the # of data lines doesn't match # of columns
     */
    public static Map<String, String> parseDataFile(String filePath) throws IOException {

        // 1) Read all lines and trim them
        List<String> lines = Files.readAllLines(Paths.get(filePath))
                .stream()
                .map(String::trim)
                .collect(Collectors.toList());

        // 2) Remove trailing blank lines
        while (!lines.isEmpty() && lines.get(lines.size() - 1).isEmpty()) {
            lines.remove(lines.size() - 1);
        }

        // 3) Ensure there's at least one line for the header
        if (lines.isEmpty()) {
            throw new IOException("The file is empty or only contained blank lines: " + filePath);
        }

        // 4) The first line is the CSV header
        String headerLine = lines.get(0);
        String[] headers = headerLine.split(",");
        for (int i = 0; i < headers.length; i++) {
            headers[i] = headers[i].trim();
        }

        // 5) We expect exactly 1 line of data per header,
        //    so the total lines (minus the header line) must match headers.length
        int dataLines = lines.size() - 1;  // number of lines after the header
        if (dataLines != headers.length) {
            throw new IOException("Mismatch between header count and value count in file: " + filePath
                    + "\nHeaders: " + headers.length
                    + "  Data lines: " + dataLines);
        }

        // 6) Build the map from the remaining lines
        Map<String, String> dataMap = new HashMap<>();
        for (int i = 1; i < lines.size(); i++) {
            dataMap.put(headers[i - 1], lines.get(i));
        }

        return dataMap;
    }
}
