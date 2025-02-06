import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class DataParser {

    /**
     * Parses the given file into a Map where keys are header names and values are the corresponding values.
     *
     * @param filePath The path to the data file.
     * @return A map of parameter names to their values.
     * @throws IOException if the file cannot be read.
     */
    public static Map<String, String> parseDataFile(String filePath) throws IOException {
        List<String> lines = Files.readAllLines(Paths.get(filePath))
                .stream()
                .filter(line -> !line.trim().isEmpty())
                .collect(Collectors.toList());

        if (lines.isEmpty()) {
            throw new IOException("The file is empty: " + filePath);
        }

        // Split header and trim whitespace
        String headerLine = lines.get(0);
        String[] headers = headerLine.split(",");
        for (int i = 0; i < headers.length; i++) {
            headers[i] = headers[i].trim();
        }

        // Expecting one value per header line
        if (lines.size() - 1 != headers.length) {
            throw new IOException("Mismatch between header count and value count in file: " + filePath);
        }

        Map<String, String> dataMap = new HashMap<>();
        for (int i = 1; i < lines.size(); i++) {
            dataMap.put(headers[i - 1], lines.get(i).trim());
        }
        return dataMap;
    }
}
