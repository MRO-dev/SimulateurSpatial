import org.hipparchus.util.FastMath;
import org.orekit.data.DataContext;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.ParseException;
import java.util.List;

public class ConsommationErgolsManual {

    private static final double MU = 3.986004418E14;
    // Default file names for ergols input and consumption output
    private static String ergolsFile = "Ergols.txt";
    private static String consommationErgolsFile = "ConsommationErgols.txt";

    // Mode parameter (default is "blue1")
    private static String modeParameter = "blue1";

    // In a static block, override the default file names and mode parameter using system properties
    static {
        ergolsFile = System.getProperty("ergolsFile", ergolsFile);
        consommationErgolsFile = System.getProperty("consommationErgolsFile", consommationErgolsFile);
        modeParameter = System.getProperty("modeParameter", modeParameter);

        // Optionally, if you want different files for a given mode (for example, "blue2")
        if ("blue2".equalsIgnoreCase(modeParameter)) {
            // Override the file names for blue2 mode:
            ergolsFile = "Ergols2.txt";
            consommationErgolsFile = "ConsommationErgols2.txt";
        }
        // Optionally, if you want different files for a given mode (for example, "blue2")
        if ("red1".equalsIgnoreCase(modeParameter)) {
            // Override the file names for blue2 mode:
            ergolsFile = "Ergols3.txt";
            consommationErgolsFile = "ConsommationErgols3.txt";
        }
        // Optionally, if you want different files for a given mode (for example, "blue2")
        if ("red2".equalsIgnoreCase(modeParameter)) {
            // Override the file names for blue2 mode:
            ergolsFile = "Ergols4.txt";
            consommationErgolsFile = "ConsommationErgols4.txt";
        }
    }// Earth's gravitational parameter (m^3/s^2)

    public static void main(String[] args) {
        double start = (double) System.currentTimeMillis();

        try {

            if (args.length >= 2) {
                ergolsFile = args[0];
                consommationErgolsFile = args[1];
            }

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

    // Méthode pour calculer la masse finale après une manœuvre impulsionnelle
    public static double calculateFinalMass(double initialMass, double deltaV, double ISP, double g0) {
        // Tsiolkovsky rocket equation
        return initialMass * FastMath.exp(-(deltaV / (ISP * g0)));
    }

}
