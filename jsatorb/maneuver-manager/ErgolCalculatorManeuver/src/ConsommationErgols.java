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

public class ConsommationErgols {

    private static final double MU = 3.986004418E14;  // Paramètre gravitationnel de la Terre (m^3/s^2)

    public static void main(String[] args) throws IOException, ParseException {
        try {
            double start = (double) System.currentTimeMillis();
            // Chargement des données Orekit
            File orekitData = new File("orekit-data");
            DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
            manager.addProvider(new DirectoryCrawler(orekitData));

            // Lecture des données d'entrée depuis le fichier
            double DRYMASS = Double.parseDouble(Files.readAllLines(Paths.get("Ergols.txt")).get(1));
            System.out.println("Masse à vide : " + DRYMASS + " kg");
            double ERGOL = Double.parseDouble(Files.readAllLines(Paths.get("Ergols.txt")).get(2));
            System.out.println("Masse d'ergols : " + ERGOL + " kg");
            double ISP = Double.parseDouble(Files.readAllLines(Paths.get("Ergols.txt")).get(3));
            System.out.println("Impulsion spécifique (ISP) : " + ISP + " s");
            double SMA = Double.parseDouble(Files.readAllLines(Paths.get("Ergols.txt")).get(5)) * 1000.0;
            System.out.println("SMA initial : " + SMA / 1000.0 + " km");
            double SMA_2 = Double.parseDouble(Files.readAllLines(Paths.get("Ergols.txt")).get(6)) * 1000.0; // SMA final en mètres
            System.out.println("SMA final : " + SMA_2 / 1000.0 + " km");

            // Calcul de la masse initiale
            double initialMass = DRYMASS + ERGOL;

            // Calculer les Delta-Vs pour le transfert de Hohmann
            double[] deltaVs = calculateDeltaVs(SMA, SMA_2);
            double DV1 = deltaVs[0];
            double DV2 = deltaVs[1];

            // Affichage des Delta-V
            System.out.println("Delta-V1 : " + DV1 + " m/s");
            System.out.println("Delta-V2 : " + DV2 + " m/s");

            // Accélération due à la gravité
            double g0 = 9.80665;

            // Calcul de la masse après la première manœuvre
            double finalMassAfterFirstManeuver = calculateFinalMass(initialMass, DV1, ISP, g0);
            double ergolConsumedFirstManeuver = initialMass - finalMassAfterFirstManeuver;

            // Calcul de la masse après la deuxième manœuvre
            double finalMassAfterSecondManeuver = calculateFinalMass(finalMassAfterFirstManeuver, DV2, ISP, g0);
            double ergolConsumedSecondManeuver = finalMassAfterFirstManeuver - finalMassAfterSecondManeuver;

            // Calcul total de la consommation d'ergols
            double totalErgolConsumed = ergolConsumedFirstManeuver + ergolConsumedSecondManeuver;

            // Affichage de la consommation d'ergols
            System.out.println("Consommation d'ergols pour la première manœuvre : " + ergolConsumedFirstManeuver + " kg");
            System.out.println("Consommation d'ergols pour la deuxième manœuvre : " + ergolConsumedSecondManeuver + " kg");
            System.out.println("Consommation totale d'ergols : " + totalErgolConsumed + " kg");

            // Enregistrement de la consommation totale dans un fichier
            FileWriter writer = new FileWriter("ConsommationErgols.txt", true);
            BufferedWriter bufferedWriter = new BufferedWriter(writer);
            bufferedWriter.write(String.valueOf(totalErgolConsumed));
            bufferedWriter.close();
        } catch (Exception e) {
            // Gestion des exceptions
            System.err.println("Une erreur s'est produite : " + e.getMessage());
            e.printStackTrace();
            try (BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter("ConsommationErgols.txt", true))) {
                bufferedWriter.write(String.valueOf(0.0));
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }
    }

    // Méthode pour calculer la masse finale après une manœuvre
    public static double calculateFinalMass(double initialMass, double deltaV, double ISP, double g0) {
        return initialMass * FastMath.exp(-Math.abs(deltaV) / (ISP * g0));
    }

    // Méthode pour calculer les Delta-Vs pour le transfert de Hohmann
    public static double[] calculateDeltaVs(double SMA1, double SMA2) {
        double mu = MU;
        double r1 = SMA1;
        double r2 = SMA2;

        // Vitesses orbitales circulaires initiale et finale
        double v_circ1 = Math.sqrt(mu / r1);
        double v_circ2 = Math.sqrt(mu / r2);

        // Semi-grand axe de l'orbite de transfert
        double a_transfer = (r1 + r2) / 2.0;

        // Vitesses sur l'orbite de transfert au périgée et à l'apogée
        double v_transf_perigee = Math.sqrt(2 * mu * r2 / (r1 * (r1 + r2)));
        double v_transf_apogee = Math.sqrt(2 * mu * r1 / (r2 * (r1 + r2)));

        // Calcul des Delta-V en magnitudes absolues
        double deltaV1 = Math.abs(v_transf_perigee - v_circ1);
        double deltaV2 = Math.abs(v_circ2 - v_transf_apogee);

        return new double[]{deltaV1, deltaV2};
    }
}
