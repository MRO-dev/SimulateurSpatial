//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by FernFlower decompiler)
//

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.ParseException;
import java.util.List;
import java.util.Locale;
import org.hipparchus.util.FastMath;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.MathUtils;
import org.json.JSONObject;
import org.orekit.data.DataContext;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.forces.maneuvers.ImpulseManeuver;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.orekit.propagation.events.DateDetector;
import org.orekit.propagation.events.EventDetector;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;

public class KeplerianManeuver {
    private static final double TIME_TOLERANCE_SECONDS = 1e-3;
    static String APSIDE_DATE;
    static String endDateString;
    // 1 millisecond tolerance
    public KeplerianManeuver() {
    }

    public static void main(String[] args) throws NullPointerException, ClassCastException, IOException, ParseException {
        double start = (double)System.currentTimeMillis();
        List<String> allLines = Files.readAllLines(Paths.get("LastManeuverDate.txt"));

        // Initialize APSIDE_DATE to null or empty
        String extractedApsideDate = null;

        // Iterate over the lines in reverse order
        for (int i = allLines.size() - 1; i >= 0; i--) {
            String line = allLines.get(i).trim();
            // Skip empty lines and lines starting with the prefix
            if (!line.isEmpty() && !line.startsWith("Apside reached at:")) {
                extractedApsideDate = line;
                break;
            }
        }

        // Validate that a date was found
        if (extractedApsideDate == null) {
            throw new IOException("No valid apside date found in file.");
        }

        // Assign the extracted date to APSIDE_DATE
        APSIDE_DATE = extractedApsideDate;
        String fileContent = new String(Files.readAllBytes(Paths.get("time-persistence.json")));

        // Step 2: Parse it using org.json
        JSONObject jsonObject = new JSONObject(fileContent);

        // Step 3: Retrieve the "enddate" field as a string
        endDateString = jsonObject.optString("enddate");
        if (endDateString == null || endDateString.isEmpty()) {
            throw new IllegalArgumentException("No valid 'enddate' found in time-persistence.json");
        }
        System.out.println("End date string: " + endDateString);

        try {
            int num = 1;
            String DATE = (String)Files.readAllLines(Paths.get("Data.txt")).get(num);
            num = 2;
            double SMA = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 3;
            double ECC = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 4;
            double INC = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 5;
            double RAAN = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 6;
            double PA = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 7;
            double ANO = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            double sma = SMA * 1000.0;
            double ecc = ECC;
            double inc = FastMath.toRadians(INC);
            double raan = FastMath.toRadians(RAAN);
            double aop = FastMath.toRadians(PA);
            double ano = FastMath.toRadians(ANO);
            double MU = 3.986004418E14;
            Frame GCRF = FramesFactory.getEME2000();
            num = 8;
            double DRYMASS = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            double dryMass = DRYMASS;
            num = 10;
            double THRUST = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 11;
            double ISP = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 12;
            double ERGOL = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 13;
            double DURA = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 14;
            String maneuverType = (String)Files.readAllLines(Paths.get("Data.txt")).get(num);
            num = 15;
            double DV = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            System.out.println("maneuverType" + maneuverType);
            double m0 = dryMass + ERGOL;
            //System.out.println("integratorType" + integratorType);

            try {
                int item = 1;
                double manoeuverRelativeDate = Double.parseDouble((String)Files.readAllLines(Paths.get("Maneuv.txt")).get(item));
                item = 2;
                double DVx = Double.parseDouble((String)Files.readAllLines(Paths.get("Maneuv.txt")).get(item));
                item = 3;
                double DVy = Double.parseDouble((String)Files.readAllLines(Paths.get("Maneuv.txt")).get(item));
                item = 4;
                double DVz = Double.parseDouble((String)Files.readAllLines(Paths.get("Maneuv.txt")).get(item));
                item = 5;
                double durationOfManoeuver = Double.parseDouble((String)Files.readAllLines(Paths.get("Maneuv.txt")).get(item));
                Vector3D direction = new Vector3D(DVx, DVy, DVz);
                System.out.println("Vecteur direction X= " + direction.getX() + " Y= " + direction.getY() + " Z= " + direction.getZ());

                try {
                    File home = new File("/app/maneuver-manager/");
                    File orekitData = new File(Paths.get("orekit-data").toString());
                    if (!orekitData.exists()) {
                        System.err.format(Locale.US, "Failed to find %s folder%n", orekitData.getAbsolutePath());
                        System.err.format(Locale.US, "You need to download %s from %s, unzip it in %s and rename it 'orekit-data' for this tutorial to work%n", "orekit-data-master.zip", "https://gitlab.orekit.org/orekit/orekit-data/-/archive/master/orekit-data-master.zip", home.getAbsolutePath());
                        System.exit(1);
                    }

                    DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
                    manager.addProvider(new DirectoryCrawler(orekitData));
                    Frame eme2000 = FramesFactory.getEME2000();
                    AbsoluteDate dateTLE = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
                    System.out.println("DATE :" + DATE);
                    System.out.println("dateTLE :" + dateTLE);
                    Orbit iniOrbit = new KeplerianOrbit(sma, ecc, inc, aop, raan, ano, PositionAngle.MEAN, eme2000, dateTLE, 3.986004418E14);
                    double initMass = ERGOL + dryMass;
                    SpacecraftState initialState = new SpacecraftState(iniOrbit, initMass);
                    System.out.println("iniOrbit:" + iniOrbit);
                    System.out.println("initMass :" + initMass);
                    AbsoluteDate apsideDateObj = new AbsoluteDate(APSIDE_DATE, TimeScalesFactory.getUTC());
                    AbsoluteDate manoeuverStartDate = dateTLE.shiftedBy(manoeuverRelativeDate);
                    System.out.println("manoeuverStartDate :" + manoeuverStartDate);
                    System.out.println("apsideDateObj :" + apsideDateObj);
                    if (!isEqualOrAfterWithTolerance(manoeuverStartDate, apsideDateObj, TIME_TOLERANCE_SECONDS)) {
                        throw new IllegalArgumentException("Initial date must be equal to or later than the apside date within the allowed tolerance.");
                    } else {
                        System.out.println("Initial date is equal to or after the apside date within tolerance.");
                    }
                    AbsoluteDate endHorizonDate = new AbsoluteDate(endDateString, TimeScalesFactory.getUTC());
                    AbsoluteDate manoeuverEndDate = manoeuverStartDate.shiftedBy(durationOfManoeuver);
                    OrbitType orbitType = OrbitType.CIRCULAR;
                    if (ecc > 0.01) {
                        orbitType = OrbitType.KEPLERIAN;
                    } else {
                        orbitType = OrbitType.CIRCULAR;
                    }
                    KeplerianPropagator propagator = new KeplerianPropagator(iniOrbit);

                    System.out.println("Initial orbit:");
                    System.out.println("Date:" + ((Orbit)iniOrbit).getDate());
                    System.out.println("SMA:" + ((Orbit)iniOrbit).getA() / 1000.0 + " km");
                    System.out.println("MASS:" + initialState.getMass() + " kg");
                    double mass_limite = Math.pow(10, -3);
                    if (initialState.getMass()<mass_limite) {
                        throw new IllegalArgumentException("La masse initiale ne peut pas être négligeable.");
                    }
                    System.out.println("Mean:" + FastMath.toDegrees((new KeplerianOrbit(iniOrbit)).getMeanAnomaly()) + " °");
                    System.out.println("RAAN:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)iniOrbit).getRightAscensionOfAscendingNode(), Math.PI)));
                    System.out.println("ECC:" + ((Orbit)iniOrbit).getE());
                    System.out.println("AOP:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)iniOrbit).getPerigeeArgument(), Math.PI)));
                    System.out.println("Mean Motion :" + 86400.0 * initialState.getOrbit().getKeplerianMeanMotion() / 6.283185307179586);
                    System.out.println("Incli :" + FastMath.toDegrees(initialState.getI()));
                    System.out.println("Fram :" + initialState.getFrame());
                    System.out.println("FramP :" + propagator.getFrame());
                    System.out.println("PV:" + initialState.getOrbit().getPVCoordinates());
                    System.out.println("Kep:" + initialState.getOrbit());
                    SpacecraftState finalState = null;
                    if (maneuverType.equals("Impulse")) {
                        //Vecteur impulsion
                        Vector3D directionImp = new Vector3D(DVx * DV, DVy * DV, DVz * DV);
                        //directionImp = new Vector3D(1, 0, 0).normalize().scalarMultiply(DV);
                        //Trigger date de maneuvre
                        EventDetector trigger = new DateDetector(manoeuverStartDate);
                        ImpulseManeuver maneuver = new ImpulseManeuver(trigger, directionImp, ISP);
                        //PVCoordinates newPV = new PVCoordinates(initialState.getPVCoordinates().getPosition(),
                                //initialState.getPVCoordinates().getVelocity().add(directionImp));

                        // mfinale avec tsiolkovski
                        //double g0 = 9.80665;  // g en m/s²
                        //double finalMass = initialState.getMass() * Math.exp(-DV / (ISP * g0));
                        double g0 = 9.80665;  // Accélération gravitationnelle (m/s²)
                        double finalMass = initialState.getMass() * Math.exp(-DV / (ISP * g0));
                        double DV_m_s = DV ;
                        finalMass = initialState.getMass() * Math.exp(-DV_m_s / (ISP * g0));
                        propagator = new KeplerianPropagator(iniOrbit);
                        propagator.addEventDetector(maneuver);

                        System.out.println("\nAfter the maneuver:");
                        finalState = propagator.propagate(manoeuverStartDate.shiftedBy(durationOfManoeuver));
                        finalState = new SpacecraftState(finalState.getOrbit(), finalMass);
                        // Etat orbital avant 2eme maneuvre
                        /*System.out.println("\nBefore the second maneuver:");
                        System.out.println("SMA: " + finalState.getA() / 1000.0 + " km");
                        System.out.println("Eccentricity: " + finalState.getE());
                        // 2eme maneuvre
                        DV = 0.58399 * 1000;  // Appliquer le ∆V pour la deuxième manœuvre à 0.58399 km/s (en m/s)

                    // Nouvelle pos et nouvelle velocity
                        Vector3D newPosition = finalState.getPVCoordinates().getPosition();  // Position du satellite
                        Vector3D newVelocity = finalState.getPVCoordinates().getVelocity();  // Vitesse du satellite

                         // Vecteur vitesse et pos
                        System.out.println("Velocity vector before second maneuver: " + newVelocity);
                        System.out.println("Position vector before second maneuver: " + newPosition);

                         // Vitesse tangentielle
                        double Vx = newVelocity.getX();
                        double Vy = newVelocity.getY();
                        double Vz = newVelocity.getZ();

                        // Si le prod scalaire est < 0, on adapte le signe
                        if (Vector3D.dotProduct(newPosition, newVelocity) < 0) {
                            Vx = -Vx;
                            Vy = -Vy;
                            Vz = -Vz;
                        }

                        //Vitesse corrigée
                        Vector3D correctedVelocity = new Vector3D(Vx, Vy, Vz);

                        // Direction tangentielle
                        Vector3D angularMomentum = Vector3D.crossProduct(newPosition, correctedVelocity);  // Moment angulaire (vecteur normal)
                        Vector3D tangentDirection = Vector3D.crossProduct(angularMomentum, newPosition);   // Tangente à l'orbite

                        // Direction Tangentielle normalisee
                        tangentDirection = tangentDirection.normalize().scalarMultiply(DV);

                        // Direction poussée
                        System.out.println("Corrected direction of impulsion vector (second maneuver): " + tangentDirection);

                        // Date déclenchement maneuvre
                        manoeuverStartDate = finalState.getDate().shiftedBy(3899.66);  // Date après la première manœuvre

                        // trigger
                        trigger = new DateDetector(manoeuverStartDate);
                        maneuver = new ImpulseManeuver(trigger, tangentDirection, ISP);

                        // masse initiale pour la 2eme maneuvre
                        double initialMassForSecondManeuver = finalState.getMass();

                        // masse finale après tsiolkovski
                        double finalMassForSecondManeuver = initialMassForSecondManeuver * Math.exp(-DV / (ISP * g0));

                        // propagateur pour la 2eme maneuvre
                        propagator = new KeplerianPropagator(finalState.getOrbit());
                        propagator.addEventDetector(maneuver);

                        // état après maneuvre
                        finalState = propagator.propagate(manoeuverStartDate.shiftedBy(durationOfManoeuver));

                        // maj état final
                        finalState = new SpacecraftState(finalState.getOrbit(), finalMassForSecondManeuver);*/

                        // paramètres orbitaux après la deuxième manœuvre
                        System.out.println("\nAfter the last maneuver:");
                        System.out.println("SMA: " + finalState.getA() / 1000.0 + " km");
                        System.out.println("Eccentricity: " + finalState.getE());
                        System.out.println("Final mass (after last maneuver): " + finalState.getMass());





                    }else if (maneuverType.equals("Continuous")){
                        //Idée faire une série de maneubre impulsionnelle par pas de temps
                        double thrustStep = 1.0;  // Durée de chaque itération (en secondes)
                        double isp = ISP;  // Impulsion spécifique en secondes
                        double g0 = 9.80665;  // Accélération gravitationnelle (m/s²)
                        double thrust = THRUST;  // Poussée en Newtons

                        // Masse initiale et masse à vide
                        double initialMass = initialState.getMass();  // Masse initiale du satellite (incluant l'ergol)
                        double currentMass = initialMass;
                        dryMass = DRYMASS;  // Masse à vide (sans ergol)

                        double totalDuration = durationOfManoeuver;

                        // etat courant
                        SpacecraftState currentState = initialState;

                        // petits incréments
                        while (currentState.getDate().compareTo(manoeuverEndDate) < 0) {

                            // Init DV
                            double DVxPerStep = (DVx * DV) / totalDuration;  // Petit Delta-V en X par pas de temps
                            double DVyPerStep = (DVy * DV) / totalDuration;  // Petit Delta-V en Y par pas de temps
                            double DVzPerStep = (DVz * DV) / totalDuration;  // Petit Delta-V en Z par pas de temps

                            Vector3D thrustPerStep = new Vector3D(DVxPerStep, DVyPerStep, DVzPerStep);
                            Vector3D updatedVelocity = currentState.getPVCoordinates().getVelocity().add(thrustPerStep);

                            // conso ergol
                            double ergolConsumedByStep = (thrust * thrustStep) / (isp * g0);
                            currentMass -= ergolConsumedByStep;

                            // si total ergol consomme : break
                            if (currentMass < dryMass) {
                                currentMass = dryMass;
                                break;  // Arrêter la manœuvre si la masse à vide est atteinte
                            }

                            // maj pos et coord
                            PVCoordinates updatedPV = new PVCoordinates(currentState.getPVCoordinates().getPosition(), updatedVelocity);
                            Orbit updatedOrbit = new KeplerianOrbit(updatedPV, currentState.getFrame(), currentState.getDate(), MU);

                            // nouvel etat
                            currentState = new SpacecraftState(updatedOrbit, currentMass);

                            // propagation
                            KeplerianPropagator updatedPropagator = new KeplerianPropagator(updatedOrbit);

                            // propa avec notre pas de calcul
                            currentState = updatedPropagator.propagate(currentState.getDate().shiftedBy(thrustStep));

                            // verif masse
                            System.out.println("Masse actuelle après mise à jour : " + currentMass + " kg");

                            /*


                            // Paramètres pour la manœuvre continue
double thrustStep = 1.0;  // Durée de chaque itération (en secondes)
double isp = ISP;  // Impulsion spécifique en secondes
double g0 = 9.80665;  // Accélération gravitationnelle (m/s²)
double thrust = THRUST;  // Poussée en Newtons

// Masse initiale et masse à vide
double initialMass = initialState.getMass();  // Masse initiale du satellite (incluant l'ergol)
double currentMass = initialMass;
dryMass = DRYMASS;  // Masse à vide (sans ergol)

double totalDuration = durationOfManoeuver;  // Durée totale de la manœuvre

// État initial du satellite
SpacecraftState currentState = initialState;

// Date de début et de fin de la manœuvre continue
AbsoluteDate manoeuverEndDate = manoeuverStartDate.shiftedBy(totalDuration);

// Boucle de propagation avec petits incréments de temps (thrustStep)
while (currentState.getDate().compareTo(manoeuverEndDate) < 0) {

    // Calcul du petit Delta-V à chaque pas de temps
    double DVxPerStep = (DVx * DV) / totalDuration;  // Petit Delta-V en X par pas de temps
    double DVyPerStep = (DVy * DV) / totalDuration;  // Petit Delta-V en Y par pas de temps
    double DVzPerStep = (DVz * DV) / totalDuration;  // Petit Delta-V en Z par pas de temps

    // Vecteur poussée à chaque pas de temps
    Vector3D thrustPerStep = new Vector3D(DVxPerStep, DVyPerStep, DVzPerStep);

    // Mettre à jour la vitesse du satellite en ajoutant la poussée
    Vector3D updatedVelocity = currentState.getPVCoordinates().getVelocity().add(thrustPerStep);

    // Calcul de la consommation d'ergol à chaque pas de temps
    double ergolConsumedByStep = (thrust * thrustStep) / (isp * g0);
    currentMass -= ergolConsumedByStep;

    // Si la masse atteint la masse à vide (plus d'ergol), on arrête la manœuvre
    if (currentMass < dryMass) {
        currentMass = dryMass;
        break;  // Arrêter la manœuvre si la masse à vide est atteinte
    }

    // Mettre à jour les coordonnées de position et de vitesse du satellite
    PVCoordinates updatedPV = new PVCoordinates(currentState.getPVCoordinates().getPosition(), updatedVelocity);
    Orbit updatedOrbit = new KeplerianOrbit(updatedPV, currentState.getFrame(), currentState.getDate(), MU);

    // Créer un nouvel état du satellite avec l'orbite mise à jour et la nouvelle masse
    currentState = new SpacecraftState(updatedOrbit, currentMass);

    // Créer un propagateur Keplerien pour l'état mis à jour
    KeplerianPropagator updatedPropagator = new KeplerianPropagator(updatedOrbit);

    // Propager l'orbite avec le pas de temps (thrustStep)
    currentState = updatedPropagator.propagate(currentState.getDate().shiftedBy(thrustStep));

    // Afficher la masse actuelle après mise à jour
    System.out.println("Masse actuelle après mise à jour : " + currentMass + " kg");
}

// Masse restante d'ergol
double ergolRestant = currentMass - dryMass;

// État final après la manœuvre continue
finalState = new SpacecraftState(currentState.getOrbit(), currentMass);

                             */
                        }

                        // m_ergol restante
                        double ergolRestant = currentMass - dryMass;

                        // etat final
                        finalState = new SpacecraftState(currentState.getOrbit(), currentMass);


                    } else {
                        System.out.println("Maneuver : NONE");
                    }

                        Orbit finalOrbit = new KeplerianOrbit(finalState.getOrbit());
                        System.out.println("Final orbit:");
                        System.out.println("Date:" + finalState.getDate());
                        System.out.println("SMA:" + finalState.getA() / 1000.0 + " km");
                        System.out.println("MeanAnomaly = " + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getMeanAnomaly(), Math.PI)) + " deg");
                        System.out.println("manoeuverStartDate:" + manoeuverStartDate + " .");
                        System.out.println("manoeuverEndDate:" + manoeuverEndDate + " .");
                        System.out.println("RAAN:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getRightAscensionOfAscendingNode(), Math.PI)));
                        System.out.println("ECC:" + ((Orbit)finalOrbit).getE());
                        System.out.println("AOP:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getPerigeeArgument(), Math.PI)));
                        System.out.println("Mean Motion :" + 86400.0 * finalState.getOrbit().getKeplerianMeanMotion() / 6.283185307179586);
                        System.out.println("Incli :" + FastMath.toDegrees(finalState.getI()));
                        System.out.println("state :" + propagator.getInitialState().getDate());
                        System.out.println("PV:" + finalState.getOrbit().getPVCoordinates());
                        System.out.println("Kep:" + finalState.getOrbit());
                        System.out.println("Quaternion: at" + finalState.getDate() + " " + finalState.getAttitude().getOrientation().getRotation().getQ0() + ", " + finalState.getAttitude().getOrientation().getRotation().getQ1() + ", " + finalState.getAttitude().getOrientation().getRotation().getQ2() + ", " + finalState.getAttitude().getOrientation().getRotation().getQ3());
                        double carburant_limite = Math.pow(10, -3);
                        if ((finalState.getMass() - dryMass) < carburant_limite || ERGOL<carburant_limite) {
                            throw new IllegalArgumentException("La masse finale ne peut pas être inférieure à la masse à vide.");
                       }

                        FileWriter writer = new FileWriter("Result.txt", true);
                        BufferedWriter bufferedWriter = new BufferedWriter(writer);
                        bufferedWriter.newLine();
                        bufferedWriter.write("Orbital parameters post-maneuver :");
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.6f", dryMass));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.6f", m0 - finalState.getMass()));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.6f", finalState.getMass() - dryMass));
                        bufferedWriter.newLine();
                        bufferedWriter.write(((Orbit)finalOrbit).getDate().toString());
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.6f", (finalState.getPVCoordinates().getPosition().getNorm() - 6378137.0) / 1000.0));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.12f", finalState.getA() / 1000.0));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.7f", finalState.getE()));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.4f", FastMath.toDegrees(finalState.getI())));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.4f", FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getRightAscensionOfAscendingNode(), Math.PI))));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.8f", FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getPerigeeArgument(), Math.PI))));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.8f", FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getMeanAnomaly(), Math.PI))));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.14f", 86400.0 * finalState.getOrbit().getKeplerianMeanMotion() / 6.283185307179586));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.6f", finalState.getKeplerianPeriod()));
                        bufferedWriter.newLine();
                        bufferedWriter.close();
                    writer = new FileWriter("LastManeuverDate.txt", true);
                    bufferedWriter = new BufferedWriter(writer);
                    bufferedWriter.newLine();
                    bufferedWriter.write("Date post - maneuvre");
                    bufferedWriter.newLine();
                    bufferedWriter.write(String.valueOf(manoeuverEndDate));
                    bufferedWriter.newLine();
                    bufferedWriter.close();
                        double end = (double)System.currentTimeMillis();
                        double duration = (end - start) / 1000.0;
                        System.out.println("Execution time:" + duration);


                } catch (OrekitException var92) {
                    OrekitException e = var92;
                    System.err.println(e.getLocalizedMessage());
                    System.exit(1);
                }
            } catch (IOException var93) {
                IOException e = var93;
                e.printStackTrace();
            }
        } catch (IOException var94) {
            IOException e = var94;
            e.printStackTrace();
        }
    }
    private static boolean isEqualOrAfterWithTolerance(AbsoluteDate dateToCheck, AbsoluteDate referenceDate, double toleranceSeconds) {
        // Calculate the difference: positive if dateToCheck is after referenceDate
        double difference = dateToCheck.durationFrom(referenceDate);

        // If difference is greater than or within the negative tolerance, consider it as valid
        return difference >= -toleranceSeconds;
    }
}
