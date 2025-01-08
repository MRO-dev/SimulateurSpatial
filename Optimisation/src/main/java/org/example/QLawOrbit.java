package org.example;

import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.Gradient;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.orekit.forces.maneuvers.ConstantThrustManeuver;
import org.orekit.orbits.EquinoctialOrbit;
import org.orekit.orbits.FieldEquinoctialOrbit;
import org.orekit.orbits.FieldOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.FieldSpacecraftState;
import org.orekit.propagation.SpacecraftState;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.FieldAbsoluteDate;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;

import java.util.ArrayList;
import java.util.List;

public class QLawOrbit {

    private final Orbit initialOrbit;

    private final List<ConstantThrustManeuver> sol;
    private final double                       maxThrust;
    private final double                       ISP;
    private final double[]                     finalOrbitState = new double[6];
    private final int                          norbitMax       = 20; // Nombre max d'orbites
    private final int                          nb_sub          = 180;

    private final DSFactory factory;
    // Nombre de points par orbite pour le calcul de la Qlaw

    private double mass;

    private ArrayList<SpacecraftState> currentSpacecraftStateList = new ArrayList<>();

    private final double g0 = Constants.G0_STANDARD_GRAVITY;

    // Qlaw constructor

    public QLawOrbit(final Orbit initialOrbit, final Orbit finalOrbit, final double maxThrust, final double ISP) {
        this.initialOrbit = initialOrbit;
        this.maxThrust    = maxThrust;
        this.ISP          = ISP;
        this.mass         = 1000; // Default mass
        this.sol          = new ArrayList<>();
        OrbitType.EQUINOCTIAL.mapOrbitToArray(finalOrbit, PositionAngle.TRUE, this.finalOrbitState, null);
        factory = new DSFactory(6, 1);

    }

//        public QLaw(final Orbit initialOrbit, final Orbit finalOrbit, final double maxThrust, final double ISP,
//                    final double isp, final double massInit) {
//
//            new QLaw(initialOrbit, finalOrbit, maxThrust, ISP);
//            this.mass = massInit;
//
//        }

    // Qlaw solver, return the list of maneuver to apply in order to arrive to the targeted orbit.
    public Output solve() {
        // Initialisation
        sol.clear();
        currentSpacecraftStateList.clear();
        int                  maneuverCounter   = 0;
        FieldOrbit<Gradient> currentFieldOrbit = createCurrentFieldOrbit(factory, initialOrbit);
        Gradient fieldMass =
                currentFieldOrbit.getA().getField().getOne().multiply(mass);
        FieldSpacecraftState<Gradient> currentState = new FieldSpacecraftState(currentFieldOrbit, fieldMass);
        currentSpacecraftStateList.add(currentState.toSpacecraftState());

        // Convert initial orbit to equinoctial orbit in GCRF Frame -->  Choix de cette Frame car Inertielle
        // TO DO

        // Convert reference orbit to field reference orbit and create a FieldSpacecraftState attached to this orbit

        // Boucle  de calcul

        // Première Boucle sur les orbites. Pour le moment pas de critère d'arret sur les paramètres. Je fixe un nbre d'orbite Max
        for (int nOrbit = 0; nOrbit < norbitMax; nOrbit++) {

            // Calcul de la période de l'orbite. Je la discrétise en 180 points avec Timestep. et je calcul le DV sur cette durée.
            double period   = currentFieldOrbit.getKeplerianPeriod().getValue();
            int    timestep = FastMath.toIntExact(FastMath.round(period / nb_sub));

            // Print initialisation
            //  System.out.println(
            //        "\nOrbite numéro : " + norbit + "\n Periode :" + period + " s \n timestep : " + timestep);

            // Boucle interne sur tous les timestep de l'orbite
            for (int t = 0; t < period; t = t + timestep) {

                //                if (t != 0) {
                //                    fieldMass    =
                //                            currentFieldOrbit.getA().getField().getOne().multiply(mass);
                //                    currentState = new FieldSpacecraftState(currentFieldOrbit, fieldMass);
                //                    // Ajoute l'état transitoire du SpacecraftState
                //                    currentSpacecraftStateList.add(currentState.toSpacecraftState());
                //                }

                // Calcule Q et Find the best Direction (in cartesian Frame) for the maneuver.
                final Gradient Q         = createQLaw(currentFieldOrbit);
                final Vector3D DVUnitary = computeBestDirection(Q, currentFieldOrbit);

                // Problème ICI Définir correctement le critère de CV pour savoir si on pousse ou non ?
                // Pour le moment je choisis des coeff pour que ça pousse tout le temps
                // --> Ce qui correspond au critère "Minimiser la durée de transfert" en poussant tjr dans la meilleure direction
                final double Qdot      = 0;
                final double rendement = 1;
                if (Qdot < rendement) {

                    //                    // Problème de repère?
                                        currentFieldOrbit = createCurrentFieldOrbit(factory, addImpulse(currentState, DVUnitary));
                    System.out.println(currentFieldOrbit.getADot());
                    //
                    //                    //Ajoute la constante thrust maneuver à la liste début à la current date - 0.5*timestep (centrage manoeuvre), durée du timestep, DvXYZ à mettre en unitaire + TNW?
                                        sol.add(new ConstantThrustManeuver(currentFieldOrbit.toOrbit().getDate().shiftedBy(-0.5 * timestep),
                                                                           timestep, maxThrust, ISP,
                                                                           DVUnitary)); // a changer pour avoir dans repere TNW
                    //
                    //                    // Modifie la masse après la manoeuvre
                                       maneuverCounter += 1;
                                        mass = newMass(currentSpacecraftStateList, sol, maneuverCounter);

                    //   System.out.println("\n orbite numéro" + norbit + "\n time" + timestep + " Norme de DV" + DvXYZ.getNorm()
                    //                            + "Norme DVMAX" + DVmax + "Nouvelle masse" + mass);
                }

                  currentFieldOrbit = currentFieldOrbit.shiftedBy(timestep);

                //  Print resultats pour chaque timestep
                //   System.out.println("\n orbite numéro" + norbit + "\n time" + timestep + " Norme de DV" +DvXYZ.getNorm());

                //                System.out.println(
                //                        "\n \n \n Paramètres \n Demi grand axe A : " + currentFieldOrbit.getA().getValue() + " ATarget :"
                //                                + finalOrbit.getA() + "\n Ex  :"
                //                                + currentFieldOrbit.getEquinoctialEx().getValue() + " ExTarget : "
                //                                + finalOrbit.getEquinoctialEx() + "\n Ey :"
                //                                + currentFieldOrbit.getEquinoctialEy().getValue() + " EyTarget : "
                //                                + finalOrbit.getEquinoctialEy() + " \n Hx : " + currentFieldOrbit.getHx().getValue()
                //                                + " HxTarget : " + finalOrbit.getHx() + " \n Hy : " + currentFieldOrbit.getHx().getValue()
                //                                + " HyTarget : " + finalOrbit.getHy() + " \n Longitude Vraie : "
                //                                + currentFieldOrbit.getLv().getValue());

            }

        }

        return new Output(sol);

    }

    // Step 2 : Create the equinoctial Q-Law function for 5 parameters (except L)

    // Step 3 : Calculate the proximity coefficient of the Q-Law and its first derivative

    // Step 4 : Find alpha Beta (thrust angle) for which GradQ is maximum

    // Step 5 : Apply this thrust angle to our SpaceCraft and propagate till the next timestep

    // Step 6 : iterate

    private FieldOrbit<Gradient> createCurrentFieldOrbit(final DSFactory factory, final Orbit currentOrbit) {

        final Gradient a   = new Gradient(factory.variable(0, currentOrbit.getA()));
        final Gradient ex  = new Gradient(factory.variable(1, currentOrbit.getEquinoctialEx()));
        final Gradient ey  = new Gradient(factory.variable(2, currentOrbit.getEquinoctialEy()));
        final Gradient hx  = new Gradient(factory.variable(3, currentOrbit.getHx()));
        final Gradient hy  = new Gradient(factory.variable(4, currentOrbit.getHy()));
        final Gradient lv  = new Gradient(factory.variable(5, currentOrbit.getLv()));
        final Gradient one = a.getField().getOne();

        // Creation du field pour AbsoluteDate /!!\ Pas sur ! a faire verifier
        final FieldAbsoluteDate<Gradient> fieldDate = new FieldAbsoluteDate<>(one.getField(), currentOrbit.getDate());

        return new FieldEquinoctialOrbit<>(a, ex, ey, hx, hy, lv, PositionAngle.TRUE,
                                           currentOrbit.getFrame(), fieldDate,
                                           one.multiply(currentOrbit.getMu()));

    }

    private double newMass(final ArrayList<SpacecraftState> currentSpacecraftStateList,
                           final List<ConstantThrustManeuver> sol, final int maneuverCounter) {
        return currentSpacecraftStateList.get(maneuverCounter - 1).getMass()
                + sol.get(maneuverCounter - 1).getFlowRate() * sol.get(maneuverCounter - 1).getDuration();

    }

    private Gradient createQLaw(final FieldOrbit<Gradient> currentOrbit) {
        // Parameters name

        Gradient a   = currentOrbit.getA();
        Gradient ex  = currentOrbit.getEquinoctialEx();
        Gradient ey  = currentOrbit.getEquinoctialEy();
        Gradient hx  = currentOrbit.getHx();
        Gradient hy  = currentOrbit.getHy();
        Gradient lv  = currentOrbit.getLv();
        Gradient mu  = currentOrbit.getMu();
        Gradient one = a.getField().getOne();
        Gradient p   = a.multiply(one.subtract(currentOrbit.getE().pow(2)));

        //delta
        Gradient da  = a.subtract(finalOrbitState[0]);
        Gradient dex = ex.subtract(finalOrbitState[1]);
        Gradient dey = ey.subtract(finalOrbitState[2]);
        Gradient dhx = hx.subtract(finalOrbitState[3]);
        Gradient dhy = hy.subtract(finalOrbitState[4]);
        Gradient dlv = lv.subtract(finalOrbitState[5]);

        // Constants
        Gradient k = one
                .multiply(100); // Slope of the exponential barrier of penalty Function around critical function
        Gradient rp     = a.multiply(one.subtract(currentOrbit.getE())); // Current periapsis Radius in meter
        Gradient RMIN   = one.multiply(6.8e6); // Lowest permitted periapsis radius ie Approximately Earth Radius in meter
        Gradient coeffP = k.multiply((one.subtract(rp.divide(RMIN))));
        Gradient P      = coeffP.exp(); // Q Law Penalty Function

        // Weight Coefficient
        Gradient Wp  = one; // Weight of Penalty Function
        Gradient Wa  = one; // Weight for targeting a final value
        Gradient Wex = one; // Weight for targeting ex final value
        Gradient Wey = one; // Weight for targeting ey final value
        Gradient Whx = one; // Weight for targeting hx final value
        Gradient Why = one; // Weight for targeting hy final value

        // Scaling function  Sa for convergence
        Gradient aMinusAT          = (a.subtract(finalOrbitState[0])).abs();
        Gradient aMinusATDivide3aT = aMinusAT.divide(one.multiply(finalOrbitState[0]).multiply(3));
        Gradient coeffPower4       = aMinusATDivide3aT.pow(4);
        Gradient Sa                = (one.add(coeffPower4)).sqrt();

        // useful Coefficient
        Gradient sqrtADivideMu               = (a.divide(mu)).sqrt();
        Gradient sqrtPDivideMu               = (p.divide(mu)).sqrt();
        Gradient NormExPlusEy                = ((ex.pow(2)).add((ey.pow(2)))).sqrt();
        Gradient EyPlusSqrtOneMinusExSquare  = ey.add(one.subtract((ex.pow(2)))).sqrt();
        Gradient ExPlusSqrtOneMinusEySquare  = ex.add(one.subtract(ey.pow(2))).sqrt();
        Gradient OnePlusHxSquarePlusHySquare = one.add(hx.pow(2)).add(hy.pow(2));
        Gradient aCoefficient                = (NormExPlusEy.add(one)).divide(one.subtract(NormExPlusEy));

        // Maximum rate of change of Equinoctial Parameters for unitary Thrust
        Gradient aDotMax  = a.multiply(2).multiply(sqrtADivideMu).multiply((aCoefficient.sqrt()));
        Gradient exDotMax = sqrtPDivideMu.multiply(2); // Also Equal to eyDotMax
        Gradient hxDotMax =
                sqrtPDivideMu.multiply(0.5).multiply((OnePlusHxSquarePlusHySquare).divide(ExPlusSqrtOneMinusEySquare));
        Gradient hyDotMax =
                sqrtPDivideMu.multiply(0.5).multiply((OnePlusHxSquarePlusHySquare).divide(EyPlusSqrtOneMinusExSquare));

        // Q components
        Gradient qa  = Sa.multiply(Wa).multiply((da.divide(aDotMax)).pow(2));
        Gradient qex = Wex.multiply((dex.divide(exDotMax)).pow(2));
        Gradient qey = Wey.multiply((dey.divide(exDotMax)).pow(2));
        Gradient qhx = Whx.multiply((dhx.divide(hxDotMax)).pow(2));
        Gradient qhy = Why.multiply((dhy.divide(hyDotMax)).pow(2));

        Gradient Q = (one.add((Wp.multiply(P)))).multiply(qa.add(qex).add(qey).add(qhx).add(qhy));

        return Q;

    }

    private Vector3D computeBestDirection(final Gradient Q, final FieldOrbit FieldOrbit) {
        double[] gradQ = Q.getGradient();

        // Initialisation avant la boucle de fx fy et fz
        double fx = 0;
        double fy = 0;
        double fz = 0;
        // Calcul de la Jacobienne de L'orbite dans le repere cartesien
        Gradient[][] jacobianCartesian = new Gradient[6][6];
        FieldOrbit.getJacobianWrtCartesian(PositionAngle.TRUE, jacobianCartesian);
        // Calcul de fx fy et fz .
        for (int i = 0; i <= 5; i++) {
            fx += jacobianCartesian[i][3].getValue() * gradQ[i];
            fy += jacobianCartesian[i][4].getValue() * gradQ[i];
            fz += jacobianCartesian[i][5].getValue() * gradQ[i];

        }

        // Calcul de la direction optimale de poussée. Dans le repère cartésien? Pas sur de mon assertion
        double ux = -fx;
        double uy = -fy;
        double uz = -fz;

        return new Vector3D(ux, uy, uz).normalize();

    }

    private EquinoctialOrbit addImpulse(final FieldSpacecraftState<Gradient> state,
                                        final Vector3D impulseVector) {

        // Update velocity
        final Vector3D position                = state.getPVCoordinates().getPosition().toVector3D();
        final Vector3D velocityWithoutManeuver = state.getPVCoordinates().getVelocity().toVector3D();

        // Applique le DV dans le repère inertiel
        final Vector3D velocityWithManeuver = velocityWithoutManeuver.add(impulseVector);

        // Calcule les nouvelles coordonnées de l'orbite avec l'incrément de Vitesse (même position mais vitesse = V+DV)
        final PVCoordinates pvCoordinates = new PVCoordinates(position, velocityWithManeuver);
        final EquinoctialOrbit orbitWithManeuver =
                new EquinoctialOrbit(pvCoordinates, state.getFrame(), state.getDate().
                                                                           toAbsoluteDate(), state.getMu().getValue());

        return orbitWithManeuver;
    }

    private double getL2Error() {
        // Last Spacecraft State
        int index = currentSpacecraftStateList.size() - 1;
        // Calculus of error between last State and Target State , rescaled with the value of the Last State.
        double da =
                FastMath.pow((currentSpacecraftStateList.get(index).getA() - finalOrbitState[0]) / finalOrbitState[0],
                             2);
        double dex = FastMath.pow(
                (currentSpacecraftStateList.get(index).getEquinoctialEx() - finalOrbitState[1]) / finalOrbitState[1], 2);
        double dey = FastMath.pow(
                (currentSpacecraftStateList.get(index).getEquinoctialEy() - finalOrbitState[2]) / finalOrbitState[2], 2);
        double dhx =
                FastMath.pow((currentSpacecraftStateList.get(index).getHx() - finalOrbitState[3]) / finalOrbitState[3],
                             2);
        double dhy =
                FastMath.pow((currentSpacecraftStateList.get(index).getHy() - finalOrbitState[4]) / finalOrbitState[4],
                             2);
        // Still no convegence criterion for the true longitude
        //   double dlv= FastMath.pow(currentSpacecraftStateList.get(index).getLv()-finalOrbitState[5]/finalOrbitState[5],2);
        double L2Error = FastMath.sqrt(da + dex + dey + dhx + dhy);

        return L2Error;
    }

    class Output {
        final List<ConstantThrustManeuver> maneuvers;

        private Output(final List<ConstantThrustManeuver> maneuvers) {
            this.maneuvers = maneuvers;
        }

        // A modifier pour renvoyer un fichier csv avec les manoeuvres.
        public List<ConstantThrustManeuver> getManeuvers() {
            return maneuvers;
        }

        public double getDVTotal() {
            double DV = 0;
            for (int i = 0; i < maneuvers.size(); i++) {
                DV += (maneuvers.get(i).getThrustVector().getNorm() * maneuvers.get(i).getDuration())
                        / currentSpacecraftStateList.get(i).getMass();
            }
            return DV;
        }

        public AbsoluteDate getEndDate() {

            return maneuvers.get(maneuvers.size() - 1).getEndDate();
        }

        public double getTransfertTime() {
            return getEndDate().durationFrom(initialOrbit.getDate());
        }

        public ArrayList<Double> getMassEvolution() {
            ArrayList<Double> massList = new ArrayList<>();
            for (int i = 0; i < currentSpacecraftStateList.size(); i++) {
                massList.add(currentSpacecraftStateList.get(i).getMass());

            }

            return massList;
        }

        public double getDeltaM() {
            return getMassEvolution().get(0) - getMassEvolution().get(currentSpacecraftStateList.size() - 1);
        }

        // A modifier pour avoir en fichier csv.
        public List<SpacecraftState> getSpacecraftEvolution() {

            return currentSpacecraftStateList;
        }

        public ArrayList<Double> getL2ErrorEvolution() {
            ArrayList<Double> L2ErrorEvolution = new ArrayList<>();

            for (int index = 0; index < currentSpacecraftStateList.size(); index++) {

                // Calcuclus of error between last State and Target State , rescaled with the value of the Last State.
                double da = FastMath.pow(
                        (currentSpacecraftStateList.get(index).getA() - finalOrbitState[0]) / finalOrbitState[0], 2);
                double dex = FastMath.pow(
                        (currentSpacecraftStateList.get(index).getEquinoctialEx() - finalOrbitState[1]) / finalOrbitState[1],
                        2);
                double dey = FastMath.pow(
                        (currentSpacecraftStateList.get(index).getEquinoctialEy() - finalOrbitState[2]) / finalOrbitState[2],
                        2);
                double dhx = FastMath.pow(
                        (currentSpacecraftStateList.get(index).getHx() - finalOrbitState[3]) / finalOrbitState[3], 2);
                double dhy = FastMath.pow(
                        (currentSpacecraftStateList.get(index).getHy() - finalOrbitState[4]) / finalOrbitState[4], 2);
                // Still no convegence criterion for the true longitude
                //   double dlv= FastMath.pow(currentSpacecraftStateList.get(index).getLv()-finalOrbitState[5]/finalOrbitState[5],2);
                L2ErrorEvolution.add(FastMath.sqrt(da + dex + dey + dhx + dhy));

            }

            return L2ErrorEvolution;
        }

        public double getL2Error() {

            return getL2Error();
        }

    }

}







