package org.example;

import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.Gradient;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.orekit.forces.maneuvers.ConstantThrustManeuver;
import org.orekit.orbits.EquinoctialOrbit;
import org.orekit.orbits.FieldEquinoctialOrbit;
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
import java.util.Collections;
import java.util.List;

public class QLawWithConstraint {

    private final SpacecraftState initialSpacecraftState;

    private final List<ConstantThrustManeuver> sol;
    private final double                       maxThrust;
    private final double                       maximumManeuverRatioPerOrbit;
    private final double                       thrustEffectivity;

    private final double demiAxisCriterion;
    private final double eccentricityCriterion;
    private final double inclinationCriterion;
    private final double perigeeArgumentCriterion;
    private final double rightAscensionoftheAscendingNodeCriterion;

    private final double   ISP;
    private final double[] targetOrbitState     = new double[6];
    private final double[] targetOrbitKeplerian = new double[6];

    private final DSFactory factory;

    private ArrayList<SpacecraftState> currentSpacecraftStateList = new ArrayList<>();

    private final double g0 = Constants.G0_STANDARD_GRAVITY;

    // Qlaw constructor

    public QLawWithConstraint(final Orbit initialOrbit, final Orbit targetOrbit, final double maxThrust,
                              final double ISP) {
        this(initialOrbit, targetOrbit, maxThrust, ISP, 1000, 0., 1, 10e4, 1000, 1000, 1000, 1000);

    }

    public QLawWithConstraint(final Orbit initialOrbit, final Orbit targetOrbit, final double maxThrust,
                              final double ISP, final double massInit) {
        this(initialOrbit, targetOrbit, maxThrust, ISP, massInit, 0., 1, 10e4, 1000, 1000, 1000, 1000);

    }

    public QLawWithConstraint(final Orbit initialOrbit, final Orbit targetOrbit, final double maxThrust,
                              final double ISP, final double massInit, final double thrustEffectivity) {
        this(initialOrbit, targetOrbit, maxThrust, ISP, massInit, thrustEffectivity, 1, 10e4, 1000, 1000, 1000, 1000);

    }

    public QLawWithConstraint(final Orbit initialOrbit, final Orbit targetOrbit, final double maxThrust,
                              final double ISP, final double massInit, final double thrustEffectivity,
                              final double maximumManeuverRatioPerOrbit) {
        this(initialOrbit, targetOrbit, maxThrust, ISP, massInit, thrustEffectivity, maximumManeuverRatioPerOrbit, 10e4,
             1000, 1000, 1000, 1000);

    }

    public QLawWithConstraint(final Orbit initialOrbit, final Orbit targetOrbit, final double maxThrust,
                              final double ISP, final double massInit, final double thrustEffectivity,
                              final double maximumManeuverRatioPerOrbit, final double demiAxisCriterion) {
        this(initialOrbit, targetOrbit, maxThrust, ISP, massInit, thrustEffectivity, maximumManeuverRatioPerOrbit,
             demiAxisCriterion, 1000, 1000, 1000, 1000);

    }

    public QLawWithConstraint(final Orbit initialOrbit, final Orbit targetOrbit, final double maxThrust,
                              final double ISP, final double massInit, final double thrustEffectivity,
                              final double maximumManeuverRatioPerOrbit, final double demiAxisCriterion,
                              final double eccentricityCriterion) {
        this(initialOrbit, targetOrbit, maxThrust, ISP, massInit, thrustEffectivity, maximumManeuverRatioPerOrbit,
             demiAxisCriterion, eccentricityCriterion, 1000, 1000, 1000);

    }

    public QLawWithConstraint(final Orbit initialOrbit, final Orbit targetOrbit, final double maxThrust,
                              final double ISP, final double massInit, final double thrustEffectivity,
                              final double maximumManeuverRatioPerOrbit, final double demiAxisCriterion,
                              final double eccentricityCriterion, final double inclinationCriterion) {
        this(initialOrbit, targetOrbit, maxThrust, ISP, massInit, thrustEffectivity, maximumManeuverRatioPerOrbit,
             demiAxisCriterion, eccentricityCriterion, inclinationCriterion, 1000, 1000);

    }

    public QLawWithConstraint(final Orbit initialOrbit, final Orbit targetOrbit, final double maxThrust,
                              final double ISP, final double massInit, final double thrustEffectivity,
                              final double maximumManeuverRatioPerOrbit, final double demiAxisCriterion,
                              final double eccentricityCriterion, final double inclinationCriterion,
                              final double perigeeArgumentCriterion) {
        this(initialOrbit, targetOrbit, maxThrust, ISP, massInit, thrustEffectivity, maximumManeuverRatioPerOrbit,
             demiAxisCriterion, eccentricityCriterion, inclinationCriterion, perigeeArgumentCriterion, 1000);

    }

    public QLawWithConstraint(final Orbit initialOrbit, final Orbit targetOrbit, final double maxThrust,
                              final double ISP,
                              final double massInit, final double thrustEffectivity,
                              final double maximumManeuverRatioPerOrbit, final double demiAxisCriterion,
                              final double eccentricityCriterion, final double inclinationCriterion,
                              final double perigeeArgumentCriterion,
                              final double rightAscensionoftheAscendingNodeCriterion) {
        this.initialSpacecraftState = new SpacecraftState(initialOrbit, massInit);// Default mass

        this.maxThrust                                 = maxThrust;
        this.ISP                                       = ISP;
        this.sol                                       = new ArrayList<>();
        this.thrustEffectivity                         = thrustEffectivity;
        this.maximumManeuverRatioPerOrbit              = maximumManeuverRatioPerOrbit;
        this.demiAxisCriterion                         = demiAxisCriterion;
        this.eccentricityCriterion                     = eccentricityCriterion;
        this.inclinationCriterion                      = inclinationCriterion;
        this.perigeeArgumentCriterion                  = perigeeArgumentCriterion;
        this.rightAscensionoftheAscendingNodeCriterion = rightAscensionoftheAscendingNodeCriterion;
        OrbitType.EQUINOCTIAL.mapOrbitToArray(targetOrbit, PositionAngle.TRUE, this.targetOrbitState, null);
        OrbitType.KEPLERIAN.mapOrbitToArray(targetOrbit, PositionAngle.TRUE, this.targetOrbitKeplerian, null);
        factory = new DSFactory(6, 1);

    }

    // Qlaw solver, return the list of maneuver to apply in order to arrive to the targeted orbit.
    public Output solve() {
        final int nOrbitMax = 1500; // Nombre max d'orbites
        final int nbSub     = 180;
        // Initialisation
        sol.clear();
        currentSpacecraftStateList.clear();

        FieldSpacecraftState<Gradient> currentFieldSpacecraftState =
                createCurrentFieldSpacecraftState(factory, initialSpacecraftState);
        currentSpacecraftStateList.add(currentFieldSpacecraftState.toSpacecraftState());
        int     nOrbit               = 0;
        boolean ConvergenceCriterion = createConvergenceCriterion(currentFieldSpacecraftState);

        // Boucle  de calcul

        // Première Boucle sur les orbites. Pour le moment pas de critère d'arret sur les paramètres. Je fixe un nbre d'orbite Max
        while (nOrbit < nOrbitMax && ConvergenceCriterion) {

            // Calcul de la période de l'orbite. Je la discrétise en 180 points avec Timestep. et je calcul le DV sur cette durée.
            double period   = currentFieldSpacecraftState.getKeplerianPeriod().getValue();
            int    timestep = FastMath.toIntExact(FastMath.round(period / nbSub));
            double DVMax =
                    (maxThrust * timestep) / currentFieldSpacecraftState.getMass().getValue();

            ArrayList<Double>  QValue;
            ArrayList<Integer> IndexBestQValue;

            // Boucle Interne sur tous les timestep de l'orbite pour récuperer dQ/dV à chaque instant. Calculé une fois au début de chaque orbite uniquement car hypothèse poussée faible.
            QValue = getQOrbitValue(currentFieldSpacecraftState, period, timestep);

            // Récupère les meilleurs moments pour manoeuvrer en fonction de l'efficience et du poussée max par orbite. Et les range dans l'ordre chronologique
            IndexBestQValue = getBestManeuvers(QValue, nbSub);

            for (int i = 0; i < IndexBestQValue.size(); i++) {
                currentFieldSpacecraftState =
                        propagateToNextManeuver(currentFieldSpacecraftState, i, timestep, IndexBestQValue);

                final Gradient Q = createQLaw(currentFieldSpacecraftState);

                Vector3D DV = computeBestDirection(Q, currentFieldSpacecraftState);

                boolean DVCriterion = getDVCriterion(DV, DVMax);

                if (!DVCriterion) {
                    Vector3D DVTrue = computeBestMagnitude(DV, Q, DVMax);

                    currentFieldSpacecraftState = addImpulse(currentFieldSpacecraftState, DVTrue);

                    ConstantThrustManeuver newManeuver =
                            createConstantThrustManeuver(currentFieldSpacecraftState, timestep, DVTrue);

                    sol.add(newManeuver);
                    currentSpacecraftStateList.add(currentFieldSpacecraftState.toSpacecraftState());
                }

            }
            // shift ensuite jusqu'à la fin de l'orbite pour recommencer ensuite
            currentFieldSpacecraftState =
                    PropagateToNextOrbit(currentFieldSpacecraftState, IndexBestQValue, nbSub, timestep);

            currentSpacecraftStateList.add(currentFieldSpacecraftState.toSpacecraftState());

            ConvergenceCriterion = createConvergenceCriterion(currentFieldSpacecraftState);
            nOrbit += 1;

        }
        System.out.println(nOrbit);
        return new

                Output(sol);

    }

    private boolean getDVCriterion(final Vector3D DV, final double DVMax) {
        return DV.getNorm() < DVMax * 0.01;
    }

    private Vector3D computeBestMagnitude(final Vector3D DV, final Gradient Q, final double DVMax) {
        final double ratio = Q.getValue() / DV.getNorm();
        if (ratio >= 1) {
            return DV.normalize().scalarMultiply(DVMax);
        }
        else {
            return DV.normalize().scalarMultiply(ratio * DVMax);
        }
    }

    private boolean createConvergenceCriterion(final FieldSpacecraftState currentSpacecraftState) {
        double[] currentOrbitKeplerian = new double[6];
        OrbitType.KEPLERIAN.mapOrbitToArray(currentSpacecraftState.toSpacecraftState().getOrbit(), PositionAngle.TRUE,
                                            currentOrbitKeplerian, null);

        double da = FastMath.abs(currentOrbitKeplerian[0] - targetOrbitKeplerian[0]);
        double de = FastMath.abs(currentOrbitKeplerian[1] - targetOrbitKeplerian[1]);
        double di = FastMath.toDegrees(FastMath.abs(currentOrbitKeplerian[2] - targetOrbitKeplerian[2]));
        double dw = FastMath.toDegrees(FastMath.abs(currentOrbitKeplerian[3] - targetOrbitKeplerian[3]));
        double dW = FastMath.toDegrees(FastMath.abs(currentOrbitKeplerian[4] - targetOrbitKeplerian[4]));

        return (da > demiAxisCriterion || de > eccentricityCriterion || di > inclinationCriterion
                || dw > perigeeArgumentCriterion || dW > rightAscensionoftheAscendingNodeCriterion);

    }

    private FieldSpacecraftState<Gradient> PropagateToNextOrbit(
            FieldSpacecraftState<Gradient> currentFieldSpacecraftState, final ArrayList<Integer> IndexBestQValue,
            final int nbSub, final int timestep) {
        if (nbSub - IndexBestQValue.get(IndexBestQValue.size() - 1) > 0) {
            return currentFieldSpacecraftState.shiftedBy(
                    (nbSub - IndexBestQValue.get(IndexBestQValue.size() - 1)) * timestep);
        }
        else {
            return currentFieldSpacecraftState;
        }
    }

    private FieldSpacecraftState propagateToNextManeuver(FieldSpacecraftState currentFieldSpacecraftState, final int i,
                                                         final int timestep, final ArrayList<Integer> IndexBestQValue) {
        if (i == 0) {
            currentFieldSpacecraftState = currentFieldSpacecraftState.shiftedBy(IndexBestQValue.get(i) * timestep);
        }
        else {
            currentFieldSpacecraftState = currentFieldSpacecraftState.shiftedBy(
                    timestep * (IndexBestQValue.get(i) - IndexBestQValue.get(i - 1)));

        }
        return currentFieldSpacecraftState;
    }

    private ArrayList<Integer> getBestManeuvers(final ArrayList<Double> QValue, final int nbSub) {
        final ArrayList<Integer> IndexBestQValue = new ArrayList<>();
        final double             QMax            = Collections.max(QValue);
        int                      QBestIndex;
        for (int i = 0; i < FastMath.round(nbSub * maximumManeuverRatioPerOrbit); i++) {

            double QBestValue = Collections.max(QValue);
            if (QBestValue / QMax >= thrustEffectivity) {
                QBestIndex = QValue.indexOf(QBestValue);
                IndexBestQValue.add(QBestIndex);
                QValue.set(QBestIndex, 0.);

            }
            else {
                break;
            }

        }
        Collections.sort(IndexBestQValue);
        return IndexBestQValue;
    }

    private ArrayList<Double> getQOrbitValue(FieldSpacecraftState<Gradient> currentFieldSpacecraftState,
                                             final double period, final int timestep) {
        FieldSpacecraftState testFieldSpacecraftState =
                createCurrentFieldSpacecraftState(factory, currentFieldSpacecraftState.toSpacecraftState());
        final ArrayList<Double> QValue = new ArrayList<>();
        for (int t = 0; t < period; t = t + timestep) {

            // créer un nouveau Field
            //SpacecraftState pour tester les meilleurs DV à appliquer en calculant le Qdot à chaque point de l'orbite

            // Calcule Q et Find the best Direction (in cartesian Frame) for the maneuver.
            final Gradient Qtest = createQLaw(testFieldSpacecraftState);

            // Plus direction normalisée
            final Vector3D DVtest = computeBestDirection(Qtest, testFieldSpacecraftState);
            // DVTest.getNorm <=> dQ/DV
            QValue.add(DVtest.getNorm());
            // Calcule Qmax et Qmin, DV.getNorm() <=> dQ/DV ce que l'on recherche

            testFieldSpacecraftState = testFieldSpacecraftState.shiftedBy(timestep);
        }
        return QValue;
    }

    private ConstantThrustManeuver createConstantThrustManeuver(FieldSpacecraftState currentFieldSpacecraftState,
                                                                double timestep, Vector3D DV) {

        final double thrust = DV.getNorm() * currentFieldSpacecraftState.toSpacecraftState().getMass() / timestep;
        return new ConstantThrustManeuver(
                currentFieldSpacecraftState.toSpacecraftState().getDate().shiftedBy(-0.5 * timestep), timestep, thrust,
                ISP, DV.normalize());
    }

    private FieldSpacecraftState<Gradient> createCurrentFieldSpacecraftState(final DSFactory factory,
                                                                             final SpacecraftState initialState) {

        final Gradient a         = new Gradient(factory.variable(0, initialState.getA()));
        final Gradient ex        = new Gradient(factory.variable(1, initialState.getEquinoctialEx()));
        final Gradient ey        = new Gradient(factory.variable(2, initialState.getEquinoctialEy()));
        final Gradient hx        = new Gradient(factory.variable(3, initialState.getHx()));
        final Gradient hy        = new Gradient(factory.variable(4, initialState.getHy()));
        final Gradient lv        = new Gradient(factory.variable(5, initialState.getLv()));
        final Gradient one       = a.getField().getOne();
        final Gradient fieldMass = one.multiply(initialState.getMass());

        final FieldAbsoluteDate<Gradient> fieldDate = new FieldAbsoluteDate<>(one.getField(), initialState.getDate());

        final FieldEquinoctialOrbit<Gradient> orbit = new FieldEquinoctialOrbit<>(a, ex, ey, hx, hy, lv, PositionAngle.TRUE,
                                                                                  initialState.getFrame(), fieldDate,
                                                                                  one.multiply(initialState.getMu()));

        return new FieldSpacecraftState<>(orbit, fieldMass);

    }

    private Gradient createQLaw(final FieldSpacecraftState<Gradient> fieldSpacecraftState) {
        // Parameters name

        Gradient a   = fieldSpacecraftState.getA();
        Gradient ex  = fieldSpacecraftState.getEquinoctialEx();
        Gradient ey  = fieldSpacecraftState.getEquinoctialEy();
        Gradient hx  = fieldSpacecraftState.getHx();
        Gradient hy  = fieldSpacecraftState.getHy();
        Gradient lv  = fieldSpacecraftState.getLv();
        Gradient mu  = fieldSpacecraftState.getMu();
        Gradient one = a.getField().getOne();
        Gradient p   = a.multiply(one.subtract(fieldSpacecraftState.getE().pow(2)));

        //delta
        Gradient da  = a.subtract(targetOrbitState[0]);
        Gradient dex = ex.subtract(targetOrbitState[1]);
        Gradient dey = ey.subtract(targetOrbitState[2]);
        Gradient dhx = hx.subtract(targetOrbitState[3]);
        Gradient dhy = hy.subtract(targetOrbitState[4]);
        // Gradient dlv = lv.subtract(finalOrbitState[5]);

        // Constants
        Gradient k = one
                .multiply(100); // Slope of the exponential barrier of penalty Function around critical function
        Gradient rp     = a.multiply(one.subtract(fieldSpacecraftState.getE())); // Current periapsis Radius in meter
        Gradient RMIN   = one.multiply(6.8e6); // Lowest permitted periapsis radius ie Approximately Earth Radius in meter
        Gradient coeffP = k.multiply((one.subtract(rp.divide(RMIN))));
        Gradient P      = coeffP.exp(); // Q Law Penalty Function

        // Weight Coefficient
        Gradient Wp  = one; // Weight of Penalty Function
        Gradient Wa  = one; // Weight for targeting a final value
        Gradient Wex = one; // Weight for targeting ex final value
        Gradient Wey = one; // Weight for targeting ey final value
        Gradient Whx = one.multiply(0); // Weight for targeting hx final value
        Gradient Why = one.multiply(0); // Weight for targeting hy final value

        // Scaling function  Sa for convergence
        Gradient aMinusAT          = (a.subtract(targetOrbitState[0])).abs();
        Gradient aMinusATDivide3aT = aMinusAT.divide(one.multiply(targetOrbitState[0]).multiply(3));
        Gradient coeffPower4       = aMinusATDivide3aT.pow(4);
        Gradient Sa                = (one.add(coeffPower4)).sqrt();

        // useful Coefficient
        Gradient sqrtADivideMuMultiply2      = (a.divide(mu)).sqrt().multiply(2);
        Gradient sqrtPDivideMu               = (p.divide(mu)).sqrt();
        Gradient sqrtPDivideMuDivide2        = sqrtPDivideMu.multiply(0.5);
        Gradient NormExPlusEy                = ((ex.pow(2)).add((ey.pow(2)))).sqrt();
        Gradient EyPlusSqrtOneMinusExSquare  = ey.add(one.subtract((ex.pow(2)))).sqrt();
        Gradient ExPlusSqrtOneMinusEySquare  = ex.add(one.subtract(ey.pow(2))).sqrt();
        Gradient OnePlusHxSquarePlusHySquare = one.add(hx.pow(2)).add(hy.pow(2));
        Gradient aCoefficient                = (NormExPlusEy.add(one)).divide(one.subtract(NormExPlusEy));

        // Maximum rate of change of Equinoctial Parameters for unitary Thrust
        Gradient aDotMax  = a.multiply(sqrtADivideMuMultiply2).multiply((aCoefficient.sqrt()));
        Gradient exDotMax = sqrtPDivideMu.multiply(2); // Also Equal to eyDotMax
        Gradient hxDotMax =
                sqrtPDivideMuDivide2.multiply((OnePlusHxSquarePlusHySquare).divide(ExPlusSqrtOneMinusEySquare));
        Gradient hyDotMax =
                sqrtPDivideMuDivide2.multiply((OnePlusHxSquarePlusHySquare).divide(EyPlusSqrtOneMinusExSquare));

        // Q components
        Gradient qa  = Sa.multiply(Wa).multiply((da.divide(aDotMax)).pow(2));
        Gradient qex = Wex.multiply((dex.divide(exDotMax)).pow(2));
        Gradient qey = Wey.multiply((dey.divide(exDotMax)).pow(2));
        Gradient qhx = Whx.multiply((dhx.divide(hxDotMax)).pow(2));
        Gradient qhy = Why.multiply((dhy.divide(hyDotMax)).pow(2));

        return (one.add((Wp.multiply(P)))).multiply(qa.add(qex).add(qey).add(qhx).add(qhy));

    }

    private Vector3D computeBestDirection(final Gradient Q, final FieldSpacecraftState fieldSpacecraftState) {
        double[] gradQ = Q.getGradient();

        // Initialisation avant la boucle de fx fy et fz
        double fx = 0;
        double fy = 0;
        double fz = 0;
        // Calcul de la Jacobienne de L'orbite dans le repere cartesien
        Gradient[][] jacobianCartesian = new Gradient[6][6];
        fieldSpacecraftState.getOrbit().getJacobianWrtCartesian(PositionAngle.TRUE, jacobianCartesian);
        // Calcul de fx fy et fz .
        for (int i = 0; i <= 5; i++) {
            fx += jacobianCartesian[i][3].getValue() * gradQ[i];
            fy += jacobianCartesian[i][4].getValue() * gradQ[i];
            fz += jacobianCartesian[i][5].getValue() * gradQ[i];

        }

        // Calcul de la direction optimale de poussée. Dans le repère cartésien
        double ux = -fx;
        double uy = -fy;
        double uz = -fz;
        // Plus normalisée <=> dQ/dV quand je prends la norme.
        return new Vector3D(ux, uy, uz);

    }

    private FieldSpacecraftState<Gradient> addImpulse(final FieldSpacecraftState<Gradient> fieldSpacecraftState,
                                                      final Vector3D impulseVector) {

        // Update velocity
        final Vector3D position                = fieldSpacecraftState.getPVCoordinates().getPosition().toVector3D();
        final Vector3D velocityWithoutManeuver = fieldSpacecraftState.getPVCoordinates().getVelocity().toVector3D();

        // Applique le DV dans le repère inertiel
        final Vector3D velocityWithManeuver = velocityWithoutManeuver.add(impulseVector);

        // Calcule les nouvelles coordonnées de l'orbite avec l'incrément de Vitesse (même position mais vitesse = V+DV)
        final PVCoordinates pvCoordinates = new PVCoordinates(position, velocityWithManeuver);
        final EquinoctialOrbit orbitWithManeuver =
                new EquinoctialOrbit(pvCoordinates, fieldSpacecraftState.getFrame(), fieldSpacecraftState.getDate().
                                                                                                         toAbsoluteDate(),
                                     fieldSpacecraftState.getMu().getValue());
        final double newMass =
                fieldSpacecraftState.getMass().getValue() * (FastMath.exp(-(impulseVector.getNorm()) / (g0 * ISP)));
        final SpacecraftState currentState = new SpacecraftState(orbitWithManeuver, newMass);

        return createCurrentFieldSpacecraftState(factory, currentState);
    }

    private double getL2Error() {
        // Last Spacecraft State
        int index = currentSpacecraftStateList.size() - 1;
        // Calculus of error between last State and Target State , rescaled with the value of the Last State.

        double da =
                FastMath.pow(
                        (currentSpacecraftStateList.get(index).getA() - targetOrbitState[0]) / (initialSpacecraftState.getA()
                                - targetOrbitState[0]),
                        2);
        double dex = FastMath.pow(
                (currentSpacecraftStateList.get(index).getEquinoctialEx() - targetOrbitState[1]) / (
                        initialSpacecraftState.getEquinoctialEx() - targetOrbitState[1]), 2);
        double dey = FastMath.pow(
                (currentSpacecraftStateList.get(index).getEquinoctialEy() - targetOrbitState[2]) / (
                        initialSpacecraftState.getEquinoctialEy() - targetOrbitState[2]), 2);
        double dhx =
                FastMath.pow((currentSpacecraftStateList.get(index).getHx() - targetOrbitState[3]) / (
                                     initialSpacecraftState.getHx() - targetOrbitState[3]),
                             2);
        double dhy =
                FastMath.pow((currentSpacecraftStateList.get(index).getHy() - targetOrbitState[4]) / (
                                     initialSpacecraftState.getHy() - targetOrbitState[4]),
                             2);
        // Still no convegence criterion for the true longitude
        //   double dlv= FastMath.pow(currentSpacecraftStateList.get(index).getLv()-finalOrbitState[5]/finalOrbitState[5],2);

        return FastMath.sqrt(da + dex + dey + dhx + dhy);
    }

    class Output {
        final List<ConstantThrustManeuver> maneuvers;

        private Output(final List<ConstantThrustManeuver> maneuvers) {
            this.maneuvers = maneuvers;
        }

        public List<ConstantThrustManeuver> getManeuvers() {
            return maneuvers;
        }

        public double getDVTotal() {
            //                        double Vfinal = currentSpacecraftStateList.get(currentSpacecraftStateList.size()-1).getPVCoordinates().getVelocity().getNorm();
            //                        double Vinit = currentSpacecraftStateList.get(0).getPVCoordinates().getVelocity().getNorm();
            //
            //                        return FastMath.abs(Vfinal-Vinit);

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
            return getEndDate().durationFrom(initialSpacecraftState.getDate());
        }

        public ArrayList<Double> getMassEvolution() {
            ArrayList<Double> massList = new ArrayList<>();
            for (SpacecraftState spacecraftState : currentSpacecraftStateList) {
                massList.add(spacecraftState.getMass());

            }

            return massList;
        }

        public double getDeltaM() {
            return currentSpacecraftStateList.get(0).getMass() - currentSpacecraftStateList.get(
                    currentSpacecraftStateList.size() - 1).getMass();
        }

        public List<SpacecraftState> getSpacecraftEvolution() {

            return currentSpacecraftStateList;
        }

        public Orbit getOrbitFinalState() {
            Orbit orbitFinalState = currentSpacecraftStateList.get(currentSpacecraftStateList.size() - 1).getOrbit();
            return OrbitType.KEPLERIAN.convertType(orbitFinalState);

        }

        public ArrayList<Double> getL2ErrorEvolution() {
            ArrayList<Double> L2ErrorEvolution = new ArrayList<>();

            for (SpacecraftState spacecraftState : currentSpacecraftStateList) {

                // Calculus of error between last State and Target State , rescaled with the value of the Last State.
                double da = FastMath.pow(
                        (spacecraftState.getA() - targetOrbitState[0]) / (initialSpacecraftState.getA()
                                - targetOrbitState[0]),
                        2);
                double dex = FastMath.pow(
                        (spacecraftState.getEquinoctialEx() - targetOrbitState[1]) / (
                                initialSpacecraftState.getEquinoctialEx() - targetOrbitState[1]),
                        2);
                double dey = FastMath.pow(
                        (spacecraftState.getEquinoctialEy() - targetOrbitState[2]) / (
                                initialSpacecraftState.getEquinoctialEy() - targetOrbitState[2]),
                        2);
                double dhx = FastMath.pow(
                        (spacecraftState.getHx() - targetOrbitState[3]) / (initialSpacecraftState.getHx()
                                - targetOrbitState[3]), 2);
                double dhy = FastMath.pow(
                        (spacecraftState.getHy() - targetOrbitState[4]) / (initialSpacecraftState.getHy()
                                - targetOrbitState[4]), 2);
                // Still no convegence criterion for the true longitude
                //   double dlv= FastMath.pow(currentSpacecraftStateList.get(index).getLv()-finalOrbitState[5]/finalOrbitState[5],2);
                L2ErrorEvolution.add(FastMath.sqrt(da + dex + dey + dhx + dhy));

            }

            return L2ErrorEvolution;
        }

        public double getL2Error() {

            return QLawWithConstraint.this.getL2Error();
        }

    }

}







