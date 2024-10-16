package org.example.leoplaw.qlaw;

import org.example.constraint.Constraint;
import org.example.constraint.ConstraintsProcessor;
import org.example.constraint.ManeuverWindow;
import org.example.constraint.Maneuverability;
import org.example.leoplaw.LEOPLaw;
import org.example.leoplaw.tolerance.ToleranceBuilder;
import org.hipparchus.Field;
import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.Gradient;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.nonstiff.DormandPrince853FieldIntegrator;
import org.hipparchus.util.FastMath;
import org.orekit.forces.gravity.HolmesFeatherstoneAttractionModel;
import org.orekit.forces.gravity.potential.GravityFieldFactory;
import org.orekit.forces.gravity.potential.NormalizedSphericalHarmonicsProvider;
import org.orekit.forces.maneuvers.ConstantThrustManeuver;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.EquinoctialOrbit;
import org.orekit.orbits.FieldOrbit;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.FieldPropagator;
import org.orekit.propagation.FieldSpacecraftState;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.FieldKeplerianPropagator;
import org.orekit.propagation.numerical.FieldNumericalPropagator;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.Constants;
import org.orekit.utils.Fieldifier;
import org.orekit.utils.IERSConventions;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.TimeSpanMap;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.BiFunction;

public abstract class QLaw<T extends ToleranceBuilder.Tolerance> implements LEOPLaw {

    protected final SpacecraftState              initialSpacecraftState;
    protected final double[]                     targetOrbitKeplerian   = new double[6];
    protected final double[]                     targetOrbitEquinoctial = new double[6];
    protected final T                            tolerance;
    private final   Orbit                        targetOrbit;
    private final   List<ConstantThrustManeuver> sol                    = new ArrayList<>();
    private final   double                       maxThrust;
    private final   double                       maximumManeuverRatioPerOrbit;
    private final   double                       thrustEffectivity;
    private final   FieldPropagator<Gradient>    fieldPropagator;

    private final List<Constraint> constraintList;

    private final double                     ISP;
    private final DSFactory                  factory;
    private final double                     g0                         = Constants.G0_STANDARD_GRAVITY;
    private final ArrayList<SpacecraftState> currentSpacecraftStateList = new ArrayList<>();

    protected QLaw(final SpacecraftState initialSpacecraftState, final Orbit targetOrbit, final double maxThrust,
                   final double ISP, final double thrustEffectivity,
                   final double maximumManeuverRatioPerOrbit,
                   final List<Constraint> constraintList, final T tolerance,
                   final BiFunction<Field<Gradient>, SpacecraftState, FieldPropagator<Gradient>> fieldPropagator) {
        this.initialSpacecraftState = initialSpacecraftState;
        this.tolerance              = tolerance;
        this.targetOrbit            = targetOrbit;

        OrbitType.KEPLERIAN.mapOrbitToArray(targetOrbit, PositionAngle.TRUE,
                                            this.targetOrbitKeplerian, null);

        OrbitType.EQUINOCTIAL.mapOrbitToArray(targetOrbit, PositionAngle.TRUE,
                                              this.targetOrbitEquinoctial, null);

        this.maxThrust                    = maxThrust;
        this.ISP                          = ISP;
        this.thrustEffectivity            = thrustEffectivity;
        this.maximumManeuverRatioPerOrbit = maximumManeuverRatioPerOrbit;
        this.constraintList               = constraintList;

        factory = new DSFactory(6, 1);
        Gradient gradient = new Gradient(factory.constant(0));

        this.fieldPropagator = fieldPropagator.apply(gradient.getField(), initialSpacecraftState);

    }

    static protected FieldSpacecraftState<Gradient> createCurrentFieldSpacecraftState(final Field<Gradient> field,
                                                                                      final SpacecraftState initialState) {
        final FieldOrbit<Gradient> fieldOrbit = Fieldifier.fieldify(field, initialState.getOrbit());
        final Gradient             fieldMass  = field.getOne().multiply(initialState.getMass());

        return new FieldSpacecraftState<>(fieldOrbit, fieldMass);
    }

    // KeplerianPropagator
    public static FieldPropagator<Gradient> keplerianPropagator(Field<Gradient> field, SpacecraftState state) {
        FieldSpacecraftState<Gradient> fieldState = createCurrentFieldSpacecraftState(field, state);
        return new FieldKeplerianPropagator<>(fieldState.getOrbit(), fieldState.getMu());

    }

    public Output solve() {
        // Maximum number of orbit and subdivision of orbit to calculate potential maneuvers.
        final int nOrbitMax = 3000;
        final int nbSub     = 180;
        int       nOrbit    = 0;
        // Initialisation
        sol.clear();
        currentSpacecraftStateList.clear();

        FieldSpacecraftState<Gradient> currentFieldSpacecraftState =
                createCurrentFieldSpacecraftState(factory, initialSpacecraftState);
        fieldPropagator.resetInitialState(currentFieldSpacecraftState);

        currentSpacecraftStateList.add(currentFieldSpacecraftState.toSpacecraftState());

        //Criterion convergence
        boolean ConvergenceCriterion = createConvergenceCriterion(currentFieldSpacecraftState);

        // First Loop. Keep calculating best maneuver while number maximum of orbit is not reached or convergence criterion on keplerian parameters is not satisfied
        while (nOrbit < nOrbitMax && ConvergenceCriterion) {
            int index = 0;
            // Compute the period of the current orbit. and discretize this orbit in nbSub time interval (timestep).
            //  Compute the maximum possible DV during timestep
            double period   = currentFieldSpacecraftState.getKeplerianPeriod().getValue();
            int    timestep = FastMath.toIntExact(FastMath.round(period / nbSub));
            double DVMax =
                    (maxThrust * timestep) / currentFieldSpacecraftState.getMass().getValue();

            // Array to stock the Value of dQ/DV and their date in order to choose the best One.
            ArrayList<Double>  QValue;
            ArrayList<Integer> IndexBestQValue;

            // Loop over every timestep of the orbit in order to get dQ/dV at each timestep.
            // QValue calculated only at the beginning of each orbit because of low thrust hypothesis
            // we assume that even if we choose to thrust the thrust is low enough to consider that Qvalue didn't get so much perturbed.

            QValue = getQOrbitValue(currentFieldSpacecraftState, period, timestep);

            // Get the best moments (index) to compute a Maneuver considering the efficiency criterion and the maximum thrust amount per orbit
            // And sort these moments in a chronological order
            IndexBestQValue = getBestManeuvers(QValue, nbSub);

            // Loop to apply the chosen maneuvers to our SpacecraftState and eventually modulate the spacecraftState thrust.
            while (index < IndexBestQValue.size() && ConvergenceCriterion) {

                currentFieldSpacecraftState =
                        propagateToNextManeuver(currentFieldSpacecraftState, index, timestep, IndexBestQValue);

                final Gradient Q = createQLaw(currentFieldSpacecraftState);

                // Compute the best direction to thrust
                Vector3D DV = computeBestDirection(Q, currentFieldSpacecraftState);

                // If DV is too low we choose to ignore this maneuver
                boolean DVCriterion = getDVCriterion(DV, DVMax);
                // We recheck the convergence criterion after propagating to the next maneuver
                ConvergenceCriterion = createConvergenceCriterion(currentFieldSpacecraftState);

                // After choosing the best direction, we look for the best magnitude to apply
                // (DVmax if we are far to reach our target and the best fraction of DVmax otherwise)
                if (!DVCriterion && ConvergenceCriterion) {
                    Vector3D DVTrue = computeBestMagnitude(DV, Q, DVMax);

                    currentFieldSpacecraftState = addImpulse(currentFieldSpacecraftState, DVTrue);

                    ConstantThrustManeuver newManeuver =
                            createConstantThrustManeuver(currentFieldSpacecraftState, timestep, DVTrue);

                    // Store every maneuver and the evolution of the spacecraft state around the time.
                    sol.add(newManeuver);
                    currentSpacecraftStateList.add(currentFieldSpacecraftState.toSpacecraftState());

                    // Refresh the convergence criterion after each maneuver.
                    ConvergenceCriterion = createConvergenceCriterion(currentFieldSpacecraftState);
                }

                index += 1;
            }

            // After we made every possible maneuvers during this orbit, we shift our SpacecraftState to the next orbit.
            currentFieldSpacecraftState =
                    PropagateToNextOrbit(currentFieldSpacecraftState, IndexBestQValue, nbSub, timestep);

            currentSpacecraftStateList.add(currentFieldSpacecraftState.toSpacecraftState());

            ConvergenceCriterion = createConvergenceCriterion(currentFieldSpacecraftState);
            nOrbit += 1;

        }

        //  System.out.println("parametres de CV" + semiAxisCriterion +eccentricityCriterion + inclinationCriterion + perigeeArgumentCriterion + rightAscensionOfTheAscendingNodeCriterion);

        return new Output(sol);

    }

    // earth Numerical Propagator taking into account J2 effect
    public static FieldPropagator<Gradient> earthNumericalPropagatorJ2(Field<Gradient> field, SpacecraftState state) {
        FieldSpacecraftState<Gradient> fieldState = createCurrentFieldSpacecraftState(field, state);
        final double                   dP         = 0.001;
        final double                   minStep    = 0.001;
        final double                   maxStep    = 60;

        final double[][] tolerance = NumericalPropagator.tolerances(dP, state.getOrbit(), OrbitType.EQUINOCTIAL);

        DormandPrince853FieldIntegrator<Gradient> fieldIntegrator =
                new DormandPrince853FieldIntegrator<>(field, minStep, maxStep, tolerance[0], tolerance[1]);

        FieldNumericalPropagator<Gradient> numericalPropagator = new FieldNumericalPropagator<>(field, fieldIntegrator);

        NormalizedSphericalHarmonicsProvider potential = GravityFieldFactory.getConstantNormalizedProvider(2, 0);

        HolmesFeatherstoneAttractionModel J2model =
                new HolmesFeatherstoneAttractionModel(FramesFactory.getITRF(IERSConventions.IERS_2010, true),
                                                      potential);

        numericalPropagator.addForceModel(J2model);

        return numericalPropagator;
    }

    // Return 0 if there is no criterion convergence for demi axis , 1 otherwise

    abstract protected FieldSpacecraftState<Gradient> createCurrentFieldSpacecraftState(final DSFactory factory,
                                                                                        final SpacecraftState initialState);

    abstract protected Gradient createQLaw(final FieldSpacecraftState<Gradient> fieldSpacecraftState);

    abstract protected boolean createConvergenceCriterion(final FieldSpacecraftState<Gradient> currentSpacecraftState);

    // Compute the value of dQ/dV every timestep over one orbit.
    private ArrayList<Double> getQOrbitValue(FieldSpacecraftState<Gradient> currentFieldSpacecraftState,
                                             final double period, final int timestep) {
        final TimeSpanMap<ManeuverWindow> maneuverMap;
        final ArrayList<Double>           QValue = new ArrayList<>();

        FieldSpacecraftState<Gradient> testFieldSpacecraftState =
                createCurrentFieldSpacecraftState(factory, currentFieldSpacecraftState.toSpacecraftState());
        fieldPropagator.resetInitialState(testFieldSpacecraftState);

        if (!constraintList.isEmpty()) {
            maneuverMap =
                    new ConstraintsProcessor(constraintList).computeManeuverWindowMap(
                            testFieldSpacecraftState.getOrbit().toOrbit(),
                            testFieldSpacecraftState.getDate().toAbsoluteDate(),
                            testFieldSpacecraftState.getDate().shiftedBy(period).toAbsoluteDate());

            for (int t = 0; t < period; t = t + timestep) {

                // Calcule Q et Find the best Direction (in cartesian Frame) for the maneuver.
                // If the spacecraft state is in a "no-thrust slot" we force the value of dQ/Dv to 0 in order
                // to ensure that this slot won't be chosen.

                final Maneuverability maneuverability =
                        maneuverMap.get(testFieldSpacecraftState.getDate().toAbsoluteDate()).getManeuverability();

                if (maneuverability == Maneuverability.CANNOT) {
                    QValue.add(-1.);
                }
                else {
                    final Gradient Qtest = createQLaw(testFieldSpacecraftState);

                    // DVTest.getNorm <=> dQ/dV
                    final Vector3D DVtest = computeBestDirection(Qtest, testFieldSpacecraftState);

                    QValue.add(DVtest.getNorm());
                }
                testFieldSpacecraftState =
                        fieldPropagator.propagate(testFieldSpacecraftState.getDate().shiftedBy(timestep));
            }

        }
        else {
            for (int t = 0; t < period; t = t + timestep) {

                final Gradient Qtest = createQLaw(testFieldSpacecraftState);

                final Vector3D DVtest = computeBestDirection(Qtest, testFieldSpacecraftState);

                QValue.add(DVtest.getNorm());

                testFieldSpacecraftState = fieldPropagator.propagate(testFieldSpacecraftState.getDate().shiftedBy(timestep));
            }
        }

        return QValue;
    }

    // Compute the best direction to thrust to the SpacecraftState considering the value of Q
    private Vector3D computeBestDirection(final Gradient Q, final FieldSpacecraftState<Gradient> fieldSpacecraftState) {
        double[] gradQ = Q.getGradient();

        // Initialisation before the loop of fx fy and fz
        double fx = 0;
        double fy = 0;
        double fz = 0;

        // Compute the Cartesian Jacobian of the orbit.
        Gradient[][] jacobianCartesian = new Gradient[6][6];
        fieldSpacecraftState.getOrbit().getJacobianWrtCartesian(PositionAngle.TRUE, jacobianCartesian);

        // Calcul de fx fy et fz .
        for (int i = 0; i <= 5; i++) {
            fx += jacobianCartesian[i][3].getValue() * gradQ[i];
            fy += jacobianCartesian[i][4].getValue() * gradQ[i];
            fz += jacobianCartesian[i][5].getValue() * gradQ[i];

        }

        //Compute the optimal thrust direction in cartesian coordinate system
        double ux = -fx;
        double uy = -fy;
        double uz = -fz;

        return new Vector3D(ux, uy, uz);

    }

    // Compute the best moment to perform a maneuver.
    private ArrayList<Integer> getBestManeuvers(final ArrayList<Double> QValue, final int nbSub) {
        final ArrayList<Integer> IndexBestQValue = new ArrayList<>();
        double                   QMax            = Collections.max(QValue);

        // Case where the slot constraint is superior to the orbital period, to ensure that no manoeuvers will be taken
        // into account during this orbital period.
        if (QMax <= 0) {
            return IndexBestQValue;
        }
        else {
            int QBestIndex;
            for (int i = 0; i <= FastMath.round(nbSub * maximumManeuverRatioPerOrbit); i++) {

                double QBestValue = Collections.max(QValue);
                if (QBestValue / QMax >= thrustEffectivity) {
                    QBestIndex = QValue.indexOf(QBestValue);
                    IndexBestQValue.add(QBestIndex);
                    QValue.set(QBestIndex, -1.);

                }
                else {
                    break;
                }

            }
            Collections.sort(IndexBestQValue);
            return IndexBestQValue;
        }
    }

    // Propagate the currentSpacecraftState to the next admissible maneuver.
    private FieldSpacecraftState<Gradient> propagateToNextManeuver(
            FieldSpacecraftState<Gradient> currentFieldSpacecraftState, final int i,
            final int timestep, final ArrayList<Integer> IndexBestQValue) {
        fieldPropagator.resetInitialState(currentFieldSpacecraftState);
        if (i == 0) {
            currentFieldSpacecraftState = fieldPropagator.propagate(
                    currentFieldSpacecraftState.getDate().shiftedBy(IndexBestQValue.get(i) * timestep));
        }
        else {
            currentFieldSpacecraftState = fieldPropagator.propagate(currentFieldSpacecraftState.getDate().shiftedBy(
                    timestep * (IndexBestQValue.get(i) - IndexBestQValue.get(i - 1))));

        }
        return currentFieldSpacecraftState;
    }

    // Simple Function to control the value of DV. if DV is too low this function will allow us to ignore the associated maneuver.
    private boolean getDVCriterion(final Vector3D DV, final double DVMax) {

        return DV.getNorm() < DVMax * 0.01;
    }

    // Function to compute the best fraction of DV to apply in order to reach the target State.
    private Vector3D computeBestMagnitude(final Vector3D DV, final Gradient Q, final double DVMax) {
        final double ratio = Q.getValue() / DV.getNorm();
        if (ratio >= 1) {
            return DV.normalize().scalarMultiply(DVMax);
        }
        else {
            return DV.normalize().scalarMultiply(ratio * DVMax);
        }
    }

    // Add the impulse vector to our SpacecraftState
    private FieldSpacecraftState<Gradient> addImpulse(final FieldSpacecraftState<Gradient> fieldSpacecraftState,
                                                      final Vector3D impulseVector) {

        // Update velocity
        final Vector3D position                = fieldSpacecraftState.getPVCoordinates().getPosition().toVector3D();
        final Vector3D velocityWithoutManeuver = fieldSpacecraftState.getPVCoordinates().getVelocity().toVector3D();

        // Apply  DV in inertial coordinate system
        final Vector3D velocityWithManeuver = velocityWithoutManeuver.add(impulseVector);

        // Compute new coordinate of the orbit with the new velocity increment (same position but velocity = V+DV)
        final PVCoordinates pvCoordinates = new PVCoordinates(position, velocityWithManeuver);
        final EquinoctialOrbit orbitWithManeuver =
                new EquinoctialOrbit(pvCoordinates, fieldSpacecraftState.getFrame(), fieldSpacecraftState.getDate().
                                                                                                         toAbsoluteDate(),
                                     fieldSpacecraftState.getMu().getValue());

        // Update the mass of the spacecraftState using Tsiolkovsky equation.
        final double newMass =
                fieldSpacecraftState.getMass().getValue() * (FastMath.exp(-(impulseVector.getNorm()) / (g0 * ISP)));

        final SpacecraftState currentState = new SpacecraftState(orbitWithManeuver, newMass);

        return createCurrentFieldSpacecraftState(factory, currentState);
    }

    // Compute the ConstantThrustManeuver associated to our currentSpacecraftState and the DV we want to apply to it.
    private ConstantThrustManeuver createConstantThrustManeuver(FieldSpacecraftState<Gradient> currentFieldSpacecraftState,
                                                                double timestep, Vector3D DV) {

        final double thrust = DV.getNorm() * currentFieldSpacecraftState.toSpacecraftState().getMass() / timestep;

        // We center the maneuver this is why it is shifted by -0.5 * timestep.
        return new ConstantThrustManeuver(
                currentFieldSpacecraftState.toSpacecraftState().getDate().shiftedBy(-0.5 * timestep), timestep, thrust,
                ISP, DV.normalize());
    }

    // Propagate the currentSpacecraftState to the next orbit.
    private FieldSpacecraftState<Gradient> PropagateToNextOrbit(
            FieldSpacecraftState<Gradient> currentFieldSpacecraftState, final ArrayList<Integer> IndexBestQValue,
            final int nbSub, final int timestep) {
        fieldPropagator.resetInitialState(currentFieldSpacecraftState);

        // Case where no manoeuvers are made during one orbital period (due for example to a constrained slot manoeuver).
        if (IndexBestQValue.isEmpty()) {
            return fieldPropagator.propagate(currentFieldSpacecraftState.getDate().shiftedBy(nbSub * timestep));

        }
        else if (nbSub - IndexBestQValue.get(IndexBestQValue.size() - 1) > 0) {
            return fieldPropagator.propagate(
                    currentFieldSpacecraftState.getDate().shiftedBy(
                            nbSub - IndexBestQValue.get(IndexBestQValue.size() - 1) * timestep));

        }
        else {
            return currentFieldSpacecraftState;
        }
    }

    /** @return L2 error */
    // getL2Error of our spacecraftState
    private double getL2Error() {
        // Last Spacecraft State
        int index = currentSpacecraftStateList.size() - 1;
        // Calculus of error between last State and Target State , rescaled with the value of the Last State.
        double[] currentOrbitKeplerian = new double[6];
        OrbitType.KEPLERIAN.mapOrbitToArray(currentSpacecraftStateList.get(index).getOrbit(), PositionAngle.TRUE,
                                            currentOrbitKeplerian, null);

        // Set the value of da to 0 if there is no demi axis convergence criterion
        double da =
                FastMath.pow((currentOrbitKeplerian[0] - targetOrbitKeplerian[0]) / (initialSpacecraftState.getA()
                                     - targetOrbitKeplerian[0]),
                             2);

        // Set the value of de to 0 if there is no eccentricity convergence criterion
        double de = FastMath.pow(
                (currentOrbitKeplerian[1] - targetOrbitKeplerian[1]), 2);

        // Set the value of di to 0 if there is no inclination convergence criterion
        double di = FastMath.pow(
                (currentOrbitKeplerian[2] - targetOrbitKeplerian[2]) / FastMath.PI, 2);

        // Set the value of de to 0 if there is no perigee argument convergence criterion
        double dw =
                FastMath.pow((currentOrbitKeplerian[3] - targetOrbitKeplerian[3]) / FastMath.PI,
                             2);

        // Set the value of dW to 0 if there is no right ascension Argument convergence criterion
        double dW =
                FastMath.pow((currentOrbitKeplerian[4] - targetOrbitKeplerian[4]) / FastMath.PI,
                             2);

        // Still no convegence criterion for the true longitude
        //   double dlv= FastMath.pow(currentSpacecraftStateList.get(index).getLv()-finalOrbitState[5]/finalOrbitState[5],2);

        return FastMath.sqrt(da + de + di + dw + dW);
    }

    public class Output {
        final List<ConstantThrustManeuver> maneuvers;

        protected Output(final List<ConstantThrustManeuver> maneuvers) {
            this.maneuvers = maneuvers;
        }

        private double computeDeltaPhi(double phiTarget, double phiFinal) {
            // if one of the phi = 0 it returns statement 0.0 when using signum, so we change the value in order to have always a valid statement
            if (phiTarget == 0.0) {
                phiTarget = 0.0001;
            }
            if (phiFinal == 0.0) {
                phiFinal = 0.0001;
            }
            double deltaPhi;
            if (FastMath.signum(phiTarget) == FastMath.signum(phiFinal)) {
                deltaPhi = phiFinal - phiTarget;

            }
            // Case phiT negative and phiFinal Positive
            else if (FastMath.signum(phiTarget) == -1.0) {
                phiTarget = 2 * FastMath.PI + phiTarget;
                deltaPhi  = (phiFinal - phiTarget);

            }
            else {
                phiFinal = 2 * FastMath.PI + phiFinal;
                deltaPhi = phiFinal - phiTarget;
            }
            while (deltaPhi < 0) {
                deltaPhi = deltaPhi + 2 * FastMath.PI;
            }

            return deltaPhi;
        }

        // Compute the minimum AnomalyShift.
        private double computeMinimumAnomalyShift(final double initialPeriod, final double targetPeriod) {
            return initialPeriod * 2 * FastMath.PI / targetPeriod;
        }

        //Same method than before but for a lower orbit
        private List<ConstantThrustManeuver> bestPhasingMethodToDescend(double deltaPhi,
                                                                        final double targetMinimumAnomalyShift,
                                                                        final double targetPeriod,
                                                                        final KeplerianOrbit targetOrbitKeplerian,
                                                                        final double phiFinal) {
            double                       coastingTime;
            double                       aCoasting;
            AbsoluteDate                 coastingDate     = new AbsoluteDate();
            int                          indexNumber      = 0;
            double                       mu               = initialSpacecraftState.getMu();
            List<ConstantThrustManeuver> newManeuversList = new ArrayList<>(maneuvers);
            int                          revolutionNumber = 1;

            while (deltaPhi < targetMinimumAnomalyShift) {

                coastingTime = initialSpacecraftState.getKeplerianPeriod();

                newManeuversList = shiftConstantThrustManeuver(newManeuversList, coastingTime, 0);

                double phiTarget = targetOrbitKeplerian.shiftedBy(revolutionNumber * coastingTime).getTrueAnomaly();
                deltaPhi = computeDeltaPhi(phiTarget, phiFinal);
                revolutionNumber += 1;

            }
            if (revolutionNumber > 1) {
                System.out.println("Impossible de réaliser le phasage en une orbite, premier phasage effectué à t0");
                System.out.println("Nombre de révolution à l'orbite initiale" + revolutionNumber);
                System.out.println(" temps d'attente à l'orbite initiale de"
                                           + revolutionNumber * initialSpacecraftState.getKeplerianPeriod());
            }

            if (deltaPhi > targetMinimumAnomalyShift) {

                coastingTime = targetPeriod + (deltaPhi * targetPeriod) / (2 * FastMath.PI);
                aCoasting    = FastMath.pow(mu * FastMath.pow(coastingTime / (2 * FastMath.PI), 2), 1. / 3);
                System.out.println("Acoasting est =" + aCoasting);
                System.out.println("coasting time est de " + coastingTime);

                for (int i = 0; i < currentSpacecraftStateList.size() - 1; i++) {

                    // Check the moment where you cross the coasting orbit
                    if (FastMath.min(currentSpacecraftStateList.get(i).getA(), currentSpacecraftStateList.get(i + 1).getA())
                            < aCoasting
                            && FastMath.max(currentSpacecraftStateList.get(i).getA(),
                                            currentSpacecraftStateList.get(i + 1).getA()) > aCoasting) {

                        // Check the best moment to coast between the two spacecraftState.
                        if (FastMath.abs(currentSpacecraftStateList.get(i).getA() - aCoasting) >
                                FastMath.abs(currentSpacecraftStateList.get(i + 1).getA() - aCoasting)) {
                            coastingDate = currentSpacecraftStateList.get(i + 1).getDate();
                        }
                        else {
                            coastingDate = currentSpacecraftStateList.get(i).getDate();
                        }

                        break;
                    }

                }
                for (int i = 0; i < maneuvers.size(); i++) {
                    // Check the moment where the maneuvers are at the coastingDate
                    if (maneuvers.get(i).getStartDate().isBefore(coastingDate) &&
                            maneuvers.get(i + 1).getStartDate().isAfter(coastingDate)) {

                        // Choose the closest moment between the 2 maneuvers.
                        if (FastMath.abs(maneuvers.get(i).getStartDate().durationFrom(coastingDate)) > FastMath.abs(
                                maneuvers.get(i + 1).getStartDate().durationFrom(coastingDate))) {
                            indexNumber = i + 1;
                        }
                        else {
                            indexNumber = i;
                        }
                    }
                }
                newManeuversList =
                        shiftConstantThrustManeuver(newManeuversList, coastingTime, indexNumber);
            }

            return newManeuversList;

        }

        // Rewrite the maneuvers calculated with the Qlaw with the coasting time to apply to phase the final SpacecraftState
        private List<ConstantThrustManeuver> shiftConstantThrustManeuver(final List<ConstantThrustManeuver> maneuvers,
                                                                         final double orbitCoastingPeriod,
                                                                         final int indexNumber) {

            List<ConstantThrustManeuver> newConstantManeuverList = new ArrayList<>(maneuvers);

            // Index Number is the index where the coasting is supposed to start
            for (int i = indexNumber; i < newConstantManeuverList.size(); i++) {

                // Create the exact same maneuver than before but shifted by the orbitCoastingPeriod
                ConstantThrustManeuver newConstantThrustManeuver =
                        new ConstantThrustManeuver(maneuvers.get(i).getStartDate().shiftedBy(orbitCoastingPeriod),
                                                   maneuvers.get(i).getDuration(), maneuvers.get(i).getThrust(),
                                                   maneuvers.get(i).getISP(), maneuvers.get(i).getThrustVector());

                // Then replace the ith maneuver by the new maneuver shifted by the orbitCoasting
                newConstantManeuverList.set(i, newConstantThrustManeuver);

            }

            return newConstantManeuverList;

        }

        // Look for the best method to phase the SpacecraftState to an upper orbit
        private List<ConstantThrustManeuver> bestPhasingMethodToRise(double deltaPhi, final double targetMinimumAnomalyShift,
                                                                     final double targetPeriod,
                                                                     KeplerianOrbit targetOrbitKeplerian,
                                                                     final double phiFinal) {
            double                       coastingTime;
            double                       aCoasting;
            AbsoluteDate                 coastingDate     = new AbsoluteDate();
            int                          indexNumber      = 0;
            double                       mu               = initialSpacecraftState.getMu();
            List<ConstantThrustManeuver> newManeuversList = new ArrayList<>(maneuvers);
            int                          revolutionNumber = 1;

            while (deltaPhi < targetMinimumAnomalyShift) {

                coastingTime = initialSpacecraftState.getKeplerianPeriod();
                // coastingTime = deltaPhi / (targetOrbit.getKeplerianMeanMotion() - initialSpacecraftState.getOrbit()
                //                                                                                         .getKeplerianMeanMotion());
                newManeuversList = shiftConstantThrustManeuver(newManeuversList, coastingTime, 0);
                double phiTarget = targetOrbitKeplerian.shiftedBy(revolutionNumber * coastingTime).getTrueAnomaly();
                deltaPhi = computeDeltaPhi(phiTarget, phiFinal);
                revolutionNumber += 1;
            }

            if (revolutionNumber > 1) {
                System.out.println("Impossible de réaliser le phasage en une orbite, premier phasage effectué à t0");
                System.out.println("nombre de révolution à l'orbite initiale" + revolutionNumber);
            }

            if (deltaPhi > targetMinimumAnomalyShift) {

                coastingTime = deltaPhi * targetPeriod / (2 * FastMath.PI);
                aCoasting    = FastMath.pow(mu * FastMath.pow(coastingTime / (2 * FastMath.PI), 2), 1. / 3);
                System.out.println("Acoasting est =" + aCoasting);
                System.out.println("coasting time est de " + coastingTime);

                // Check the best moment to coast between the two spacecraftState.
                for (int i = 0; i < currentSpacecraftStateList.size() - 1; i++) {
                    if (FastMath.min(currentSpacecraftStateList.get(i).getA(), currentSpacecraftStateList.get(i + 1).getA())
                            < aCoasting
                            && FastMath.max(currentSpacecraftStateList.get(i).getA(),
                                            currentSpacecraftStateList.get(i + 1).getA()) > aCoasting) {
                        coastingDate = currentSpacecraftStateList.get(i).getDate();

                        // Check the best moment to coast between the two spacecraftState.
                        if (FastMath.abs(currentSpacecraftStateList.get(i).getA() - aCoasting) >
                                FastMath.abs(currentSpacecraftStateList.get(i + 1).getA() - aCoasting)) {
                            coastingDate = currentSpacecraftStateList.get(i + 1).getDate();
                        }
                        else {
                            coastingDate = currentSpacecraftStateList.get(i).getDate();
                        }

                        break;
                    }

                }
                for (int i = 0; i < maneuvers.size() - 1; i++) {
                    if (maneuvers.get(i).getStartDate().isBefore(coastingDate) && maneuvers.get(i + 1).getStartDate()
                                                                                           .isAfter(coastingDate)) {
                        // Choose the closest moment between the 2 maneuvers.
                        if (FastMath.abs(maneuvers.get(i).getStartDate().durationFrom(coastingDate)) > FastMath.abs(
                                maneuvers.get(i + 1).getStartDate().durationFrom(coastingDate))) {
                            indexNumber = i + 1;
                        }
                        else {
                            indexNumber = i;
                        }

                    }

                }

                newManeuversList =
                        shiftConstantThrustManeuver(newManeuversList, coastingTime, indexNumber);

            }

            return newManeuversList;

        }

        public List<ConstantThrustManeuver> getManeuvers() {
            return maneuvers;
        }

        // get the DVTotal for the Qlaw.
        public double getDVTotal() {

            return g0 * ISP * FastMath.log(currentSpacecraftStateList.get(0).getMass() / (currentSpacecraftStateList.get(
                    currentSpacecraftStateList.size() - 1).getMass()));
        }

        // Get the evolution of the spacecraftState mass during the Qlaw.
        public ArrayList<Double> getMassEvolution() {
            ArrayList<Double> massList = new ArrayList<>();
            for (SpacecraftState spacecraftState : currentSpacecraftStateList) {
                massList.add(spacecraftState.getMass());

            }

            return massList;
        }

        // get the Mass consumed during the Qlaw evolution.
        public double getDeltaM() {
            return initialSpacecraftState.getMass() - currentSpacecraftStateList.get(
                    currentSpacecraftStateList.size() - 1).getMass();
        }

        // Get the evolution of the spacecraftState during the Qlaw.
        public List<SpacecraftState> getSpacecraftEvolution() {

            return currentSpacecraftStateList;
        }

        // Get the evolution of the L2Error during the Qlaw
        public ArrayList<Double> getL2ErrorEvolution() {
            ArrayList<Double> L2ErrorEvolution = new ArrayList<>();

            for (SpacecraftState spacecraftState : currentSpacecraftStateList) {

                // Calculus of error between last State and Target State , rescaled with the value of the Last State.
                double[] currentOrbitKeplerian = new double[6];
                OrbitType.KEPLERIAN.mapOrbitToArray(spacecraftState.getOrbit(), PositionAngle.TRUE,
                                                    currentOrbitKeplerian, null);

                // Set the value of da to 0 if there is no demi axis convergence criterion
                double da =
                        FastMath.pow((currentOrbitKeplerian[0] - targetOrbitKeplerian[0]) / (initialSpacecraftState.getA()
                                             - targetOrbitKeplerian[0]),
                                     2);

                // Set the value of de to 0 if there is no eccentricity convergence criterion
                double de = FastMath.pow(
                        (currentOrbitKeplerian[1] - targetOrbitKeplerian[1]), 2);

                // Set the value of di to 0 if there is no inclination convergence criterion
                double di = FastMath.pow(
                        (currentOrbitKeplerian[2] - targetOrbitKeplerian[2]) / FastMath.PI, 2);

                //                // Set the value of de to 0 if there is no perigee argument convergence criterion
                double dw =
                        FastMath.pow((currentOrbitKeplerian[3] - targetOrbitKeplerian[3]) / (FastMath.PI / 2),
                                     2);

                //                // Set the value of dW to 0 if there is no right ascension Argument convergence criterion
                double dW =
                        FastMath.pow((currentOrbitKeplerian[4] - targetOrbitKeplerian[4]) / FastMath.PI,
                                     2);

                L2ErrorEvolution.add(FastMath.sqrt(da + de + di + dw + dW));

            }

            return L2ErrorEvolution;
        }

        // Get the number of revolution made during the Qlaw.
        public int getNumberofRevolution() {
            int revolutionNumber = 0;
            for (int i = 0; i < currentSpacecraftStateList.size() - 1; i++) {
                if (FastMath.toDegrees(currentSpacecraftStateList.get(i).getLv()) < initialSpacecraftState.getLv()
                        && FastMath.toDegrees(currentSpacecraftStateList.get(i + 1).getLv())
                        >= initialSpacecraftState.getLv()) {
                    revolutionNumber += 1;
                }
            }
            return revolutionNumber;
        }

        // Compute the new ManeuverList to apply to phase the final SpacecraftState.
        public List<ConstantThrustManeuver> getPhasingManeuvers() {

            // Compute deltaPhi between our Target Orbit and the obtained Orbit from the Qlaw.
            // Transform spacecraftState into keplerian Orbit. Mandatory to get true anomaly.
            KeplerianOrbit targetOrbitKeplerian =
                    new KeplerianOrbit(getTargetOrbitFinalState());
            KeplerianOrbit finalOrbitKeplerian =
                    new KeplerianOrbit(getOrbitFinalState());

            double phiFinal      = finalOrbitKeplerian.getTrueAnomaly();
            double phiTarget     = targetOrbitKeplerian.getTrueAnomaly();
            double deltaPhi;
            double initialPeriod = initialSpacecraftState.getKeplerianPeriod();
            double targetPeriod  = targetOrbit.getKeplerianPeriod();
            double targetMinimumAnomalyShift;

            // Case where we are looking for a full orbital period to coast, so we calculate phiF -phiT because it's up to orbit Target to "reach" the final one.
            if (FastMath.round(FastMath.toDegrees(phiFinal) - FastMath.toDegrees(phiTarget)) == 0) {
                System.out.println("No needed maneuvers");
                return maneuvers;
            }

            // Make the dinstinction between different True Anomaly value as phi is between [-Pi,Pi] in order to always get the shorter arc and hence the right deltaPhi

            // Not sure for this statement seems to work still perform test
            if (initialSpacecraftState.getA() > targetOrbitKeplerian.getA()) {
                deltaPhi                  = computeDeltaPhi(phiTarget, phiFinal);
                targetMinimumAnomalyShift = computeMinimumAnomalyShift(targetPeriod, initialPeriod);
                System.out.println("cas a>aT");
                return bestPhasingMethodToDescend(deltaPhi, targetMinimumAnomalyShift, targetPeriod, targetOrbitKeplerian,
                                                  phiFinal);
            }
            else {
                deltaPhi                  = computeDeltaPhi(phiTarget, phiFinal);
                targetMinimumAnomalyShift = computeMinimumAnomalyShift(initialPeriod, targetPeriod);
                System.out.println("Cas a < aT");
                return bestPhasingMethodToRise(deltaPhi, targetMinimumAnomalyShift, targetPeriod, targetOrbitKeplerian,
                                               phiFinal);
            }
            // Compute deltaPhiMin which is the minimum dephasage possible

            // At this orbit (initial one) you know you will perform one revolution in one period while the target true anomaly will shift for this amount

        }

        public Orbit getTargetOrbitFinalState() {
            return targetOrbit.shiftedBy(getTransferTime());

        }

        // get The transfer time
        public double getTransferTime() {

            return getEndDate().durationFrom(initialSpacecraftState.getDate());
        }

        // Get the endDate of the Qlaw (end of maneuvers).
        public AbsoluteDate getEndDate() {
            return maneuvers.get(maneuvers.size() - 1).getEndDate();
        }

        // Get the orbit reached after the end of the Qlaw.
        public Orbit getOrbitFinalState() {
            Orbit orbitFinalState = currentSpacecraftStateList.get(currentSpacecraftStateList.size() - 1).getOrbit();
            return OrbitType.KEPLERIAN.convertType(orbitFinalState);

        }

        // getL2Error()
        public double getL2Error() {

            return this.getL2Error();
        }
    }

}

