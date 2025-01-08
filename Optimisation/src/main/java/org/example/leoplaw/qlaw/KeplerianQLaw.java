package org.example.leoplaw.qlaw;

import org.example.constraint.Constraint;
import org.example.leoplaw.LEOPLaw;
import org.example.leoplaw.tolerance.ToleranceBuilder;
import org.hipparchus.Field;
import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.Gradient;
import org.hipparchus.util.FastMath;
import org.orekit.orbits.FieldKeplerianOrbit;
import org.orekit.orbits.FieldOrbit;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.FieldPropagator;
import org.orekit.propagation.FieldSpacecraftState;
import org.orekit.propagation.SpacecraftState;
import org.orekit.time.FieldAbsoluteDate;

import java.util.List;
import java.util.function.BiFunction;

public class KeplerianQLaw extends QLaw<ToleranceBuilder.KeplerianTolerance> implements LEOPLaw {

    // Qlaw constructor

    public KeplerianQLaw(final Orbit initialOrbit, final Orbit targetOrbit,
                         final double maxThrust,
                         final double ISP,
                         final double initialMass, final ToleranceBuilder.KeplerianTolerance tolerance,
                         final double thrustEffectivity,
                         final double maximumManeuverRatioPerOrbit,
                         final List<Constraint> constraintList,
                         final BiFunction<Field<Gradient>, SpacecraftState, FieldPropagator<Gradient>> fieldPropagator) {

        super(new SpacecraftState(initialOrbit, initialMass), targetOrbit, maxThrust, ISP, thrustEffectivity,
              maximumManeuverRatioPerOrbit, constraintList, tolerance, fieldPropagator);

    }

    // Create a FieldSpacecraftState from a DSfactory and a SpacecraftState
    protected FieldSpacecraftState<Gradient> createCurrentFieldSpacecraftState(final DSFactory factory,
                                                                               final SpacecraftState initialState) {
        double[] initialKeplerianOrbit = new double[6];

        OrbitType.KEPLERIAN.mapOrbitToArray(initialState.getOrbit(), PositionAngle.TRUE,
                                            initialKeplerianOrbit, null);

        final Gradient a         = new Gradient(factory.variable(0, initialState.getA()));
        final Gradient e         = new Gradient(factory.variable(1, initialState.getE()));
        final Gradient i         = new Gradient(factory.variable(2, initialState.getI()));
        final Gradient w         = new Gradient(factory.variable(3, initialKeplerianOrbit[3]));
        final Gradient W         = new Gradient(factory.variable(4, initialKeplerianOrbit[4]));
        final Gradient lv        = new Gradient(factory.variable(5, initialState.getLv()));
        final Gradient one       = a.getField().getOne();
        final Gradient fieldMass = one.multiply(initialState.getMass());

        final FieldAbsoluteDate<Gradient> fieldDate = new FieldAbsoluteDate<>(one.getField(), initialState.getDate());

        final FieldKeplerianOrbit<Gradient> orbit = new FieldKeplerianOrbit<>(a, e, i, w, W, lv, PositionAngle.TRUE,
                                                                              initialState.getFrame(), fieldDate,
                                                                              one.multiply(initialState.getMu()));

        return new FieldSpacecraftState<>(orbit, fieldMass);

    }

    // Compute the QLaw from Varga paper.
    protected Gradient createQLaw(final FieldSpacecraftState<Gradient> fieldSpacecraftState) {

        FieldOrbit<Gradient> fieldOrbit = fieldSpacecraftState.getOrbit();
        FieldKeplerianOrbit<Gradient> fieldKepOrbit =
                (FieldKeplerianOrbit<Gradient>) OrbitType.KEPLERIAN.convertType(fieldOrbit);

        // Parameters name

        Gradient a  = fieldKepOrbit.getA();
        Gradient e  = fieldKepOrbit.getE();
        Gradient i  = fieldKepOrbit.getI();
        Gradient w  = fieldKepOrbit.getPerigeeArgument();
        Gradient W  = fieldKepOrbit.getRightAscensionOfAscendingNode();
        Gradient mu = fieldSpacecraftState.getMu();
        w = w.divide(w.getValue()).multiply(w.getValue() % 2 * FastMath.PI);

        Gradient one = a.getField().getOne();
        Gradient p   = a.multiply(one.subtract(fieldSpacecraftState.getE().pow(2)));

        if (fieldSpacecraftState.getE().getValue() == 0) {
            p = a;
        }

        //delta
        Gradient da = a.subtract(targetOrbitKeplerian[0]);
        Gradient de = e.subtract(targetOrbitKeplerian[1]);
        Gradient di = i.subtract(targetOrbitKeplerian[2]);
        Gradient dw = circularDifference(w, targetOrbitKeplerian[3] % (2 * FastMath.PI));
        Gradient dW = circularDifference(W, targetOrbitKeplerian[4] % (2 * FastMath.PI));

        // Constants
        Gradient k = one
                .multiply(100); // Slope of the exponential barrier of penalty Function around critical function

        Gradient rp = a.multiply(one.subtract(fieldSpacecraftState.getE())); // Current periapsis Radius in meter
        if (fieldSpacecraftState.getE().getValue() == 0) {
            rp = a;
        }

        Gradient RMIN =
                one.multiply(6.8e6); // Lowest permitted periapsis radius ie Approximately Earth Radius in meter
        Gradient coeffP = k.multiply((one.subtract(rp.divide(RMIN))));
        Gradient P      = coeffP.exp(); // Q Law Penalty Function

        // Weight Coefficient
        Gradient Wp    = one; // Weight of Penalty Function
        Gradient Wa    = computeWaCoefficient(fieldSpacecraftState);
        Gradient We    = computeWeCoefficient(fieldSpacecraftState);
        Gradient Wi    = computeWiCoefficient(fieldSpacecraftState);
        Gradient Wpa   = computeWpaCoefficient(fieldSpacecraftState);
        Gradient Wraan = computeWraanCoefficient(fieldSpacecraftState);

        // Scaling function  Sa for convergence
        Gradient aMinusAT          = (a.subtract(targetOrbitEquinoctial[0])).abs();
        Gradient aMinusATDivide3aT = aMinusAT.divide(one.multiply(targetOrbitEquinoctial[0]).multiply(3));
        Gradient coeffPower4       = aMinusATDivide3aT.pow(4);
        Gradient Sa                = (one.add(coeffPower4)).sqrt();

        // useful Coefficient
        Gradient h    = p.multiply(mu).sqrt();
        Gradient sub1 = one.subtract(e.pow(2)).divide(e.pow(3).multiply(2.));
        Gradient sub2 = sub1.pow(2).add(1. / 27.).sqrt();
        // subtract au lieu de add?
        Gradient cosThetaxx = sub2.add(sub1).cbrt().subtract(sub2.subtract(sub1).cbrt()).subtract(e.pow(-1));
        Gradient sinThetaxx = one.subtract(cosThetaxx.pow(2)).sqrt();
        Gradient rxx        = p.divide(one.add(e.multiply(cosThetaxx)));

        // maximum parameter change rate
        Gradient axx = a.pow(3).divide(mu).multiply(one.add(e)).divide(one.subtract(e)).sqrt().multiply(2.);
        Gradient exx = p.divide(h).multiply(2.);
        Gradient ixx =
                p.divide(h.multiply(one.subtract(e.multiply(w.sin()).pow(2)).sqrt().subtract(e.multiply(w.cos().abs()))));
        Gradient wxx =
                p.multiply(cosThetaxx).pow(2).add(p.add(rxx).multiply(sinThetaxx).pow(2)).sqrt().divide(e.multiply(h));

        Gradient Wxx = p.divide(h.multiply(i.sin()).multiply(
                one.subtract(e.multiply(w.cos()).pow(2)).sqrt().subtract(e.multiply(w.sin().abs()))));

        // Q components
        Gradient qa = (da.divide(axx)).pow(2).multiply(Sa).multiply(Wa);
        Gradient qe = (de.divide(exx)).pow(2).multiply(We);
        Gradient qi = (di.divide(ixx)).pow(2).multiply(Wi);
        Gradient qw = (dw.divide(wxx)).pow(2).multiply(Wpa);
        Gradient qW = (dW.divide(Wxx)).pow(2).multiply(Wraan);

        return (one.add((Wp.multiply(P)))).multiply(qa.add(qe).add(qi).add(qw).add(qW));

    }

    private Gradient computeWaCoefficient(final FieldSpacecraftState<Gradient> fieldSpacecraftState) {
        Gradient one = fieldSpacecraftState.getA().getField().getOne();
        double   da  = FastMath.abs(fieldSpacecraftState.getA().getValue() - targetOrbitKeplerian[0]);

        Gradient Wa = one.multiply(tolerance.isFirstParameterCriterion());
        if (da < tolerance.getTolerances()[0] / 2) {

            return fieldSpacecraftState.getA().getField().getZero();
        }
        else {
            return Wa;
        }
    }

    private Gradient computeWeCoefficient(final FieldSpacecraftState<Gradient> fieldSpacecraftState) {
        Gradient one = fieldSpacecraftState.getA().getField().getOne();
        double   de  = FastMath.abs(fieldSpacecraftState.getE().getValue() - targetOrbitKeplerian[1]);

        Gradient We = one.multiply(tolerance.isSecondParameterCriterion());
        if (de < tolerance.getTolerances()[1] / 2) {

            return fieldSpacecraftState.getA().getField().getZero();
        }
        else {
            return We;
        }
    }

    private Gradient computeWiCoefficient(final FieldSpacecraftState<Gradient> fieldSpacecraftState) {
        Gradient one = fieldSpacecraftState.getA().getField().getOne();
        double di = FastMath.toDegrees(
                FastMath.abs(fieldSpacecraftState.getI().getValue() - targetOrbitKeplerian[2]) % (2 * FastMath.PI));

        Gradient Wi = one.multiply(tolerance.isThirdParameterCriterion());
        if (di < tolerance.getTolerances()[2] / 2) {

            return fieldSpacecraftState.getA().getField().getZero();
        }
        else {
            return Wi;
        }
    }

    private Gradient computeWpaCoefficient(final FieldSpacecraftState<Gradient> fieldSpacecraftState) {
        Gradient one = fieldSpacecraftState.getA().getField().getOne();
        KeplerianOrbit currentKeplerianOrbit =
                (KeplerianOrbit) OrbitType.KEPLERIAN.convertType(fieldSpacecraftState.getOrbit().toOrbit());

        double dw = FastMath.toDegrees(FastMath.abs(
                currentKeplerianOrbit.getPerigeeArgument() % (2 * FastMath.PI) - targetOrbitKeplerian[3] % (2
                        * FastMath.PI)) % (2 * FastMath.PI));

        Gradient Wpa = one.multiply(tolerance.isFourthParameterCriterion());
        if (dw < tolerance.getTolerances()[3] / 2) {

            return fieldSpacecraftState.getA().getField().getZero();
        }
        else {
            return Wpa;
        }
    }

    private Gradient computeWraanCoefficient(final FieldSpacecraftState<Gradient> fieldSpacecraftState) {
        Gradient one = fieldSpacecraftState.getA().getField().getOne();
        KeplerianOrbit currentKeplerianOrbit =
                (KeplerianOrbit) OrbitType.KEPLERIAN.convertType(fieldSpacecraftState.getOrbit().toOrbit());

        double dW = FastMath.toDegrees(FastMath.abs(
                currentKeplerianOrbit.getRightAscensionOfAscendingNode() % (2 * FastMath.PI) - targetOrbitKeplerian[4] % (2
                        * FastMath.PI)) % (2 * FastMath.PI));

        Gradient Wraan = one.multiply(tolerance.isFifthParameterCriterion());
        if (dW < tolerance.getTolerances()[4] / 2) {

            return fieldSpacecraftState.getA().getField().getZero();
        }
        else {
            return Wraan;
        }
    }

    protected boolean createConvergenceCriterion(final FieldSpacecraftState<Gradient> currentSpacecraftState) {
        double[] currentOrbitKeplerian = new double[6];
        OrbitType.KEPLERIAN.mapOrbitToArray(currentSpacecraftState.toSpacecraftState().getOrbit(), PositionAngle.TRUE,
                                            currentOrbitKeplerian, null);

        double da = FastMath.abs(currentOrbitKeplerian[0] - targetOrbitKeplerian[0]) * tolerance.isFirstParameterCriterion();
        double de =
                FastMath.abs(currentOrbitKeplerian[1] - targetOrbitKeplerian[1]) * tolerance.isSecondParameterCriterion();
        double di = FastMath.toDegrees(
                FastMath.abs(currentOrbitKeplerian[2] % FastMath.PI - targetOrbitKeplerian[2] % FastMath.PI))
                * tolerance.isThirdParameterCriterion();
        double dw = FastMath.toDegrees(
                FastMath.abs(currentOrbitKeplerian[3] % (2 * FastMath.PI) - targetOrbitKeplerian[3] % (2 * FastMath.PI)))
                * tolerance.isFourthParameterCriterion();
        double dW = FastMath.toDegrees(
                FastMath.abs(currentOrbitKeplerian[4] % (2 * FastMath.PI) - targetOrbitKeplerian[4] % (2 * FastMath.PI)))
                * tolerance.isFifthParameterCriterion();

        return (da > tolerance.getTolerances()[0] || de > tolerance.getTolerances()[1] || di > tolerance.getTolerances()[2]
                || dw > tolerance.getTolerances()[3] || dW > tolerance.getTolerances()[4]);

    }

    protected Gradient circularDifference(final Gradient w0, final double wT) {
        Gradient dv;
        Gradient one  = w0.getField().getOne();
        double   diff = w0.getValue() - wT;
        if (diff % (2. * FastMath.PI) < FastMath.PI) {
            if (diff > FastMath.PI) {
                dv = one.multiply(2. * FastMath.PI + wT).subtract(w0);
            }
            else {
                dv = w0.subtract(wT);
            }
        }

        else {
            if (diff > FastMath.PI) {
                dv = w0.subtract(2. * FastMath.PI + wT);
            }

            else {
                dv = one.multiply(wT).subtract(w0);
            }
        }

        return dv;
    }
}

















