package org.example.leoplaw.qlaw;

import org.example.constraint.Constraint;
import org.example.leoplaw.LEOPLaw;
import org.example.leoplaw.tolerance.ToleranceBuilder;
import org.hipparchus.Field;
import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.Gradient;
import org.hipparchus.util.FastMath;
import org.orekit.orbits.FieldEquinoctialOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.FieldPropagator;
import org.orekit.propagation.FieldSpacecraftState;
import org.orekit.propagation.SpacecraftState;
import org.orekit.time.FieldAbsoluteDate;

import java.util.List;
import java.util.function.BiFunction;

public class EquinoctialQLaw extends QLaw<ToleranceBuilder.EquinoctialTolerance> implements LEOPLaw {

    // Qlaw constructor

    public EquinoctialQLaw(final Orbit initialOrbit, final Orbit targetOrbit,
                           final double maxThrust,
                           final double ISP,
                           final double initialMass, final ToleranceBuilder.EquinoctialTolerance tolerance,
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

    // Compute the QLaw from Varga paper.
    protected Gradient createQLaw(final FieldSpacecraftState<Gradient> fieldSpacecraftState) {
        // Parameters name
        //        FieldOrbit<Gradient> fieldOrbit = fieldSpacecraftState.getOrbit();
        //        FieldKeplerianOrbit<Gradient> fieldKepOrbit =
        //                (FieldKeplerianOrbit<Gradient>) OrbitType.KEPLERIAN.convertType(fieldOrbit);
        //        System.out.println(fieldKepOrbit.getPerigeeArgument().getValue());

        Gradient a   = fieldSpacecraftState.getA();
        Gradient ex  = fieldSpacecraftState.getEquinoctialEx();
        Gradient ey  = fieldSpacecraftState.getEquinoctialEy();
        Gradient hx  = fieldSpacecraftState.getHx();
        Gradient hy  = fieldSpacecraftState.getHy();
        Gradient mu  = fieldSpacecraftState.getMu();
        Gradient one = a.getField().getOne();
        Gradient p   = a.multiply(one.subtract(fieldSpacecraftState.getE().pow(2)));
        if (fieldSpacecraftState.getE().getValue() == 0) {
            p = a;
        }

        //delta
        Gradient da  = a.subtract(targetOrbitEquinoctial[0]);
        Gradient dex = ex.subtract(targetOrbitEquinoctial[1]);
        Gradient dey = ey.subtract(targetOrbitEquinoctial[2]);
        Gradient dhx = hx.subtract(targetOrbitEquinoctial[3]);
        Gradient dhy = hy.subtract(targetOrbitEquinoctial[4]);

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

        Gradient Wp  = one; // Weight of Penalty Function
        Gradient Wa  = computeWaCoefficient(fieldSpacecraftState);
        Gradient Wex = computeWexCoefficient(fieldSpacecraftState);
        Gradient Wey = computeWeyCoefficient(fieldSpacecraftState);
        Gradient Whx = computeWhxCoefficient(fieldSpacecraftState);
        Gradient Why = computeWhyCoefficient(fieldSpacecraftState);

        // Scaling function  Sa for convergence
        Gradient aMinusAT          = (a.subtract(targetOrbitEquinoctial[0])).abs();
        Gradient aMinusATDivide3aT = aMinusAT.divide(one.multiply(targetOrbitEquinoctial[0]).multiply(3));
        Gradient coeffPower4       = aMinusATDivide3aT.pow(4);
        Gradient Sa                = (one.add(coeffPower4)).sqrt();

        // useful Coefficient
        Gradient sqrtADivideMuMultiply2 = (a.divide(mu)).sqrt().multiply(2);
        Gradient sqrtPDivideMu          = (p.divide(mu)).sqrt();
        Gradient sqrtPDivideMuDivide2   = sqrtPDivideMu.multiply(0.5);
        Gradient NormExPlusEy           = ((ex.pow(2)).add((ey.pow(2)))).sqrt();
        if (ex.getValue() == 0 && ey.getValue() == 0) {
            NormExPlusEy = ex.multiply(0).add(ey.multiply(0));
        }
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

    private Gradient computeWaCoefficient(final FieldSpacecraftState<Gradient> fieldSpacecraftState) {
        Gradient one = fieldSpacecraftState.getA().getField().getOne();
        double   da  = FastMath.abs(fieldSpacecraftState.getA().getValue() - targetOrbitEquinoctial[0]);

        Gradient Wa = one.multiply(tolerance.isFirstParameterCriterion());
        if (da < tolerance.getTolerances()[0] / 2) {

            return fieldSpacecraftState.getA().getField().getZero();
        }
        else {
            return Wa;
        }
    }

    private Gradient computeWexCoefficient(final FieldSpacecraftState<Gradient> fieldSpacecraftState) {
        Gradient one = fieldSpacecraftState.getA().getField().getOne();
        double   dex = FastMath.abs(fieldSpacecraftState.getEquinoctialEx().getValue() - targetOrbitEquinoctial[1]);

        Gradient Wex = one.multiply(tolerance.isSecondParameterCriterion());
        if (dex < tolerance.getTolerances()[1] / 2) {

            return fieldSpacecraftState.getA().getField().getZero();
        }
        else {
            return Wex;
        }
    }

    private Gradient computeWeyCoefficient(final FieldSpacecraftState<Gradient> fieldSpacecraftState) {
        Gradient one = fieldSpacecraftState.getA().getField().getOne();
        double   dey = FastMath.abs(fieldSpacecraftState.getEquinoctialEy().getValue() - targetOrbitEquinoctial[2]);

        Gradient Wey = one.multiply(tolerance.isThirdParameterCriterion());
        if (dey < tolerance.getTolerances()[2] / 2) {

            return fieldSpacecraftState.getA().getField().getZero();
        }
        else {
            return Wey;
        }
    }

    private Gradient computeWhxCoefficient(final FieldSpacecraftState<Gradient> fieldSpacecraftState) {
        Gradient one = fieldSpacecraftState.getA().getField().getOne();
        double   dhx = FastMath.abs(fieldSpacecraftState.getHx().getValue() - targetOrbitEquinoctial[3]);

        Gradient Whx = one.multiply(tolerance.isFourthParameterCriterion());
        if (dhx < tolerance.getTolerances()[3] / 2) {

            return fieldSpacecraftState.getA().getField().getZero();
        }
        else {
            return Whx;
        }
    }

    private Gradient computeWhyCoefficient(final FieldSpacecraftState<Gradient> fieldSpacecraftState) {
        Gradient one = fieldSpacecraftState.getA().getField().getOne();
        double   dhy = FastMath.abs(fieldSpacecraftState.getHy().getValue() - targetOrbitEquinoctial[4]);

        Gradient Why = one.multiply(tolerance.isFifthParameterCriterion());
        if (dhy < tolerance.getTolerances()[4] / 2) {

            return fieldSpacecraftState.getA().getField().getZero();
        }
        else {
            return Why;
        }
    }

    protected boolean createConvergenceCriterion(final FieldSpacecraftState<Gradient> currentSpacecraftState) {

        double[] currentOrbitEquinoctial = new double[6];
        OrbitType.EQUINOCTIAL.mapOrbitToArray(currentSpacecraftState.toSpacecraftState().getOrbit(), PositionAngle.TRUE,
                                              currentOrbitEquinoctial, null);

        double da =
                FastMath.abs(currentOrbitEquinoctial[0] - targetOrbitEquinoctial[0]) * tolerance.isFirstParameterCriterion();
        double dex = FastMath.abs(currentOrbitEquinoctial[1] - targetOrbitEquinoctial[1])
                * tolerance.isSecondParameterCriterion();
        double dey = FastMath.abs(currentOrbitEquinoctial[2] - targetOrbitEquinoctial[2])
                * tolerance.isThirdParameterCriterion();
        double dhx = FastMath.abs(currentOrbitEquinoctial[3] - targetOrbitEquinoctial[3])
                * tolerance.isFourthParameterCriterion();
        double dhy =
                FastMath.abs(currentOrbitEquinoctial[4] - targetOrbitEquinoctial[4]) * tolerance.isFifthParameterCriterion();

        return (da > tolerance.getTolerances()[0] || dex > tolerance.getTolerances()[1]
                || dey > tolerance.getTolerances()[2] || dhx > tolerance.getTolerances()[3]
                || dhy > tolerance.getTolerances()[4]);

    }

}














