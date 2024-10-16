package org.example.leoplaw;

import org.example.constraint.Constraint;
import org.example.leoplaw.qlaw.EquinoctialQLaw;
import org.example.leoplaw.qlaw.KeplerianQLaw;
import org.example.leoplaw.tolerance.ToleranceBuilder;
import org.hipparchus.Field;
import org.hipparchus.analysis.differentiation.Gradient;
import org.hipparchus.util.FastMath;
import org.orekit.bodies.CelestialBody;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.EquinoctialOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.FieldPropagator;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.events.EclipseDetector;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;

public class LEOPLawBuilder {

    // Mandatory values
    private final Orbit  initialOrbit;
    private final Orbit  targetOrbit;
    private final double maxThrust;
    private final double isp;
    private final double initialMass;

    // Optional values
    private final List<Constraint> constraintList;
    private final boolean          eclipseFlag;

    private final double maximumManeuverRatioPerOrbit;
    private final double thrustEffectivity;

    public LEOPLawBuilder(final Orbit initialOrbit, final Orbit targetOrbit,
                          final double maxThrust, final double isp, final double initialMass) {

        // Threshold to avoid numerical instabilities.

        Orbit correctedInitialOrbit = computeCorrectedOrbit(initialOrbit);

        Orbit correctedTargetOrbit = computeCorrectedOrbit(targetOrbit);

        // Mandatory values
        this.initialOrbit = correctedInitialOrbit;
        this.targetOrbit  = correctedTargetOrbit;
        this.maxThrust    = maxThrust;
        this.isp          = isp;
        this.initialMass  = initialMass;

        // Default values for optional values
        this.eclipseFlag                  = false;
        this.constraintList               = new ArrayList<>();
        this.maximumManeuverRatioPerOrbit = 1.;
        this.thrustEffectivity            = 0.;
    }

    private Orbit computeCorrectedOrbit(final Orbit initialOrbit) {
        double[] initialOrbitValue = new double[6];
        OrbitType.KEPLERIAN.mapOrbitToArray(initialOrbit, PositionAngle.TRUE,
                                            initialOrbitValue, null);

        // Threshold in order to avoid numerical errors
        // Particular case we set the perigee and the raan argument to 0 in order to avoid error as they have nonsense.

        if (initialOrbitValue[3] <= 1e-5 || initialOrbit.getE() == 0) {
            initialOrbitValue[3] = 1e-4;
        }

        if (initialOrbitValue[4] <= 1e-5 || initialOrbit.getI() == 0) {
            initialOrbitValue[4] = 1e-4;
        }

        if (initialOrbit.getE() < 1e-6) {
            initialOrbitValue[1] = 1e-5;

        }
        if (initialOrbit.getI() < 1e-6) {
            initialOrbitValue[2] = 1e-5;
        }

        double ex = initialOrbitValue[1] * FastMath.cos(initialOrbitValue[3] + initialOrbitValue[4]);
        double ey = initialOrbitValue[1] * FastMath.sin(initialOrbitValue[3] + initialOrbitValue[4]);
        double hx = FastMath.tan(initialOrbitValue[2] / 2) * FastMath.cos(initialOrbitValue[4]);
        double hy = FastMath.tan(initialOrbitValue[2] / 2) * FastMath.sin(initialOrbitValue[4]);

        return new EquinoctialOrbit(initialOrbit.getA(), ex, ey, hx, hy, initialOrbit.getLv(), PositionAngle.TRUE,
                                    initialOrbit.getFrame(), initialOrbit.getDate(), initialOrbit.getMu());

    }

    public LEOPLawBuilder(final Orbit initialOrbit, final Orbit targetOrbit,
                          final double maxThrust, final double isp,
                          final double initialMass,
                          final double thrustEffectivity,
                          final double maximumManeuverRatioPerOrbit, final List<Constraint> constraintList,
                          final boolean eclipseFlag
                         ) {

        // Mandatory values
        this.initialOrbit = initialOrbit;
        this.targetOrbit  = targetOrbit;
        this.maxThrust    = maxThrust;
        this.isp          = isp;
        this.initialMass  = initialMass;

        // Default values for optional values
        this.eclipseFlag                  = eclipseFlag;
        this.maximumManeuverRatioPerOrbit = maximumManeuverRatioPerOrbit;
        this.thrustEffectivity            = thrustEffectivity;

        if (eclipseFlag) {
            this.constraintList = addEclipseConstraint(constraintList);
        }
        else {
            this.constraintList = constraintList;
        }
    }

    public LEOPLaw buildQLawEquinoctial(final ToleranceBuilder.EquinoctialTolerance tolerance,
                                        final BiFunction<Field<Gradient>, SpacecraftState, FieldPropagator<Gradient>>
                                                fieldPropagatorMapper) {

        return new EquinoctialQLaw(initialOrbit, targetOrbit, maxThrust, isp, initialMass, tolerance, thrustEffectivity,
                                   maximumManeuverRatioPerOrbit,
                                   constraintList,
                                   fieldPropagatorMapper);
    }

    public LEOPLaw buildQLawKeplerian(final ToleranceBuilder.KeplerianTolerance tolerance,
                                      final BiFunction<Field<Gradient>, SpacecraftState, FieldPropagator<Gradient>>
                                              fieldPropagatorMapper) {

        return new KeplerianQLaw(initialOrbit, targetOrbit, maxThrust, isp, initialMass, tolerance, thrustEffectivity,
                                 maximumManeuverRatioPerOrbit,
                                 constraintList,
                                 fieldPropagatorMapper);
    }

    private List<Constraint> addEclipseConstraint(final List<Constraint> constraintList) {
        final Frame         itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, true);
        final CelestialBody sun  = CelestialBodyFactory.getSun();

        final OneAxisEllipsoid earth =
                new OneAxisEllipsoid(Constants.IERS2010_EARTH_EQUATORIAL_RADIUS, Constants.IERS2010_EARTH_FLATTENING, itrf);

        EclipseDetector detector = new EclipseDetector(sun, Constants.SUN_RADIUS, earth)
                .withMaxCheck(60.0);
        Constraint eclipseConstraint = new Constraint(detector, Constraint.Type.INCLUDE);
        constraintList.add(eclipseConstraint);

        return constraintList;

    }

    public LEOPLawBuilder withMaximumManeuverRatioPerOrbit(final double maximumManeuverRatioPerOrbit) {
        // Ajouter une erreur comme quoi maximum ratio <1 et >0
        return new LEOPLawBuilder(initialOrbit, targetOrbit, maxThrust, isp,
                                  initialMass, thrustEffectivity, maximumManeuverRatioPerOrbit, constraintList,
                                  eclipseFlag);

    }

    public LEOPLawBuilder withThrustEffectivity(final double thrustEffectivity) {
        // Ajouter erreur comme quoi valeur >0 et <1
        return new LEOPLawBuilder(initialOrbit, targetOrbit, maxThrust, isp,
                                  initialMass, thrustEffectivity, maximumManeuverRatioPerOrbit, constraintList,
                                  eclipseFlag);
    }

    public LEOPLawBuilder withEclipse() {
        return new LEOPLawBuilder(initialOrbit, targetOrbit, maxThrust, isp,
                                  initialMass, thrustEffectivity, maximumManeuverRatioPerOrbit, constraintList,
                                  true);
    }

    public LEOPLawBuilder withConstraints(final List<Constraint> constraintList) {
        return new LEOPLawBuilder(initialOrbit, targetOrbit, maxThrust, isp,
                                  initialMass, thrustEffectivity, maximumManeuverRatioPerOrbit, constraintList,
                                  eclipseFlag);
    }

}
