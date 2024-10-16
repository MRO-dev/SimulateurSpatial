package org.example;

import org.example.constraint.Constraint;
import org.example.leoplaw.LEOPLawBuilder;
import org.example.leoplaw.qlaw.EquinoctialQLaw;
import org.example.leoplaw.qlaw.QLaw;
import org.example.leoplaw.tolerance.ToleranceBuilder;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.orekit.data.DataContext;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.forces.maneuvers.ConstantThrustManeuver;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.EquinoctialOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.events.ParameterDrivenDateIntervalDetector;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public class equinoctialQlawTest {

    public static void main(String[] args) {
        try {

            // configure Orekit
            final File home       = new File(System.getProperty("user.home"));
            final File orekitData = new File(home, "orekit-data");
            if (!orekitData.exists()) {
                System.err.format(Locale.US, "Failed to find %s folder%n",
                                  orekitData.getAbsolutePath());
                System.err.format(Locale.US,
                                  "You need to download %s from %s, unzip it in %s and rename it 'orekit-data' for this tutorial to work%n",
                                  "orekit-data-master.zip",
                                  "https://gitlab.orekit.org/orekit/orekit-data/-/archive/master/orekit-data-master.zip",
                                  home.getAbsolutePath());
                System.exit(1);
            }
            final DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
            manager.addProvider(new DirectoryCrawler(orekitData));

            System.out.println(TimeScalesFactory.getUTC().getName());
        }
        catch (OrekitException oe) {
            System.err.println(oe.getLocalizedMessage());
        }
        Frame        GCRF     = FramesFactory.getGCRF();
        AbsoluteDate initdate = AbsoluteDate.J2000_EPOCH;

        // Initit orbit
        double a0     = 7e6;
        double e0     = 0.03;
        double i0     = FastMath.toRadians(5); // Angle in radians
        double w0     = FastMath.toRadians(10);
        double W0     = FastMath.toRadians(10);
        double LV0    = w0 + W0;
        double thrust = 1; // (N)
        double mass   = 300; //(Kg)
        double ISP    = 3100; // (s)

        // Target Orbit
        double aT = 12e6;
        double eT = 0.1;
        //  double iT  = FastMath.toRadians(30);
        double iT  = FastMath.toRadians(5);
        double wT  = FastMath.toRadians(10);
        double WT  = FastMath.toRadians(10);
        double LVT = wT + WT;

        // Transform into equinoctial parameters Next step automatize it
        double ex0 = e0 * FastMath.cos(w0 + W0);
        double exT = eT * FastMath.cos(wT + WT);

        double ey0 = e0 * FastMath.sin(w0 + W0);
        double eyT = eT * FastMath.sin(wT + WT);

        double hx0 = FastMath.tan(i0 / 2) * FastMath.cos(W0);
        double hxT = FastMath.tan(iT / 2) * FastMath.cos(WT);

        double hy0 = FastMath.tan(i0 / 2) * FastMath.sin(W0);
        double hyT = FastMath.tan(iT / 2) * FastMath.sin(WT);

//        Orbit orbit_init =
//                new KeplerianOrbit(a0, e0, i0, w0, W0, LV0, PositionAngle.TRUE, GCRF, initdate,
//                                   Constants.WGS84_EARTH_MU);
//
//        Orbit orbit_target =
//                new KeplerianOrbit(aT, eT, iT, wT, WT, LVT, PositionAngle.TRUE, GCRF, initdate,
//                                   Constants.WGS84_EARTH_MU);

                EquinoctialOrbit orbit_init =
                        new EquinoctialOrbit(a0, ex0, ey0, hx0, hy0, LV0, PositionAngle.TRUE, GCRF, initdate,
                                             Constants.WGS84_EARTH_MU);

                EquinoctialOrbit orbit_target =
                        new EquinoctialOrbit(aT, exT, eyT, hxT, hyT, LVT, PositionAngle.TRUE, GCRF, initdate,
                                             Constants.WGS84_EARTH_MU);

        ParameterDrivenDateIntervalDetector detector =
                new ParameterDrivenDateIntervalDetector("1er intervalle", initdate.shiftedBy(60.0),
                                                        initdate.shiftedBy(3600.0));

        ParameterDrivenDateIntervalDetector detector1 =
                new ParameterDrivenDateIntervalDetector("2eme intervalle", initdate.shiftedBy(15000.0),
                                                        initdate.shiftedBy(100000.0));

        Constraint constraintDate  = new Constraint(detector, Constraint.Type.EXCLUDE);
        Constraint constraintDate1 = new Constraint(detector1, Constraint.Type.EXCLUDE);

        List<Constraint> constraintList = new ArrayList<>();

        constraintList.add(constraintDate);
        constraintList.add(constraintDate1);

        LEOPLawBuilder leopLawBuilder =
                new LEOPLawBuilder(orbit_init, orbit_target, thrust, ISP, mass);

        ToleranceBuilder toleranceBuilder = new ToleranceBuilder();
        ToleranceBuilder.EquinoctialTolerance equinoctialTolerance =
                toleranceBuilder.withFirstParam(2e4).withSecondParam(0.05).withThirdParam(0.05).buildEquinoctialTolerance();
        for (int i = 0; i < 5; i++) {
            System.out.println(equinoctialTolerance.getTolerances()[i]);

        }
        System.out.println(equinoctialTolerance.isFirstParameterCriterion());
        System.out.println(equinoctialTolerance.isSecondParameterCriterion());
        System.out.println(equinoctialTolerance.isThirdParameterCriterion());
        System.out.println(equinoctialTolerance.isFourthParameterCriterion());
        System.out.println(equinoctialTolerance.isFifthParameterCriterion());


        EquinoctialQLaw qlaw =
                (EquinoctialQLaw) leopLawBuilder.buildQLawEquinoctial(equinoctialTolerance, QLaw::keplerianPropagator);

        EquinoctialQLaw.Output Qsolve = qlaw.solve();


        System.out.println(Qsolve.getManeuvers().size());

        System.out.println(Qsolve.getEndDate());
        System.out.println("DVTOT" + Qsolve.getDVTotal());
        System.out.println(Qsolve.getDeltaM() + "Delta M");
        System.out.println(Qsolve.getTransferTime() + "  time transfert");

        System.out.println(
                "orbit target propagée à t final" + OrbitType.EQUINOCTIAL.convertType(Qsolve.getTargetOrbitFinalState()));
        System.out.println(
                "orbit target propagée à t final" + OrbitType.KEPLERIAN.convertType(Qsolve.getTargetOrbitFinalState()));

        System.out.println("Orbit propagée finale" + OrbitType.EQUINOCTIAL.convertType(Qsolve.getOrbitFinalState()));

        System.out.println("Orbit propagée finale" + (Qsolve.getOrbitFinalState()));
        System.out.println(
                "coordonnée de l'orbite propagée à t final" + Qsolve.getOrbitFinalState().getPVCoordinates().getPosition());

        // Partie phasage
        List<ConstantThrustManeuver> newManeuvers = Qsolve.getPhasingManeuvers();
        double dt =
                newManeuvers.get(newManeuvers.size() - 1).getEndDate().durationFrom(initdate);
        System.out.println(" new transfert time is : " + dt);
        System.out.println("orbit target a nouveau t phasing" + OrbitType.KEPLERIAN.convertType(orbit_target.shiftedBy(dt)));
        System.out.println(
                "coordonnée de l'orbite à t final phasing" + OrbitType.EQUINOCTIAL.convertType(orbit_target.shiftedBy(dt))
                                                                                .getPVCoordinates().getPosition());

        // Comparaison Hohman
        HohmanTransfert hohmanTransfert =
                new HohmanTransfert(a0, aT, orbit_init.getMu(), mass, ISP);
        double deltaInclination = FastMath.abs(i0 - iT);
        System.out.println(
                "Hohman transfert deltaV \n" + hohmanTransfert.getDeltaV(deltaInclination)
                        + "\n get deltaM for Hohman transfert\n"
                        + hohmanTransfert.getDeltaM(deltaInclination));

        System.out.println(hohmanTransfert.getMinimumDeltaV());

        List<Orbit> orbitListTarget = new ArrayList<>();
        orbitListTarget.add(orbit_target);
        double period   = orbit_target.getKeplerianPeriod();
        int    timestep = FastMath.toIntExact(FastMath.round(period / 260));
        int    index    = 0;
        for (int t = 0; t < period; t = t + timestep) {
            orbitListTarget.add(orbitListTarget.get(index).shiftedBy(timestep));
            index += 1;
        }

        //Create output

        final StringBuilder dataManeuver = new StringBuilder();
        for (ConstantThrustManeuver maneuver : Qsolve.getManeuvers()) {
            addManeuverToOutput(dataManeuver, initdate, maneuver);
        }

        // Write output maneuvers
        final File rootFileManeuvers   = new File(System.getProperty("user.home"));
        final File outputFileManeuvers = new File(rootFileManeuvers, "Benchmark/maneuversTestA.csv");
        BenchmarkUtils.writeData(outputFileManeuvers, dataManeuver);

        // Write output position
        final StringBuilder dataPosition = new StringBuilder();
        for (SpacecraftState spacecraftState : Qsolve.getSpacecraftEvolution()) {
            addPositionSpacecraftStateToOutput(dataPosition, spacecraftState, initdate);
        }
        // Write output position
        final File rootFilePos   = new File(System.getProperty("user.home"));
        final File outputFilePos = new File(rootFilePos, "Benchmark/PositionTestA.csv");
        BenchmarkUtils.writeData(outputFilePos, dataPosition);

        //        // Write output position TargetSpacecraft
        final StringBuilder dataTarget = new StringBuilder();
        for (Orbit orbit : orbitListTarget) {
            addPositionTargetToOutput(dataTarget, orbit);
        }
        // Write output
        final File rootFileTarget   = new File(System.getProperty("user.home"));
        final File outputFileTarget = new File(rootFileTarget, "Benchmark/PositionTargetTestA.csv");
        BenchmarkUtils.writeData(outputFileTarget, dataTarget);

        //        // Write output L2Error evolution
        final StringBuilder dataL2error = new StringBuilder();

        addL2EvolutiontoOutput(dataL2error, Qsolve);

        // Write output
        final File rootFileL2error   = new File(System.getProperty("user.home"));
        final File outputFileL2error = new File(rootFileL2error, "Benchmark/L2errorTestA.csv");
        BenchmarkUtils.writeData(outputFileL2error, dataL2error);

    }

    private static void addPositionTargetToOutput(final StringBuilder dataTarget, final Orbit orbit) {

        final Vector3D position = orbit.getPVCoordinates().getPosition();
        BenchmarkUtils.addData(dataTarget, position.toArray());
    }

    private static void addL2EvolutiontoOutput(final StringBuilder dataL2Error,
                                               final EquinoctialQLaw.Output Qsolve) {

        BenchmarkUtils.addData(dataL2Error,
                               Qsolve.getL2ErrorEvolution().stream().mapToDouble(Double::doubleValue).toArray());
    }

    private static void addManeuverToOutput(final StringBuilder data, final AbsoluteDate initialDate,
                                            final ConstantThrustManeuver maneuver) {

        final double   timeFromStart = maneuver.getStartDate().durationFrom(initialDate);
        final Vector3D thrust        = maneuver.getThrustVector();
        BenchmarkUtils.addData(data, timeFromStart, thrust.toArray());
    }

    private static void addPositionSpacecraftStateToOutput(final StringBuilder data,
                                                           final SpacecraftState spacecraftState,
                                                           final AbsoluteDate initdate) {

        final Vector3D position = spacecraftState.getPVCoordinates().getPosition();
        BenchmarkUtils.addData(data, spacecraftState.getDate().durationFrom(initdate), position.toArray());

    }

    //                SpacecraftState state = new SpacecraftState(orbit_init, mass);
    //                final double    dP    = 1;
    //                final double               minStep    = 0.001;
    //                final double               maxStep    = 500;
    //                final double               initStep   = 60;
    //                final double[][]           tolerance  = NumericalPropagator.tolerances(dP, orbit_init, OrbitType.EQUINOCTIAL);
    //                AdaptiveStepsizeIntegrator integrator = new DormandPrince853Integrator(minStep, maxStep, tolerance[0], tolerance[1]);
    //                integrator.setInitialStepSize(initStep);
    //                NumericalPropagator propagator = new NumericalPropagator(integrator);
    //                propagator.setInitialState(state);
    //
    //        // Adding sep handler for progression monitoring
    //        final double propagationDuration = Qsolve.getEndDate().durationFrom(initdate);
    //        propagator.getMultiplexer().add(3600, (s) -> {
    //            final double progress = s.getDate().durationFrom(newDate) / propagationDuration;
    //            System.out.format("Progress : %4.2f | sma (km) = %f \n", progress * 100, s.getA() / 1000);
    //        });
    //
    //        // Adding maneuvers
    //        for (int i = 0; i < Qsolve.maneuvers.size() - 1; i++) {
    //
    //            propagator.addForceModel(Qsolve.maneuvers.get(i));
    //        }
    //
    //        System.out.println("Starting propagation");
    //        state = propagator.propagate(newDate, Qsolve.getEndDate());
    //
    //        System.out.println(state.getOrbit() + "Propagation de l'orbite sans le phasing");
    //        System.out.println("new Orbit Finale" + OrbitType.KEPLERIAN.convertType(state.getOrbit()));

    // Same code for phase
    //                // Adding sep handler for progression monitoring
    //                final double propagationDuration = newManeuvers.get(newManeuvers.size()-1).getEndDate().durationFrom(initdate);
    //                propagator.getMultiplexer().add(3600, (s) -> {
    //                    final double progress = s.getDate().durationFrom(initdate) / propagationDuration;
    //                    System.out.format("Progress : %4.2f | sma (km) = %f \n", progress * 100, s.getA() / 1000);
    //                });
    //
    //                // Adding maneuvers
    //                for (int i = 0; i < newManeuvers.size() - 1; i++) {
    //
    //                    propagator.addForceModel(newManeuvers.get(i));
    //                }
    //
    //                System.out.println("Starting propagation");
    //                state = propagator.propagate(initdate, newManeuvers.get(newManeuvers.size()-1).getEndDate());
    //
    //                System.out.println(state.getOrbit() + "Propagation de l'orbite avec les manoeuvres calculées");
    //        System.out.println("new Orbit Finale" + OrbitType.KEPLERIAN.convertType(state.getOrbit()));

}



