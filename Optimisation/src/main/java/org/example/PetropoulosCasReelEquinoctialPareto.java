package org.example;

import org.example.leoplaw.LEOPLawBuilder;
import org.example.leoplaw.qlaw.EquinoctialQLaw;
import org.example.leoplaw.qlaw.QLaw;
import org.example.leoplaw.tolerance.ToleranceBuilder;
import org.hipparchus.util.FastMath;
import org.orekit.data.DataContext;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.EquinoctialOrbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;

import java.io.File;
import java.util.Locale;

public class PetropoulosCasReelEquinoctialPareto {

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
        double e0     = 0.01;
        double i0     = FastMath.toRadians(10); // Angle in radians
        double w0     = FastMath.toRadians(30);
        double W0     = FastMath.toRadians(0);
        double LV0    = w0 + W0;
        double thrust = 1; // (N)
        double mass   = 300; //(Kg)
        double ISP    = 2000; // (s)

        // Target Orbit
        double aT = 12e6;
        double eT = 0.2;
        //  double iT  = FastMath.toRadians(30);
        double iT  = FastMath.toRadians(30);
        double wT  = FastMath.toRadians(0);
        double WT  = FastMath.toRadians(20);
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


        EquinoctialOrbit orbit_init =
                new EquinoctialOrbit(a0, ex0, ey0, hx0, hy0, LV0, PositionAngle.TRUE, GCRF, initdate,
                                     Constants.WGS84_EARTH_MU);

        EquinoctialOrbit orbit_target =
                new EquinoctialOrbit(aT, exT, eyT, hxT, hyT, LVT, PositionAngle.TRUE, GCRF, initdate,
                                     Constants.WGS84_EARTH_MU);
        for (int i = 0; i <21; i=i+1) {
            double doublei = i;

            double inte= 1 - doublei/20;
            LEOPLawBuilder leopLawBuilder =
                    new LEOPLawBuilder(orbit_init, orbit_target, thrust, ISP, mass).withThrustEffectivity(doublei/20)
                                                                                   .withMaximumManeuverRatioPerOrbit(inte);

            ToleranceBuilder toleranceBuilder = new ToleranceBuilder();
            ToleranceBuilder.EquinoctialTolerance equinoctialTolerance =
                    toleranceBuilder.withFirstParam(1e4).withSecondParam(0.01).withThirdParam(0.01).withFourthParam(0.01).withFifthParam(0.01)
                                    .buildEquinoctialTolerance();

            EquinoctialQLaw qlaw =
                    (EquinoctialQLaw) leopLawBuilder.buildQLawEquinoctial(equinoctialTolerance, QLaw::keplerianPropagator);

            EquinoctialQLaw.Output Qsolve = qlaw.solve();

            System.out.println("Thrust effectivity"+doublei/20+"Max Manoeuver per orbit"+inte);
            System.out.println("DVTOT" + Qsolve.getDVTotal());
            System.out.println(Qsolve.getTransferTime() + "  time transfert");



        }







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





