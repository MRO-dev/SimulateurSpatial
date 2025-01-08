package org.example;

import org.hipparchus.util.FastMath;
import org.orekit.bodies.CelestialBody;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.data.DataContext;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.EquinoctialOrbit;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;

import java.io.File;
import java.util.Locale;

public class Eclipsetest {

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

        double a0     = 12.e6;
        double e0     = 0.1;
        double i0     = FastMath.toRadians(5); // Angle in radians
        double w0     = FastMath.toRadians(10);
        double W0     = FastMath.toRadians(10);
        double LV0    = 0;
        double thrust = 1; // (N)
        double mass   = 300; //(Kg)
        double ISP    = 3100; // (s)

        Orbit orbit_init =
                new KeplerianOrbit(a0, e0, i0, w0, W0, LV0, PositionAngle.TRUE, GCRF, initdate,
                                   Constants.WGS84_EARTH_MU);

// La meme en equinoctiale
        double ex0 = e0 * FastMath.cos(w0 + W0);
        double ey0 = e0 * FastMath.sin(w0 + W0);
        double hx0 = FastMath.tan(i0 / 2) * FastMath.cos(W0);
        double hy0 = FastMath.tan(i0 / 2) * FastMath.sin(W0);


        Orbit equinoctialOrbit = new EquinoctialOrbit(a0,ex0,ey0,hx0,hy0,FastMath.toRadians(20),PositionAngle.TRUE,GCRF,initdate,Constants.WGS84_EARTH_MU);

        CelestialBody sun = CelestialBodyFactory.getSun();

        OneAxisEllipsoid earth =
                new OneAxisEllipsoid(Constants.IERS2010_EARTH_EQUATORIAL_RADIUS, Constants.IERS2010_EARTH_FLATTENING, GCRF);



        // Création du propagateur
        KeplerianPropagator propagator   = new KeplerianPropagator(orbit_init);
        double              stepSize     = 60.0; // pas de temps de 60 secondes
        double              minElevation = 0.0; // élévation minimale pour déterminer une éclipse
        System.out.println(
                "orbit initiale" + OrbitType.KEPLERIAN.convertType(orbit_init));
        ;
        System.out.println(OrbitType.KEPLERIAN.convertType(orbit_init.shiftedBy(100000)));
        System.out.println("Equinoctiale orbit"+OrbitType.KEPLERIAN.convertType(equinoctialOrbit.shiftedBy(100000)) );






        // Partie Eclipse
//        EclipseDetector detector = new EclipseDetector(sun,
//                Constants.SUN_RADIUS,
//                                                        earth)
//                .withMaxCheck(60.0)
//                .withThreshold(1.0e-3);
//        Constraint eclipseConstraint = new Constraint(detector, Constraint.Type.INCLUDE);
//        List<Constraint> listConstraint = Collections.singletonList(eclipseConstraint);
//
//        // Partie avec TimeSpanWindowManeuver
//        ConstraintsProcessor       constraintsProcessor       = new ConstraintsProcessor(listConstraint);
//
//        TimeSpanMap<ManeuverWindow>
//                eclipseManeuver = constraintsProcessor.computeManeuverWindowMap(orbit_init, initdate, initdate.shiftedBy(3600));
//
//
//
//
//        // Propagation de l'orbite
//        for (int i = 0; i < 60; i++) {
//            SpacecraftState finalState = propagator.propagate(initdate.shiftedBy(60.0*i));
//
//            if (detector.g(finalState) <= 0) {
//                System.out.println("Le satellite est en éclipse à la fin de la propagation"+finalState.getDate());
//            }
//            else {
//                System.out.println("Le satellite n'est pas en éclipse à la fin de la propagation."+finalState.getDate());
//            }
//            System.out.println(eclipseManeuver.getSpan(finalState.getDate()).getData().getManeuverability()+"Resultat avec la methode de vincent");
//        }



        // Vérification si le satellite est en éclipse à la fin de la propagation

    }
}







