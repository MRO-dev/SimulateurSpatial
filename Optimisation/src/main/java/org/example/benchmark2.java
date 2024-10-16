package org.example;

import org.hipparchus.util.FastMath;
import org.orekit.data.DataContext;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.EquinoctialOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;

import java.io.File;
import java.util.Locale;

public class benchmark2 {

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
        double a0     = 7.e6;
        double e0     = 0.01;
        double i0     = FastMath.toRadians(0.05); // Angle in radians
        double w0     = FastMath.toRadians(0);
        double W0     = 0;
        double LV0    = 0;
        double thrust = 10; // (N)
        double mass   = 300; //(Kg)
        double ISP    = 1080; // (s) non utilisé pour le moment
        ;

        // Target Orbit
        double aT  = 42.e6;
        double eT  = 0.001;
        double iT  = FastMath.toRadians(60);
        double wT  = 0;
        double WT  = FastMath.toRadians(0);
        double LVT = 0;

        // Transform into equinoctial parameters Next step automatize it
        double ex0 = e0 * FastMath.cos(w0 + W0);
        double exT = eT * FastMath.cos(wT + WT);

        double ey0 = e0 * FastMath.sin(w0 + W0);
        double eyT = eT * FastMath.sin(wT + WT);

        double hx0 = FastMath.tan(i0 / 2) * FastMath.cos(W0);
        double hxT = FastMath.tan(iT / 2) * FastMath.cos(WT);

        double hy0 = FastMath.tan(i0 / 2) * FastMath.sin(W0);
        double hyT = FastMath.tan(iT / 2) * FastMath.sin(WT);
        Orbit orbit_init =
                new EquinoctialOrbit(a0, ex0, ey0, hx0, hy0, LV0, PositionAngle.TRUE, GCRF, initdate,
                                     Constants.WGS84_EARTH_MU);
        Orbit orbit_target =
                new EquinoctialOrbit(aT, exT, eyT, hxT, hyT, LVT, PositionAngle.TRUE, GCRF, initdate,
                                     Constants.WGS84_EARTH_MU);

//        List<Vector3D> solution = new ArrayList<>();
//        Vector3D nul = new Vector3D(1,0,0);
//        solution.add(nul);
//
//        List<String[]> vector4D =new ArrayList<String[]>();
//
//        String[] testliste;
//        testliste = new String[] { initdate.toString(), nul.toString()};
//        String[] testliste2;
//        testliste2 = new String[] {initdate.shiftedBy(100).toString(), nul.scalarMultiply(3).toString()};
//        vector4D.add(testliste);
//        System.out.println(vector4D.get(0));
        QLawOrbit qlaw = new QLawOrbit(orbit_init, orbit_target, thrust, ISP);


    }

}
