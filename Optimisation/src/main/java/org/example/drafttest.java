package org.example;

import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.EquinoctialOrbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.Constants;

public class drafttest {
    public static void main(String[] args) {
        Frame GCRF = FramesFactory.getGCRF();
        AbsoluteDate initdate = AbsoluteDate.J2000_EPOCH;
        EquinoctialOrbit newtarget = new EquinoctialOrbit(7968829.7678, 0.0110364, -2.04218497e-4, 4.3633299e-4, -1.7974e-9, 8328.004213267,
                                                          PositionAngle.TRUE, GCRF, initdate, Constants.WGS84_EARTH_MU);

        System.out.println("new Orbit Finale" + OrbitType.KEPLERIAN.convertType(newtarget));


}
}
