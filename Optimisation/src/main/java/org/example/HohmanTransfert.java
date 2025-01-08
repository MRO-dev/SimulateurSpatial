package org.example;

import org.hipparchus.util.FastMath;
import org.orekit.utils.Constants;

public class HohmanTransfert {

    private final double ISP;
    private final double g0 = Constants.G0_STANDARD_GRAVITY;

    private final double initialRadius;
    private final double targetRadius;
    private final double mu;

    private final double initialMass;

    public HohmanTransfert(final double initialRadius, final double targetRadius) {
        this(initialRadius, targetRadius, Constants.WGS84_EARTH_MU, 1000, 2500);

    }

    public HohmanTransfert(final double initialRadius, final double targetRadius, final double mu) {
        this(initialRadius, targetRadius, mu, 1000, 2500);

    }

    public HohmanTransfert(final double initialRadius, final double targetRadius, final double mu, final double initialMass,
                           final double ISP) {
        this.initialRadius = initialRadius;
        this.targetRadius  = targetRadius;
        this.initialMass   = initialMass;
        this.ISP           = ISP;
        this.mu            = mu;

    }

    public double getDeltaV() {
        return getDeltaV(0);
    }

    public double getDeltaV(final double deltaInclination) {
        // Get DeltaV1
        double deltaV1 =
                FastMath.sqrt((2 * mu * targetRadius / (initialRadius * (initialRadius + targetRadius)))) - FastMath.sqrt(
                        mu / initialRadius);
        // Get DeltaV2
        double deltaV2 = FastMath.sqrt(mu / targetRadius) - FastMath.sqrt(
                (2 * mu * initialRadius) / (targetRadius * (initialRadius + targetRadius)));

        //Get DeltaV3 due to inclination between the 2 orbits.
        // We choose to change the inclination at the orbit where the radius is maximum to minimize deltaV
        double deltaV3 = 0;
        if (deltaInclination != 0) {
            deltaV3 = 2 * FastMath.sqrt(mu / FastMath.max(targetRadius, initialRadius)) * FastMath.sin(deltaInclination / 2);
        }

        return deltaV1 + deltaV2 + deltaV3;
    }

    public double getInclinationDeltaV(final double deltaInclination) {
        return 2 * FastMath.sqrt(mu / FastMath.max(targetRadius, initialRadius)) * FastMath.sin(deltaInclination / 2);

    }

    public double getDeltaM() {
        return getDeltaM(0);
    }

    public double getDeltaM(final double deltaInclination) {
        double DV = getDeltaV(deltaInclination);
        return initialMass * (1 - FastMath.exp((-DV) / (g0 * ISP)));

    }

    public double getMinimumDeltaV() {
        double V1 = FastMath.sqrt(Constants.WGS84_EARTH_MU*(FastMath.abs(-initialRadius+2*targetRadius))/(initialRadius*targetRadius));
        double deltaV =V1-FastMath.sqrt(Constants.WGS84_EARTH_MU/initialRadius);
        return
               deltaV;
    }

}
