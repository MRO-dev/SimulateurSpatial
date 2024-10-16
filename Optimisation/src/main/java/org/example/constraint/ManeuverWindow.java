/* Copyright 2022-2023 CS GROUP, all rights reserved. */
package org.example.constraint;

public class ManeuverWindow {

    /** Enum defining if maneuvers can be performed during this maneuver window. */
    private final Maneuverability maneuverability;

    /**
     * Flag defining if this is a partially defined maneuver window .i.e that this maneuver window was not defined using two
     * distinct events.
     * <p>
     * If true, then this maneuver window may not have additional values at its initial and/or final date.
     */
    private final boolean isPartiallyDefined;

    /**
     * Constructor.
     *
     * @param maneuverability enum defining if maneuvers can be performed during this maneuver window
     * @param isPartiallyDefined flag defining if this is a partially defined maneuver window
     */
    public ManeuverWindow(final Maneuverability maneuverability,
                          final boolean isPartiallyDefined) {
        this.maneuverability               = maneuverability;
        this.isPartiallyDefined            = isPartiallyDefined;
    }

    /** @return get the enum defining if maneuvers can be performed during this maneuver window */
    public Maneuverability getManeuverability() {
        return maneuverability;
    }

    /**
     * @return get the flag defining if this is a partially defined maneuver window .i.e that this maneuver window was not
     * defined using two distinct events.
     * <p>
     * If true, then this maneuver window may not have additional values at its initial and/or final date.
     */
    public boolean isPartiallyDefined() {
        return isPartiallyDefined;
    }
}
