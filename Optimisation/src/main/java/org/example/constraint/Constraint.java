package org.example.constraint;/* Copyright 2022-2023 CS GROUP, all rights reserved. */

import org.orekit.propagation.events.EventDetector;

/**
 * Class for maneuver constraints.
 */
public class Constraint {

    /** Constraint {@link EventDetector detector}. */
    private final EventDetector detector;

    /** Constraint type. */
    private final Type type;

    public Constraint(final EventDetector detector, final Type type) {
        this.detector         = detector;
        this.type             = type;
    }

    /** @return the constraint detector */
    public EventDetector getDetector() {
        return detector;
    }

    /** @return the constraint type. */
    public Type getType() {
        return type;
    }

    /**
     * Enum indicating if this constraint shall be included or excluded.
     * <p>
     * For example, if given detector is designed to detect when spacecraft is between [-60°;60°] latitude and that the user
     * do not want to maneuver during this interval, then the constraint type shall be EXCLUDE.
     * <p>
     * On the opposite, if given detector was designed to detect when the spacecraft is <i><b>not</b></i> between [-60°;60°],
     * then the constraint type shall be INCLUDE.
     */
    public enum Type {
        /** Interval(s) defined by this constraint could be used for maneuvers. */
        INCLUDE,

        /** Interval(s) defined by this constraint shall not be used for maneuvers. */
        EXCLUDE
    }

}
