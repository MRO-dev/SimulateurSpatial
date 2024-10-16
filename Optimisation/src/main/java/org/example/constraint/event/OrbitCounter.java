/* Copyright 2022-2023 CS GROUP, all rights reserved. */
package org.example.constraint.event;

import org.hipparchus.ode.events.Action;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.events.EventDetector;
import org.orekit.propagation.events.EventsLogger;
import org.orekit.propagation.events.handlers.EventHandler;

import java.util.HashMap;

//TODO TO FINISH
/**
 * Event processor allowing for orbit counting.
 * <p>
 * Gives orbit number associated with a specific event.
 *
 * @author CS Group
 */
public class OrbitCounter implements EventProcessor, EventHandler {

    /** Orbit number. */
    private int orbitNumber;

    /** Map of orbit number with respect to a specific state. */
    private final HashMap<SpacecraftState, Integer> orbitNumberMap;

    /** Default constructor starting with orbit number 0. */
    public OrbitCounter() {
        this(0);
    }

    /**
     * Constructor.
     *
     * @param orbitNumber initial orbit number to start from
     */
    public OrbitCounter(final int orbitNumber) {
        this.orbitNumber    = orbitNumber;
        this.orbitNumberMap = new HashMap<>();
    }

    /** {@inheritDoc} */
    @Override
    public double process(final EventsLogger.LoggedEvent event) {
        return orbitNumberMap.get(event.getState());
    }

    /** {@inheritDoc} */
    @Override
    public Action eventOccurred(final SpacecraftState spacecraftState, final EventDetector eventDetector, final boolean b) {
        orbitNumber++;
        orbitNumberMap.put(spacecraftState, orbitNumber);
        return Action.CONTINUE;
    }
}
