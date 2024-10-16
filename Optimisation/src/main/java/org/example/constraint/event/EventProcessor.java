/* Copyright 2022-2023 CS GROUP, all rights reserved. */
package org.example.constraint.event;

import org.orekit.propagation.events.EventsLogger;

/**
 * Interface for event processor classes.
 * <p>
 * Designed to return a value associated to given {@link org.orekit.propagation.events.EventsLogger.LoggedEvent}.
 *
 * @author CS Group
 */
public interface EventProcessor {

    /**
     * Process given event and compute a value accordingly.
     *
     * @param event event to process
     *
     * @return value associated to processed event
     */
    double process(EventsLogger.LoggedEvent event);
}
