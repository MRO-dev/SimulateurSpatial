/* Copyright 2022-2023 CS GROUP, all rights reserved. */
package org.example.constraint;

import org.example.constraint.event.EventProcessor;
import org.orekit.orbits.Orbit;
import org.orekit.propagation.Propagator;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.orekit.propagation.conversion.PropagatorBuilder;
import org.orekit.propagation.events.BooleanDetector;
import org.orekit.propagation.events.EventDetector;
import org.orekit.propagation.events.EventsLogger;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.TimeSpanMap;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Process given maneuver constraints to compute available maneuver windows.
 *
 * @author CS Group
 */
public class ConstraintsProcessor {

    /** Constraints to take into accounts when computing available maneuver window(s). */
    private final List<Constraint> constraints;

    /**
     * Constructor.
     *
     * @param constraints constraints that will define maneuver windows
     */
    public ConstraintsProcessor(final List<Constraint> constraints) {
        this.constraints = constraints;
    }

    /**
     * Add new constraints to respect.
     *
     * @param additionalConstraints additional constraints to respect
     */
    public void addConstraints(final List<Constraint> additionalConstraints) {
        constraints.addAll(additionalConstraints);
    }

    /**
     * Add new constraint to respect.
     *
     * @param additionalConstraint additional constraint to respect
     */
    public void addConstraint(final Constraint additionalConstraint) {
        constraints.add(additionalConstraint);
    }

    /**
     * Compute a {@link TimeSpanMap maneuver window map} specifying when maneuvers can be done while respecting every given
     * constraint.
     *
     * @param initialPropagationDate initial propagation date
     * @param finalPropagationDate final propagation date
     *
     * @return maneuver window map
     */
    public TimeSpanMap<ManeuverWindow> computeManeuverWindowMap(final Orbit orbit,
                                                                final AbsoluteDate initialPropagationDate,
                                                                final AbsoluteDate finalPropagationDate) {

        // No constraints, the propagation interval is considered as a maneuver window
        if (constraints.isEmpty()) {
            return getDefaultManeuverWindow(initialPropagationDate, finalPropagationDate);
        }
        // Create boolean detector taking into accounts every constraint
        final EventDetector maneuverWindowDetector = createManeuverWindowDetector();

        // Monitor maneuver window detector
        final EventsLogger  eventsLogger                    = new EventsLogger();
        final EventDetector monitoredManeuverWindowDetector = eventsLogger.monitorDetector(maneuverWindowDetector);

        // Create propagator

        final KeplerianPropagator propagator = new KeplerianPropagator(orbit);

        // Add monitored maneuver window detector
        propagator.addEventDetector(monitoredManeuverWindowDetector);

        // Propagate
        propagator.propagate(initialPropagationDate, finalPropagationDate);

        // Store logged events
        final List<EventsLogger.LoggedEvent> loggedEvents = eventsLogger.getLoggedEvents();

        // Create maneuver window map
        return createManeuverWindowMap(initialPropagationDate, finalPropagationDate, loggedEvents);

    }

    /**
     * Get default maneuver window map where propagation interval is a maneuver window.
     *
     * @param initialPropagationDate initial propagation date
     * @param finalPropagationDate final propagation date
     *
     * @return default maneuver window map
     */
    private TimeSpanMap<ManeuverWindow> getDefaultManeuverWindow(final AbsoluteDate initialPropagationDate,
                                                                 final AbsoluteDate finalPropagationDate) {
        final TimeSpanMap<ManeuverWindow> maneuverWindowMap =
                new TimeSpanMap<>(getDefaultEmptyUnavailableManeuverWindow());

        maneuverWindowMap.addValidBetween(getDefaultEmptyAvailableManeuverWindow(),
                                          initialPropagationDate, finalPropagationDate);

        return maneuverWindowMap;
    }

    /**
     * Get default empty available maneuver window entry.
     *
     * @return default empty available maneuver window entry
     */
    private ManeuverWindow getDefaultEmptyAvailableManeuverWindow() {
        final boolean initialWindowIsPartiallyDefined = false;
        return new ManeuverWindow(Maneuverability.CAN, initialWindowIsPartiallyDefined);
    }

    /**
     * Get default empty unavailable maneuver window entry.
     *
     * @return default empty unavailable maneuver window entry
     */
    private ManeuverWindow getDefaultEmptyUnavailableManeuverWindow() {
        final boolean initialWindowIsPartiallyDefined = true;
        return new ManeuverWindow(Maneuverability.CANNOT, initialWindowIsPartiallyDefined);
    }

    /**
     * Create a maneuver window detector.
     *
     * @return maneuver window detector
     */
    private EventDetector createManeuverWindowDetector() {

        // There are constraints to include and exclude
        if (thereIsAtLeastOneConstraintOfType(Constraint.Type.INCLUDE) &&
                thereIsAtLeastOneConstraintOfType(Constraint.Type.EXCLUDE)) {

            // Combine all constraints that shall be included
            final EventDetector includeDetector =
                    combineConstraintsWrtType(Constraint.Type.INCLUDE);

            // Combine all constraints that shall be excluded
            final EventDetector excludeDetector =
                    BooleanDetector.notCombine(combineConstraintsWrtType(Constraint.Type.EXCLUDE));

            return BooleanDetector.andCombine(includeDetector, excludeDetector).withThreshold(1e-6);
        }

        // There are only constraints to include
        else if (thereIsAtLeastOneConstraintOfType(Constraint.Type.INCLUDE)) {
            // Combine all constraints that shall be included
            return combineConstraintsWrtType(Constraint.Type.INCLUDE);
        }
        // There are only constraints to exclude
        else {
            // Combine all constraints that shall be excluded
            return BooleanDetector.notCombine(combineConstraintsWrtType(Constraint.Type.EXCLUDE));
        }
    }

    /**
     * Combine maneuver constraints depending on their type into a single {@link EventDetector event detector}.
     *
     * @param type maneuver constraint type to combine
     *
     * @return detector combining constraints of given type
     */
    private EventDetector combineConstraintsWrtType(final Constraint.Type type) {

        // Filter constraints according to their type
        final Stream<Constraint> filteredConstraintsStream =
                constraints.stream().filter(constraint -> constraint.getType() == type);

        // Map event detectors
        final Stream<EventDetector> filteredDetectorsStream = filteredConstraintsStream.map(Constraint::getDetector);

        // Convert to set
        final Set<EventDetector> filteredDetectorsSet = filteredDetectorsStream.collect(Collectors.toSet());

        // Combine all filtered detectors according to their type
        if (type == Constraint.Type.EXCLUDE) {
            return BooleanDetector.orCombine(filteredDetectorsSet);
        }
        else {
            return BooleanDetector.andCombine(filteredDetectorsSet);
        }

    }

    private boolean thereIsAtLeastOneConstraintOfType(final Constraint.Type type) {
        final Stream<Constraint> constraintsStream = constraints.stream();

        final Stream<Constraint> filteredConstraintsStream =
                constraintsStream.filter(constraint -> constraint.getType() == type);

        return filteredConstraintsStream.findAny().isPresent();
    }

    /**
     * Create a propagator using given propagator builder.
     *
     * @param propagatorBuilder propagator builder
     *
     * @return propagator
     */
    private Propagator createPropagator(final PropagatorBuilder propagatorBuilder) {
        return propagatorBuilder.buildPropagator(propagatorBuilder.getSelectedNormalizedParameters());
    }

    /**
     * Create a {@link TimeSpanMap maneuver window map} specifying when maneuvers are allowed.
     *
     * @param initialPropagationDate initial propagation date
     * @param finalPropagationDate final propagation date
     * @param loggedEvents Events logged from the constraints detector
     *
     * @return maneuver window map specifying when maneuvers are allowedQ
     */
    private TimeSpanMap<ManeuverWindow> createManeuverWindowMap(final AbsoluteDate initialPropagationDate,
                                                                final AbsoluteDate finalPropagationDate,
                                                                final List<EventsLogger.LoggedEvent> loggedEvents) {

        // Handle specific case where there are only exclude constraints and no event was detected, propagation interval is
        // then considered as a maneuver window
        if (thereAreOnlyExcludeConstraints() && loggedEvents.isEmpty()) {
            return getDefaultManeuverWindow(initialPropagationDate, finalPropagationDate);
        }

        // Initialize maneuver window map
        final boolean initialWindowIsPartiallyDefined = true;
        final TimeSpanMap<ManeuverWindow> maneuverWindowMap =
                new TimeSpanMap<>(new ManeuverWindow(Maneuverability.CANNOT, initialWindowIsPartiallyDefined));

        // Add maneuver window
        for (int i = 0; i < loggedEvents.size(); i++) {
            final AbsoluteDate start;
            final AbsoluteDate end;
            final boolean      isPartiallyDefined;

            final EventsLogger.LoggedEvent currentLoggedEvent = loggedEvents.get(i);

            // Entering a maneuver window
            if (currentLoggedEvent.isIncreasing()) {
                if (i < loggedEvents.size() - 1) {
                    start = currentLoggedEvent.getDate();
                    end   = loggedEvents.get(i + 1).getDate();

                    isPartiallyDefined = false;
                    // Incrementing i to go to the next maneuver window interval
                    i++;
                }
                // Handling specific case at the end of the propagation
                else {
                    start              = loggedEvents.get(i).getDate();
                    end                = finalPropagationDate;
                    isPartiallyDefined = true;
                }

            }
            // Exiting a maneuver window, handling specific case at the beginning of the propagation
            else {
                start              = initialPropagationDate;
                end                = loggedEvents.get(i).getDate();
                isPartiallyDefined = true;
            }

            maneuverWindowMap.addValidBetween(new ManeuverWindow(Maneuverability.CAN, isPartiallyDefined),
                                              start, end);
        }

        return maneuverWindowMap;
    }

    private boolean thereAreOnlyExcludeConstraints() {
        return thereIsAtLeastOneConstraintOfType(Constraint.Type.EXCLUDE) &&
                !thereIsAtLeastOneConstraintOfType(Constraint.Type.INCLUDE);
    }

    /**
     * Loop through given event processors to build an array of additional values.
     *
     * @param event event to process
     * @param eventProcessors event processors
     *
     * @return array of additional values
     */
    private double[] loopThroughEventProcessors(final EventsLogger.LoggedEvent event,
                                                final List<EventProcessor> eventProcessors) {
        final int      nbOfEventProcessors = eventProcessors.size();
        final double[] additionalValues    = new double[nbOfEventProcessors];

        // Process event through each event processer
        for (int i = 0; i < nbOfEventProcessors; i++) {
            additionalValues[i] = eventProcessors.get(i).process(event);
        }
        return additionalValues;
    }

    /** @return unmodifiable constraints list */
    public List<Constraint> getConstraints() {
        return Collections.unmodifiableList(constraints);
    }

}
