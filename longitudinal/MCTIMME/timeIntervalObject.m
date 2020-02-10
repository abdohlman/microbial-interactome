classdef timeIntervalObject
    % object for storing information on treatment or non-treatment time
    % intervals
    %
    % Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
    % Copyright (C) 2012 Georg Gerber
    % timeIntervalObject.m (version 1.00)

    properties
        treat % 0 = non-treatment interval, 1 = treatment interval
        startTimeIdx % index of time-point to start at
        endTimeIdx % index of time-point to end at
    end
end

