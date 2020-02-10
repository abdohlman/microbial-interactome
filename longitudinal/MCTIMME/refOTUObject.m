classdef refOTUObject
    % object for storing attributes for reference operational taxonomic unit (refOTU)
    %
    % Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
    % Copyright (C) 2012 Georg Gerber
    % refOTUObject.m (version 1.00)
    properties
        subjNum % subject number
        OTUNum % relative number of OTU (in each subject, in filtered data)
        origOTU % absolute number of OTU (across all subjects, in unfiltered data)
    end  
end

