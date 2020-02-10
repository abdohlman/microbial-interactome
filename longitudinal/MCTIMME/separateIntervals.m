function [treatIntervals,noTreatIntervals,treatTimeIdx,noTreatTimeIdx] = separateIntervals(timeIntervals)
% helper function to separate time intervals into treatment and
% non-treatment times
%
% inputs:
% timeIntervals = cell array of timeIntervalObjects for all subjects
%
% outputs:
% treatIntervals = cell array containing all treatment intervals
% noTreatIntervals = cell array containing all non-treatment intervals
% treatTimeIdx = cell array containing indices of all time-points on
% treatment intervals
% noTreatTimeIdx = cell array containing indices of all time-points on
% non-treatment intervals
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% seperateIntervals.m (version 1.00)

numSubjects = length(timeIntervals);
treatTimeIdx = cell(numSubjects,1);
noTreatTimeIdx = cell(numSubjects,1);
numIntervals = length(timeIntervals{1});
treatIntervals = cell(numSubjects,1);
noTreatIntervals = cell(numSubjects,1);

for ds=1:numSubjects,
    treat = [];
    notreat = [];
    numTreat = 1;
    numNoTreat = 1;
    for i=1:numIntervals,
        tidx = timeIntervals{ds}{i}.startTimeIdx:timeIntervals{ds}{i}.endTimeIdx;
        switch timeIntervals{ds}{i}.treat,
            case 0
                notreat = [notreat tidx];
                noTreatIntervals{ds}{numNoTreat} = timeIntervals{ds}{i};
                numNoTreat = numNoTreat + 1;
            case 1
                treat = [treat tidx];
                treatIntervals{ds}{numTreat} = timeIntervals{ds}{i};
                numTreat = numTreat + 1;
        end;
    end;
    treatTimeIdx{ds} = treat;
    noTreatTimeIdx{ds} = notreat;
end;