function [timeIntervals] = readTimeIntervalFile(fileName,timePoints)
% read tab-delimited file containing time intervals
% each line specifies a time-interval
% column 1 = no_treat or treat (specifies whether this is a non-treatment or a
% treatment interval)
% column 2 = first time-point in the interval
% column 3 = last time-point in the interval
%
% inputs:
% fileName = name of file containing time interval data
% timePoints = vector of time-points at which data has been sampled
%
% output:
% timeIntervals = array of timeIntervalObjects
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% readTimeIntervalFile.m (version 1.00)

fid = fopen(fileName);
C = textscan(fid,'%s %u %u');
fclose(fid);

numIntervals = length(C{1});
timeIntervals = cell(numIntervals,1);
for i=1:numIntervals,
    ti = timeIntervalObject;
    switch lower(char(C{1}(i))),
        case 'no_treat'
            ti.treat = 0;
        case 'treat'
            ti.treat = 1;
    end;
    startTime = C{2}(i);
    endTime = C{3}(i);
    % find indices of starting and ending times w/in vector of time-points
    ti.startTimeIdx = find(timePoints == startTime);
    ti.endTimeIdx = find(timePoints == endTime);
    timeIntervals{i} = ti;
end;

