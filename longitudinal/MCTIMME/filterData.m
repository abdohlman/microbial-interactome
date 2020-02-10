function [D_filtered,psi,medianCounts,origOTUIdx,dataIdx,filterDataIdx] = filterData(D_raw)
% remove refOTUs that have small numbers of counts across all time-points
%
% input:
% D_raw = cell array w/ data from all subjects (generated from e.g.,
% loadDethlefsenData.m)
%
% outputs:
% D_filtered = data w/ refOTUs removed that don't pass filtering criteria
% psi = subject and time-point specific offsets based on total reads (see SI methods)
% medianCounts = median of sequencing reads across all subjects and
% time-points
% origOTUIdx,dataIdx,filterDataIdx = indices used for subsequent mapping back to original data
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% filterData.m (version 1.00)

% data filtering criteria
% require refOTU to have >= minCounts at >= minTimepoints
minCounts = 10;
minTimepoints = 5;

% initialize arrays
D_filtered = cell(length(D_raw),1);
psi = cell(length(D_raw),1);
origOTUIdx = cell(length(D_raw),1);
totalCounts = cell(length(D_raw),1);
totalAll = [];

% compute total read counts for each subject at each time-point
for pt=1:length(D_raw),
    totalCounts{pt} = sum(D_raw{pt});
    totalAll = [totalAll totalCounts{pt}];
    origOTUIdx{pt} = 1:size(D_raw{pt},1);
end;
% compute median read counts for each subject
medianCounts = median(totalAll);

% filter data according to criteria
for pt=1:length(D_raw),
    EC = (D_raw{pt} >= minCounts);
    EC = sum(EC,2)';
    filterDataIdx = find(EC >= minTimepoints);
    D_filtered{pt} = D_raw{pt}(filterDataIdx,:);
    origOTUIdx{pt} = origOTUIdx{pt}(filterDataIdx);
    % calculate psi offset for each subject and time-point
    psi{pt} = log(totalCounts{pt}) - log(medianCounts);
end;

dataIdx = zeros(length(D_raw),2);
dataIdx(1,1) = 1;
dataIdx(1,2) = size(D_filtered{1},1);
for pt=2:length(D_raw),
    dataIdx(pt,1) = 1+dataIdx(pt-1,2);
    dataIdx(pt,2) = dataIdx(pt,1)+size(D_filtered{pt},1)-1;
end;