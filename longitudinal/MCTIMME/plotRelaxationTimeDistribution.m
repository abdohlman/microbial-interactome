function [RTD] = plotRelaxationTimeDistribution(signatures,hparams,whichRelaxationParam,pickOTUList)
% plot smoothed distribution of relaxation time constants
%
% inputs:
% signatures = cell array of signatures, obtained from
% hparams = hyperParameters object, calculated by estimateHyperParameters
% individualOTUSignatures
% whichRelaxationParam = 1 for the first relaxation time, >1 for the
% subsequent relaxation times
% pickOTUList = optional cell array containing one array per subject of
% refOTUs to include in the plot
%
% output:
% RTD = array containing selected relaxation times
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% plotRelaxationTimeDistribution.m (version 1.00)

numSubjects = length(signatures);
% cell array of OTUs to include in plot
useOTUs = cell(numSubjects,1);

% if pickOTUList isn't specified, use all refOTUs
if nargin == 3,
    for ds=1:numSubjects,
        useOTUs{ds} = 1:length(signatures{ds});
    end;
else,
    useOTUs = pickOTUList;
end;

% initialize array to hold relaxation times, which will be grown
% dynamically
RTD = [];
% now put relaxation times into array
for ds=1:numSubjects,
    for k=1:length(useOTUs{ds}),
        otu = useOTUs{ds}(k);
        RTD = [RTD signatures{ds}{otu}.med_relaxation_params(whichRelaxationParam)];
    end;
end;

% compute kernel smoothing density estimate
[Y,X] = ksdensity(RTD,'support',[0 1.5*max(RTD)]);
% eliminate estimates outside of relaxation time parameter range
f = find(X>hparams.lambdaMax);
X(f) = [];
Y(f) = [];

plot(X,Y,'-b','LineWidth',4,'Color','b');
