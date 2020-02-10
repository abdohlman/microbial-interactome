function [NMI,RI] = CSGQualityMeasures(refCSGList,CSGList)
% compute normalized mutual information and Rand index as measures of
% consensus signature group (CSG)
% quality
%
% inputs:
% refCSGList = list s.t. line i gives the 'ground truth' CSG membership for refOTU i
% CSGList = list s.t. line i gives the CSG membership for refOTU in CSGs to
% be tested
%
% outputs:
% NMI = normalized mutual information between refCSGs and CSGs
% RI = Rand index for CSGs (refCSGs used as 'ground truth')
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% CSGQualityMeasures.m (version 1.00)

numRefCSGs = max(refCSGList);
numCSGs = max(CSGList);
numOTUs = length(CSGList);

% entropy of CSGs
H_CSG = 0;
% entropy of refCSG
H_refCSG = 0;

% create cell arrays w/ CSG/refCSG members
CSGs = cell(numCSGs,1);
refCSGs = cell(numRefCSGs,1);

for k=1:numCSGs,
    F = find(CSGList == k);
    CSGs{k} = F;
    % compute entropy
    if ~isempty(F),
        H_CSG = H_CSG - (length(F)/numOTUs)*log(length(F)/numOTUs);
    end;
end;

for k=1:numCSGs,
    F = find(refCSGList == k);
    refCSGs{k} = F;
    % compute entropy
    if ~isempty(F),
        H_refCSG = H_refCSG - (length(F)/numOTUs)*log(length(F)/numOTUs);
    end;
end;

% compute normalized entropy
NMI = 0;
for k=1:numCSGs,
    for j=1:numRefCSGs,
        F_refF = intersect(CSGs{k},refCSGs{j});
        if ~isempty(F_refF),
            NMI = NMI + (length(F_refF)/numOTUs)*log(numOTUs*length(F_refF)/(length(CSGs{k})*length(refCSGs{j})));
        end;
    end;
end;
NMI = NMI/((H_CSG+H_refCSG)/2);

% compute Rand index
RI = 0;
% true positives
TP = 0;
% true negative
TN = 0;
for i=1:(numOTUs-1),
    for j=(i+1):numOTUs,
        % true positive
        if CSGList(i) == CSGList(j) && refCSGList(i) == refCSGList(j),
            TP = TP + 1;
        end;
        % true negative
        if CSGList(i) ~= CSGList(j) && refCSGList(i) ~= refCSGList(j),
            TN = TN + 1;
        end;
    end;
end;
RI = (TP+TN)/(numOTUs*(numOTUs-1)/2);



