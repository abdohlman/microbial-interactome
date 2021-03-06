function [bestTimes,HM] = experimentalDesign(subjNum,times,Bt,R_inv,numTimePointDesign)
% find set of time-points that optimize information theoretic utility
% function using greedy algorithm
%
% inputs:
% subjNum = subject number on which to perform experimental design
% times = cell array of time-points for each subject
% Bt = WLS matrices for this subject, generated by experimentalDesignMatrices
% R_inv = inverse covariance matrices, generated by experimentalDesignMatrices
% numTimePointDesign = number of time-points to include in final design
%
% outputs:
% bestTimes = array of selected time-points
% HM = array of utility scores for selected time-points
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% experimentalDesign.m (version 1.00)

% get arrays sizes from Bt
numSamples = size(Bt,5);
numOTUs = size(Bt,4);
np = size(Bt,1);

% set up arrays for calculations at daily time-points
itimes = (1:max(times{subjNum}))';
useTimes = (1:length(itimes));

% initialize arrays
HM = [];
bestTimes = [];
IM = zeros(np,np,numOTUs,numSamples);
Bt_t = zeros(np,np,numOTUs,numSamples);
IM_new = zeros(np,np);
for s=1:numSamples,
    for otu=1:numOTUs,
        IM(:,:,otu,s) = R_inv(:,:,s);
    end;
end;

% perform greedy optimization
for nr = 1:numTimePointDesign,
    H = zeros(length(useTimes),1);
    % calculate utility score for each remaining time-point
    for tidx=1:length(useTimes),
        Bt_t = squeeze(sum(Bt(:,:,useTimes(tidx),:,:),3));
        for s=1:numSamples,
            for otu=1:numOTUs,
                % add up contributions from each refOTU
                IM_new = IM(:,:,otu,s) + Bt_t(:,:,otu,s);
                H(tidx) = H(tidx) + log(det(IM_new));
            end;
        end;
        H(tidx) = H(tidx)/numSamples;
        
        if mod(tidx,10)==0,
            fprintf('%i ',tidx);
        end;
    end;
    fprintf('\n');
    % chose time-point w/ maximum score
    [m,mi] = max(H);
   
    % break ties if necessary
    ve = find(H == m);
    if length(ve) > 1 && ~isempty(bestTimes),
        dpv = zeros(length(ve),1);
        for jj=1:length(dpv),
            dpv(jj) = min(abs(bestTimes - itimes(useTimes(ve(jj)))));
        end;
        [mm,mmi] = max(dpv);
        mi = ve(mmi); 
    end;
    
    bestTimes = [bestTimes itimes(useTimes(mi))];    
    
    % update matrices
    for s=1:numSamples,
        Bt_t = squeeze(sum(Bt(:,:,useTimes(mi),:,:),3));
        for otu=1:numOTUs,
            IM(:,:,otu,s) = IM(:,:,otu,s) + Bt_t(:,:,otu,s);
        end;
    end;
    % remove selected time-point from active set
    fprintf('Selected time-point %i\n',itimes(useTimes(mi)));
    useTimes = setdiff(useTimes,useTimes(mi));
    HM = [HM m];
end;
