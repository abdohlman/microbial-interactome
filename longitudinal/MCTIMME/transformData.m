function [B_D,B_W] = transformData(found,np,B,Y,psiM,gammaM,nu,enu,epsilon1,epsilon2,times1,times2)
% transform data for WLS estimate (see SI Methods)
%
% inputs:
% found = set of refOTU time-series indices that belong to a particular
% signature
% np = number of parameters Generalized Linear Model (GLM) parameters
% B = design matrix
% Y = complete data matrix
% psiM = matrix of psi estimates (subject & time-point offset)
% gammaM = matrix of gamma estimates (subject & refOTU specific offset)
% nu = array of log NBD means for refOTU time-series
% enu = array of NBD means for refOTU time-series
% epsilon1 = parameter controlling NBD variance on non-antibiotic exposure
% intervals
% epsilon2 = parameter controlling NBD variance on antibiotic exposure
% intervals
% times1 = array of non-antibiotic exposure times
% times2 = array of antibiotic exposure times
%
% outputs:
% B_D = offset vector (np x 1)
% B_W = weight matrix (np X np)
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% transformData.m (version 1.00)

numSubjects = length(B);
B_W = zeros(np,np);
B_D = zeros(np,1);

for ds=1:numSubjects,
    if ~isempty(found{ds}),
        Y_hat = ((Y{ds}(found{ds},:)-enu{ds}(found{ds},:))./enu{ds}(found{ds},:))+nu{ds}(found{ds},:)-gammaM{ds}(found{ds},:)-psiM{ds}(found{ds},:);
        wt = zeros(size(Y_hat));
        
        wt(:,times1{ds}) = (enu{ds}(found{ds},times1{ds}))./(epsilon1*enu{ds}(found{ds},times1{ds})+1);
        wt(:,times2{ds}) = (enu{ds}(found{ds},times2{ds}))./(epsilon2*enu{ds}(found{ds},times2{ds})+1);

        W = reshape(wt',size(wt,1)*size(wt,2),1);
        Y_hat = reshape(Y_hat',size(Y_hat,1)*size(Y_hat,2),1);
        
        Bt = repmat(B{ds},length(found{ds}),1);
        Bt_trans = Bt';
        Wt = repmat(W,1,np);
 
        B_W = B_W + Bt_trans*(Wt.*Bt);
        B_D = B_D + Bt_trans*(W.*Y_hat);
    end;
end;