function [SD1_mu,SD1_lambda,SD2,SD2I,SD2D] = signatureDiversityScores(psi,hparams,timeIntervals,baseDir)
% compute signature diversity scores
%
% input:
% psi = subject & time-point specific offsets, calculated by filterData
% hparams = hyperParameters object, calculated by estimateHyperParameters
% timeIntervals = cell array specifying time-intervals for each subject
% baseDir = base directory for MCMC samples
%
% output:
% SD1_mu = cell array per subject of credible intervals for SD1_mu scores
% SD1_lambda = cell array per subject of credible intervals for SD1_lambda
% scores
% SD2 = cell array per subject of credible intervals for SD2 scores
% SD2I = credible interval for SD3 independent score component
% SD2D = credible interval for SD3 depedent score component
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% signatureDiversityScores.m (version 1.00)

numSubjects = length(hparams.gamma);

% separate time intervals into non-treatment and treatment (antibiotic)
% intervals
[treatIntervals,noTreatIntervals,treatTimeIdx,noTreatTimeIdx] = separateIntervals(timeIntervals);
numMuParams = length(noTreatIntervals{1});
numDeltaLambdaParams = numMuParams-2;

% # of signature parameters includes equilibrium level,
% relaxation time, and adaptive complexity parameters
numSignatureParameters = length(timeIntervals{1})+1+1;
if numDeltaLambdaParams>0,
    numSignatureParameters = numSignatureParameters + numDeltaLambdaParams+ 1;
end;

% load assignment of each refOTU time-series in each subject to a
% signature
signatureAssignsFileN = cell(numSubjects,1);
for ds=1:numSubjects,
    signatureAssignsFileN{ds} = [baseDir '_signatureAssigns' int2str(ds) '.txt'];
    signatureAssigns{ds} = dlmread(signatureAssignsFileN{ds});
    numOTUs(ds) = length(hparams.gamma{ds});
end;

% determine # MCMC samples & # signatures per sample
hyperParamFileN = [baseDir '_hyperParams.txt'];
H = dlmread(hyperParamFileN);
numSamples = size(H,1);
numSignatures = H(:,1);

% load signature parameters from each MCMC sample
samples = cell(numSamples,1);
signatureParametersFileN = [baseDir '_signatureParameters.txt'];
fid = fopen(signatureParametersFileN,'r');
for i=1:numSamples,
    tline = fgetl(fid);
    C = sscanf(tline,'%f');
    nc = numSignatures(i)*numSignatureParameters;
    C = C(1:nc);
    samples{i} = reshape(C,numSignatureParameters,numSignatures(i))';
end;
fclose(fid);

% initialize cell arrays for signature diversity scores
SD1_mu = cell(numSubjects,1);
SD1_lambda = cell(numSubjects,1);
SD2 = cell(numSubjects,1);
SD2D = credibleInterval;
SD2I = credibleInterval;

SD2I_t = zeros(1,numSamples);
SD2D_t = zeros(1,numSamples);

for ds=1:numSubjects,
    % temporary arrays to store score for each sample, so we can then
    % compute credible intervals from the score distributions
    SD1_mu_t = zeros(1,numSamples);
    SD1_lambda_t = zeros(1,numSamples);
    SD2_t = zeros(1,numSamples);    
    for s=1:numSamples,
        % distribution of signature assignments
        signatureAssignDist = zeros(1,numSignatures(s));
        % count of refOTUs w/ c_mu > 1 
        count_c_mu = 0;
        % count of refOTUs w/ c_lambda > 1 
        count_c_lambda = 0;
        for otu=1:numOTUs(ds),
            % signature to which refOTU is assigned
            k = signatureAssigns{ds}(s,otu);
            % increment # assignments signature k
            signatureAssignDist(k) = signatureAssignDist(k) + 1;
            % get c_mu and c_lambda states from the MCMC sample
            C = samples{s};
            c_lambda = [];
            if numDeltaLambdaParams>0,
                c_mu = C(k,numSignatureParameters-1);
                c_lambda = C(k,numSignatureParameters);
            else,
                c_mu = C(k,numSignatureParameters);
            end;
            % increment counts for c_mu & c_lambda
            if c_mu > 1,
                count_c_mu = count_c_mu+1;
            end;
            if numDeltaLambdaParams>0,
                if c_lambda > 1,
                    count_c_lambda = count_c_lambda+1;
                end;
            end;
        end;
        % SD1_mu is the fraction of refOTUs w/ c_mu > 1
        SD1_mu_t(s) = count_c_mu/numOTUs(ds);
        if numDeltaLambdaParams>0,
            % SD1_lambda is the fraction of refOTUs w/ c_lambda > 1
            SD1_lambda_t(s) = count_c_lambda/numOTUs(ds);
        end;
        % SD2 is the expected equil # of signatures per 100 refOTUs
        SD2_t(s) = exp(entropy(signatureAssignDist))*100/numOTUs(ds);
        % SD2I score is the wt'd avg of SD2 scores
        SD2I_t(s) = SD2I_t(s) + SD2_t(s)*numOTUs(ds)/sum(numOTUs);
    end;
    % compute credible intervals for SD1 & SD2 scores
    SD1_mu{ds} = calcCI(SD1_mu_t);
    if numDeltaLambdaParams>0,
        SD1_lambda{ds} = calcCI(SD1_lambda_t);
    end;
    SD2{ds} = calcCI(SD2_t);
end;
% calculate SD2I credible intervals
SD2I = calcCI(SD2I_t);

% calculate SD2D score, which is the equil # of signatures per 100 refOTUs
% of the hypothetical 'combined ecosystem'
for s=1:numSamples,
    signatureAssignDist = zeros(1,numSignatures(s));
    % distribution of signature assignments
    for ds=1:numSubjects,
         for otu=1:numOTUs(ds),
            % signature to which refOTU is assigned in 'combined ecosystem'
            k = signatureAssigns{ds}(s,otu);
            % increment # assignments signature k
            signatureAssignDist(k) = signatureAssignDist(k) + 1;
         end;
    end;
    SD2D_t(s) = exp(entropy(signatureAssignDist))*100/sum(numOTUs);
end;
% calculate SD2D credible intervals
SD2D = calcCI(SD2D_t);

function [H] = entropy(p)
    pn = p/sum(p);
    f = find(p > 0);
    H = -sum(pn(f).*log(pn(f)));

function [CI] = calcCI(v)
    tv = prctile(v,[2.5 50 97.5]);
    CI = credibleInterval;
    CI.q025 = tv(1);
    CI.med = tv(2);
    CI.q975 = tv(3); 