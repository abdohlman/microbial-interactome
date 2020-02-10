function [signatures] = individualOTUSignatures(hparams,baseDir,times,timeIntervals)
% using MCMC samples as input, for each refOTU in each subject, output a median signature w/ 95%
% credible interval, and estimates of c_mu and c_lambda state probabilities
%
% inputs:
% D_filtered = cell array with data from all subjects
% hparams = hyperparameters object
% psi = subject & time-point specific offsets, calculated by filterData
% baseDir = directory containing MCMC sample output
% times = cell array of time-points for each subject
% timeIntervals = array specifying time-intervals for each subject
%
% outputs:
% signatures = cell array of signatureObjects
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% individualOTUSignatures.m (version 1.00)

% separate time intervals into non-treatment and treatment (antibiotic)
% intervals
[treatIntervals,noTreatIntervals,treatTimeIdx,noTreatTimeIdx] = separateIntervals(timeIntervals);
numMuParams = length(noTreatIntervals{1});
numGLMParams = length(timeIntervals{1});
numDeltaLambdaParams = numMuParams-2;
num_c_mu_states = 2^(numMuParams-1);
num_c_lambda_states = 2^numDeltaLambdaParams;

% # of signature parameters includes equilibrium level,
% relaxation time, and adaptive complexity parameters
numSignatureParameters = length(timeIntervals{1})+1+1;
if numDeltaLambdaParams>0,
    numSignatureParameters = numSignatureParameters + numDeltaLambdaParams+ 1;
end;

numSubjects = length(hparams.gamma);

% read in assignments of refOTU time-series to signatures
signatureAssignsFileN = cell(numSubjects,1);
itimes = cell(numSubjects,1);
itimeIntervals = cell(numSubjects,1);
for ds=1:numSubjects,
    signatureAssignsFileN{ds} = [baseDir '_signatureAssigns' int2str(ds) '.txt'];
    signatureAssigns{ds} = dlmread(signatureAssignsFileN{ds});
    numOTUs(ds) = length(hparams.gamma{ds});
    itimes{ds} = (1:max(times{ds}))';
    tis = cell(length(timeIntervals{ds}),1);
    for intv=1:length(timeIntervals{ds}),
        ti = timeIntervals{ds}{intv};
        ti2 = timeIntervalObject;
        ti2.treat = ti.treat;
        ti2.startTimeIdx = times{ds}(ti.startTimeIdx);
        ti2.endTimeIdx = times{ds}(ti.endTimeIdx);
        tis{intv} = ti2;
    end;
    itimeIntervals{ds} = tis;
end;

% get # MCMC samples, and # signatures in each MCMC sample
hyperParamFileN = [baseDir '_hyperParams.txt'];
H = dlmread(hyperParamFileN);
numSignatures = H(:,1);
numSamples = size(H,1);

% read in signature parameters from each MCMC sample
signatureParametersFileN = [baseDir '_signatureParameters.txt'];
samples = cell(numSamples,1);
fid = fopen(signatureParametersFileN,'r');
for i=1:numSamples,
    tline = fgetl(fid);
    C = sscanf(tline,'%f');
    nc = numSignatures(i)*numSignatureParameters;
    C = C(1:nc);
    samples{i} = reshape(C,numSignatureParameters,numSignatures(i))';
end;
fclose(fid);

% intialize cell arrays to store signature parameters and adaptive complexity probability
% estimates
signatures = cell(numSubjects,1);

% for each subject & refOTU, get all samples
for ds=1:numSubjects,
    % temp arrays to hold signatures
    tj = cell(numOTUs(ds),1);
    for otu=1:numOTUs(ds),
        % temp variable to store trajectory (log NBD mean) from each sample
        M = zeros(numSamples,length(itimes{ds}));
        % temp variable to store relaxation time parameters from each
        % sample
        L = zeros(numSamples,numDeltaLambdaParams+1);
        % temp variable to store c_mu setting for each sample
        CP = zeros(1,num_c_mu_states);
        ULDP = [];
        if numDeltaLambdaParams > 0,
            % temp variable to store c_lambda setting for each sample
            ULDP = zeros(1,num_c_lambda_states);
        end;
        for s=1:numSamples,
            k = signatureAssigns{ds}(s,otu);
            C = samples{s};
            % equilibrium level parameters
            Xt = C(k,1:numGLMParams);
            % relaxation time parameters
            lambda_1 = C(k,numGLMParams+1);
            delta_lambda = [];
            if ~isempty(ULDP),
                delta_lambda = C(k,(numGLMParams+2):(numGLMParams+2+numDeltaLambdaParams-1));
            end;
            
            % calculate design matrix for given time-points sampled,
            % time-intervals and relaxation time parameters
            B = designMatrix(itimes{ds},itimeIntervals{ds},lambda_1,delta_lambda,hparams.lambdaMin,hparams.lambdaMax);
            % calculate log GLM mean based on design matrix and parameters
            % for all time-points (frequency of 1 day)
            M(s,:) = (B*Xt')';
            % transform lambda values into space of relaxation times
            % parameters
            lambda1 = (exp(lambda_1)/(1+exp(lambda_1)))*(hparams.lambdaMax-hparams.lambdaMin) + hparams.lambdaMin;
            lambda2 = [];
            if ~isempty(ULDP),
                currentLambda = lambda_1;
                lambda_2 = zeros(1,length(delta_lambda));
                for dl=1:length(delta_lambda),
                    currentLambda = currentLambda+delta_lambda(dl);
                    lambda_2(dl) = (exp(currentLambda)/(1+exp(currentLambda)))*(hparams.lambdaMax-hparams.lambdaMin) + hparams.lambdaMin;
                end;
            end;
            % concatenate lambda parameters
            L(s,1) = lambda1;
            if ~isempty(ULDP),
                L(s,2:(numDeltaLambdaParams+1)) = lambda_2;
            end;
            c_lambda = [];
            if numDeltaLambdaParams>0,
                c_mu = C(k,numSignatureParameters-1);
                c_lambda = C(k,numSignatureParameters);
            else,
                c_mu = C(k,numSignatureParameters);
            end;
            % count # of times each c_mu or c_lambda state occurs
            CP(c_mu) = CP(c_mu)+1;
            if ~isempty(ULDP),
                ULDP(c_lambda) = ULDP(c_lambda)+1;
            end;
        end;
        M = prctile(M,[2.5 50 97.5]);
        L = prctile(L,[2.5 50 97.5]);
        CP = CP/(numSamples);
        if ~isempty(ULDP),
            ULDP = ULDP/(numSamples);
        end;
        
        % instantiate signature object
        ts = signatureObject;
        ts.q025_trajectory = M(1,:);
        ts.med_trajectory = M(2,:);
        ts.q975_trajectory = M(3,:);
        
        if ~isempty(ULDP),
            ts.q025_relaxation_params = L(1,:);
            ts.med_relaxation_params = L(2,:);
            ts.q975_relaxation_params = L(3,:);
        else,
            ts.q025_relaxation_params = L(1);
            ts.med_relaxation_params = L(2);
            ts.q975_relaxation_params = L(3);
        end;
        ts.c_mu_probs = CP;
        ts.c_lambda_probs = ULDP;
        
        tj{otu} = ts;
        
        os = sprintf('subject %i refOTU %i',ds,otu);
        disp(os);
    end;
    
    signatures{ds} = tj;
end;