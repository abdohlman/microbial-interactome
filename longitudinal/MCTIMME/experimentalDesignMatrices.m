function [Bt,R_inv] = experimentalDesignMatrices(subjNum,hparams,baseDir,times,timeIntervals)
% generate WLS and covariance matrices needed for automated experimental
% design algorithm (experimentalDesign.m)
%
% inputs:
% subjNum = subject number to calculate matrices
% hparams = hyperparameters object
% baseDir = directory containing MCMC sample output
% times = cell array of time-points for each subject
% timeIntervals = array specifying time-intervals for each subject
%
% outputs:
% Bt = array of information matrices, one for each time-point, refOTU, and
% MCMC sample
% R_iv = array of inverse covariance matrices, one for each sample
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% experimentalDesignMatrices.m (version 1.00)

% separate time intervals into non-treatment and treatment (antibiotic)
% intervals
[treatIntervals,noTreatIntervals,treatTimeIdx,noTreatTimeIdx] = separateIntervals(timeIntervals);
numMuParams = length(noTreatIntervals{1});
numGLMParams = length(timeIntervals{subjNum});
numDeltaLambdaParams = numMuParams-2;
num_c_mu_states = 2^(numMuParams-1);
num_c_lambda_states = 2^numDeltaLambdaParams;

% # of signature parameters includes equilibrium level,
% relaxation time, and adaptive complexity parameters
numSignatureParameters = length(timeIntervals{1})+2;
if numDeltaLambdaParams>0,
    numSignatureParameters = numSignatureParameters +numDeltaLambdaParams+ 1;
end;

% read assignments of refOTU times series to signatures for this subject
signatureAssignsFileN = sprintf('%s_signatureAssigns%i.txt',baseDir,subjNum);
signatureAssigns = dlmread(signatureAssignsFileN);

% read hyperparameter and signature parameter samples
hyperParamFileN = [baseDir '_hyperParams.txt'];
signatureParametersFileN = [baseDir '_signatureParameters.txt'];
H = dlmread(hyperParamFileN);
numSignatures = H(:,1);
numSamples = size(H,1);
samples = cell(1,1);

% set up arrays for calculations at daily time-points
itimes = (1:max(times{subjNum}))';
itimeIntervals = cell(length(timeIntervals{subjNum}),1);
inoTreatTimes = [];
iTreatTimes = [];
for intv=1:length(timeIntervals{subjNum}),
    ti = timeIntervals{subjNum}{intv};
    ti2 = timeIntervalObject;
    ti2.treat = ti.treat;
    ti2.startTimeIdx = times{subjNum}(ti.startTimeIdx);
    ti2.endTimeIdx = times{subjNum}(ti.endTimeIdx);
    if ti2.treat == 0,
        inoTreatTimes = [inoTreatTimes ti2.startTimeIdx:ti2.endTimeIdx];
    else,
        iTreatTimes = [iTreatTimes ti2.startTimeIdx:ti2.endTimeIdx];
    end;
    itimeIntervals{intv} = ti2;
end;

% frequency to compute matrices from MCMC samples (e.g., sampleFreq = 10 means use every 10th
% sample)
sampleFreq = 10;

% read in signature parameters
fid = fopen(signatureParametersFileN,'r');
i = 0;
sm = [];
for ii=1:numSamples,
    tline = fgetl(fid);
    C = sscanf(tline,'%f');
    if mod(ii,sampleFreq) == 0,
        i = i + 1;
        nc = numSignatures(ii)*numSignatureParameters;
        C = C(1:nc);
        samples{i} = reshape(C,numSignatureParameters,numSignatures(ii))';
        sm = [sm ii];
    end;
end;
fclose(fid);
numSamples = i;

% np = # of equilibrium levels + relaxation time parameters
np = numSignatureParameters-1;
if numDeltaLambdaParams>0,
    np = np - 1;
end;

numOTUs = length(hparams.gamma{subjNum});

% initialize covariance and WLS matrices
R_inv = zeros(np,np,numSamples);
Bt = zeros(np,np,length(itimes),numOTUs,numSamples);

for sr=1:numSamples,
    sample = samples{sr};
    s= sm(sr);
    
    % get hyperparameters
    epsilon1 = exp(H(s,3));
    epsilon2 = exp(H(s,4));
    beta_0 = H(s,6);
    rho_mu_1 = H(s,7);
    rho_mu_2 = H(s,8);
    rho_lambda_1 = H(s,10);
    if numDeltaLambdaParams>0,
        rho_lambda_2 = H(s,11);
    end;
    
    % set R_inv matrix
    vp = zeros(1,np);
    vp(1) = beta_0^2;
    for intv=2:length(timeIntervals{subjNum}),
        if timeIntervals{subjNum}{intv}.treat == 0,
            vp(intv) = rho_mu_1^2;
        else,
            vp(intv) = rho_mu_2^2;
        end;
    end;
    vp(numGLMParams+1) = rho_lambda_1^2;
    if numDeltaLambdaParams>0,
        vp((numGLMParams+2):(numGLMParams+2+numDeltaLambdaParams-1)) = rho_lambda_2^2;
    end;
    
    R_inv(:,:,sr) = diag(1./vp);
    
    % calculate WLS matrices for each refOTU
    for otu=1:numOTUs,
        % get index of signature to which this refOTU time-series was
        % assigned for this MCMC sample
        k = signatureAssigns(s,otu);
        % get signature parameters
        X = sample(k,1:numGLMParams)';
        lambda = sample(k,numGLMParams+1);
        delta_lambda = [];
        if numDeltaLambdaParams > 0,
            delta_lambda = sample(k,(numGLMParams+2):(numGLMParams+2+numDeltaLambdaParams-1));
        end;
        
        % calculate design matrix
        B = designMatrix(itimes,itimeIntervals,lambda,delta_lambda,hparams.lambdaMin,hparams.lambdaMax);
        % log NBD mean (adding back refOTU specific effect)
        nu = (B*X)' + hparams.gamma{subjNum}(otu);
        enu = exp(nu);
        
        % calculate weights for WLS
        wt = zeros(size(enu));      
        wt(inoTreatTimes) = (enu(inoTreatTimes))./(epsilon1*enu(inoTreatTimes)+1);
        wt(iTreatTimes) = (enu(iTreatTimes))./(epsilon2*enu(iTreatTimes)+1);
        
        % calculate information matrix
        B = InfoMatrix(X,itimes,itimeIntervals,lambda,delta_lambda,hparams);
        
        % calculate WLS matrix at each time-point
        for t=1:length(itimes),
            Bt(:,:,t,otu,sr) = Bt(:,:,t,otu,sr) + B(t,:)'*B(t,:)*wt(t);
        end;
    end;
    disp(sprintf('using MCMC sample %i',s));
end;

function [B] = InfoMatrix(X,times,timeIntervals,lambda,delta_lambda,hparams)
% calculate Fisher information matrix
% inputs:
% X = array of equilibrium levels/increments (GLM parameters)
% times = array of time-points at which data has been sampled
% timeIntervals = array of time-points that specify intervals of antibiotic
% exposure or interim periods
% lambda = transformed relaxation time on first exposure interval
% delta_lambda = change in transformed relaxation time from first -> second exposure
% interval
% output:
% B = Fisher information matrix 

numGLMParams = length(X);
numParams = numGLMParams + 1;
if ~isempty(delta_lambda),
    numParams = numParams + length(delta_lambda);
end;

B = zeros(length(times),numParams);
% running total for relaxation time parameter
currentLambda = lambda;

% loop through time intervals and compute design matrix entries
numNoTreat = 1;
for intv=1:length(timeIntervals),
    % non-treatment intervals
    if timeIntervals{intv}.treat == 0,
        % initial interval
        if intv==1,
            for i=timeIntervals{intv}.startTimeIdx:timeIntervals{intv}.endTimeIdx,
                B(i,1) = 1.0;
            end;
        % subsequent non-treatment intervals, which are described by relaxation
        % processes
        else,
            % calculate relaxation time parameter for this interval
            if numNoTreat > 2,
                currentLambda = currentLambda+delta_lambda(numNoTreat-2);
            end;
            % transform into actual relaxation time
            lambda_int = (exp(currentLambda)/(1+exp(currentLambda)))*(hparams.lambdaMax-hparams.lambdaMin) + hparams.lambdaMin;
           
            Bt = zeros(1,numParams);
            % contribution from previous non-treatment intervals to
            % equilibrium level on this interval
            for intv2=1:(intv-1),
                if timeIntervals{intv2}.treat == 0,
                    Bt(intv2) = 1.0;
                end;
            end;
            
            prevInt = timeIntervals{intv-1}.treat;
            li = timeIntervals{intv}.endTimeIdx-timeIntervals{intv}.startTimeIdx+1;
            B(timeIntervals{intv}.startTimeIdx:timeIntervals{intv}.endTimeIdx,:) = repmat(Bt,li,1);
            dt = (times(timeIntervals{intv}.startTimeIdx:timeIntervals{intv}.endTimeIdx)-times(timeIntervals{intv}.startTimeIdx));
            md = exp(-(1/lambda_int)*(dt));
            B(timeIntervals{intv}.startTimeIdx:timeIntervals{intv}.endTimeIdx,intv) = 1.0-md;
            if prevInt==1,
                B(timeIntervals{intv}.startTimeIdx:timeIntervals{intv}.endTimeIdx,intv-1) = md;
            end;
            % derivatives w.r.t. relaxation time parameter(s)
            dv = X(intv-1).*dt.*exp(-(1/lambda_int).*dt).*(lambda_int-hparams.lambdaMin)/((1+exp(currentLambda))*lambda_int^2);
            dv = dv - X(intv).*dt.*exp(-(1/lambda_int).*dt).*(lambda_int-hparams.lambdaMin)/((1+exp(currentLambda))*lambda_int^2);
            B(timeIntervals{intv}.startTimeIdx:timeIntervals{intv}.endTimeIdx,numGLMParams+1)=dv;
            if numNoTreat>2,
                B(timeIntervals{intv}.startTimeIdx:timeIntervals{intv}.endTimeIdx,(numGLMParams+2):(numGLMParams+2+length(delta_lambda)-1))=repmat(dv,1,length(delta_lambda));
            end;
        end;
        numNoTreat = numNoTreat+1;
    end;
    
    % treatment intervals
    if timeIntervals{intv}.treat == 1,
        Bt = zeros(1,numParams);
        % contribution from delta_mu on this interval
        Bt(intv) = 1.0;
        for intv2=1:(intv-1),
            % contribution from previous intervals
            if timeIntervals{intv2}.treat == 0,
                Bt(intv2) = 1.0;
            end;
        end;
        li = timeIntervals{intv}.endTimeIdx-timeIntervals{intv}.startTimeIdx+1;
        B(timeIntervals{intv}.startTimeIdx:timeIntervals{intv}.endTimeIdx,:) = repmat(Bt,li,1);
    end;
end;