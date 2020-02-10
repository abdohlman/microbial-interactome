function [B] = designMatrix(times,timeIntervals,lambda,delta_lambda,lambdaMin,lambdaMax)
% calculate design matrix
%
% inputs:
% times = array of time-points at which data has been sampled
% timeIntervals = array of timeIntervalObjects that specify intervals of antibiotic
% exposure or interim periods
% lambda = transformed relaxation time on first exposure interval
% delta_lambda = change in transformed relaxation time from first -> second exposure
% interval
% lambdaMin = lower bound in relaxation time constant
% lambdaMax = upper bound on relaxation time constant
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% designMatrix.m (version 1.00)

% number of equilibrium level coefficients
numMuCoefficients = length(timeIntervals);
% initialize coefficient matrix
B = zeros(length(times),numMuCoefficients);
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
            lambda_int = (exp(currentLambda)/(1+exp(currentLambda)))*(lambdaMax-lambdaMin) + lambdaMin;
           
            Bt = zeros(1,numMuCoefficients);
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
            md = exp(-(1/lambda_int)*(times(timeIntervals{intv}.startTimeIdx:timeIntervals{intv}.endTimeIdx)-times(timeIntervals{intv}.startTimeIdx)));
            md = md';
            B(timeIntervals{intv}.startTimeIdx:timeIntervals{intv}.endTimeIdx,intv) = 1.0-md;
            if prevInt==1,
                B(timeIntervals{intv}.startTimeIdx:timeIntervals{intv}.endTimeIdx,intv-1) = md;
            end;
        end;
        numNoTreat = numNoTreat+1;
    end;
    
    % treatment intervals
    if timeIntervals{intv}.treat == 1,
        Bt = zeros(1,numMuCoefficients);
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
