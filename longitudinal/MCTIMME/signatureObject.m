classdef signatureObject
    % object for holding signature estimates (either individual or consensus)
    %
    % Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
    % Copyright (C) 2012 Georg Gerber
    % signatureObject.m (version 1.00)
    properties
        med_trajectory % median trajectory (log NBD mean)
        med_relaxation_params % median relaxation time parameters (lambda_c and lambda_e)
        q025_trajectory % 2.5% quartile trajectory
        q025_relaxation_params % 2.5% quartile relaxation parameters
        q975_trajectory % 97.5% quartile trajectory
        q975_relaxation_params % 97.5% quartile relaxation parameters
        c_mu_probs % probability for each c_mu state
        c_lambda_probs % probability for each c_lambda_state
    end
end

