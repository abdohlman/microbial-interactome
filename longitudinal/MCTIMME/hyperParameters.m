classdef hyperParameters
    % object specifying hyperparameters for MC-TIMME model
    % see SI Methods for details regarding hyperparameters
    %
    % Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
    % Copyright (C) 2012 Georg Gerber
    % hyperParameters.m (version 1.00)

    properties
        % hyperparameters for gamma prior on the concentration parameter alpha
        omega_alpha1 = 10e-5;
        omega_alpha2 = 10e-5;
        
        gamma % subject and refOTU specific offset
        
        m_epsilon_1 % prior mean for log epsilon 1
        sigma_epsilon_1 % prior std for log epsilon 1
        m_epsilon_2 % prior mean for log epsilon 2
        sigma_epsilon_2 % prior std for log epsilon 2
        
        m_beta_0 % prior mean for beta_0
        sigma_beta_0 % prior std for beta_0
        v_rho_0 % prior variance for rho_0    
        nu_rho_0 = 5;  % dof for prior on rho_0
       
        v_rho_mu_1 % prior variance for rho_mu_1
        nu_rho_mu_1 = 5;  % dof for prior on rho_mu_1
        v_rho_mu_2 % prior variance for rho_mu_2
        nu_rho_mu_2 = 5;  % dof for prior on rho_mu_1
        
        m_beta_lambda_1 % prior mean for beta_lambda_1
        sigma_beta_lambda_1 % prior std for beta_lambda_1
        v_rho_lambda_1 % prior variance for rho_lambda_1
        nu_rho_lambda_1 = 5; % dof for prior on rho_lambda_1
        v_rho_lambda_2 % prior variance for rho_lambda_2
        nu_rho_lambda_2 = 5; % dof for prior on rho_lambda_2
        lambdaMin = 1; % minimum for relaxation time parameter (in days)
        lambdaMax = 45; % maximum for relaxation time parameter (in days)
        
        omega_pi_c_mu = 1.0; % hyperparameter for prior on c_mu
        omega_pi_c_lambda = 1.0 % hyperparameter for prior on c_lambda
    end
end

