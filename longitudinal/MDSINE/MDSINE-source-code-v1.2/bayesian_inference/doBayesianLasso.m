function [] = doBayesianLasso(cfg,outputDir,T,intervene_matrix,perturbations,data,species_names,BMD,experimentBlocks)
% This file is a part of the MDSINE program.
%
% MDSINE: Microbial Dynamical Systems INference Engine
% Infers dynamical systems models from microbiome time-series datasets, and
% predicts biologically relevant behaviors of the ecosystems.
% Copyright (C) 2015 Vanni Bucci, Georg Gerber
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% run Bayesian inference methods on data

outputFileName = [outputDir 'BAL.mat'];

params = cfg.preprocessing;%readParameterFile(preProcessDataParmsFileName);

% rng('default');

warning off;

% find first index for each experimental block
uniqueExperiments = unique(experimentBlocks);
blockIdx = zeros(length(uniqueExperiments),1);
for e=1:length(blockIdx),
    f = find(experimentBlocks == uniqueExperiments(e));
    blockIdx(e) = f(1);
end;

sigma_biomass_est = 0.5;

biomass = cell(length(experimentBlocks),1);
% infer spline representation for biomass data for each experimental block
if params.useSplines == 1,
    disp('Inferring splines for biomass data');
    for e=1:length(blockIdx),
        ep = find(experimentBlocks == uniqueExperiments(e));
        %[S_sample,BO,BD,biomassScaleFactor] = splineBiomassMCMC(T{blockIdx(e)},BMD(ep),sigma_biomass_est,splineBiomassMCMCParamsFileName);
        [BMD,sigma_biomass_est] = LogBiomass_GetStd(BMD,params);
        [S_sample,BO,BD,biomassScaleFactor] = splineBiomassMCMC(T{blockIdx(e)},BMD(ep),sigma_biomass_est,cfg.bayesianSplineBiomass, params.numReplicates);
        [biomass(ep),biomass_high,biomass_low] = estBiomassFromSplineSamples(S_sample,T{blockIdx(e)},BO,biomassScaleFactor);
    end;
else,
    % average biomass replicates
    biomass = avgBiomass(BMD,params.numReplicates);
    for s=1:length(uniqueExperiments),
        biomass{s} = log(biomass{s});
    end;
end;

% filter out OTUs with low counts
[keep_species,intervene_matrix_filtered,intervene_matrix_filtered_merge,data_counts_filtered,species_names_filtered] = filterData(intervene_matrix,data,species_names,experimentBlocks,cfg.preprocessing);

perturbations_filter = perturbations;

if params.useSplines == 1,
    % calculate normalizing constants
    [total_count_norm,biomass_norm,scaleFactor,med_counts,med_biomass] = calcNorms(keep_species,data,intervene_matrix_filtered,biomass,experimentBlocks);

    % infer spline representation for counts data
    f = cell(length(experimentBlocks),1);
    df_dt = cell(length(experimentBlocks),1);
    disp('Inferring splines for counts data');
    for e=1:length(blockIdx),
        ep = find(experimentBlocks == uniqueExperiments(e));
        %[phi_sample,B,BD] = splineCountsMCMC(T{blockIdx(e)},data_counts_filtered{e},total_count_norm{e},biomass_norm{e},scaleFactor{e},intervene_matrix_filtered{e},splineCountsMCMCParamsFileName,1);
        [phi_sample,B,BD] = splineCountsMCMC(T{blockIdx(e)},data_counts_filtered{e},total_count_norm{e},biomass_norm{e},scaleFactor{e},intervene_matrix_filtered{e},cfg.bayesianSplineCounts,1);
        [f(ep),f_high,f_low,df_dt(ep),df_dt_high,df_dt_low] = estOTUTrajDerivFromSplineSamples(phi_sample,T{blockIdx(e)},B,BD,intervene_matrix_filtered{e},scaleFactor{e});
    end;
else,
    med_counts = 1;
    med_biomass = 0;
    f = cell(length(experimentBlocks),1);
    df_dt = cell(length(experimentBlocks),1);
    for e=1:length(blockIdx),
        ep = find(experimentBlocks == uniqueExperiments(e));
        if ~isempty(perturbations),
            [f(ep),df_dt(ep),intervene_matrix_filtered{e},perturbations_filter(ep)] = estForwardDiffDeriv(T{blockIdx(e)},data(ep),keep_species{e},intervene_matrix_filtered{e},perturbations(ep));
        else,
            [f(ep),df_dt(ep),intervene_matrix_filtered{e}] = estForwardDiffDeriv(T{blockIdx(e)},data(ep),keep_species{e},intervene_matrix_filtered{e},perturbations(ep));
        end;
    end;
end;

% now re-map matrices so that species indices are the same
keep_species_total = [];
for e=1:length(blockIdx),
    keep_species_total = union(keep_species_total,keep_species{e});
end;
species_names_filtered_total = species_names(keep_species_total);
no_keep = setdiff(1:size(data{1},2),keep_species_total);
ivmf = cell(length(experimentBlocks),1);
for e=1:length(blockIdx),
    ep = find(experimentBlocks == uniqueExperiments(e));
    for s=1:length(ep),
        ft = zeros(size(f{ep(s)},1),size(data{1},2));
        ft_deriv = zeros(size(ft));
        ivt = zeros(size(ft));
        ft(:,keep_species{e}) = f{ep(s)};
        ft_deriv(:,keep_species{e}) = df_dt{ep(s)};
        ivt(:,keep_species{e}) = intervene_matrix_filtered{e};

        ft(:,no_keep) = [];
        ft_deriv(:,no_keep) = [];
        ivt(:,no_keep) = [];

        f{ep(s)} = ft;
        df_dt{ep(s)} = ft_deriv;
        ivmf{ep(s)} = ivt;
    end;
end;
intervene_matrix_filtered_merge = ivmf;

% calculate data matrices and estimate hyperparameters
if params.useSplines == 1,
    [otu_scale,X,F_hat_prime,growth_mean_p,growth_std_p,self_reg_mean_p,interact_std_p,deriv_std,perturb_std_p] = estMatricesHyperparameters(f,df_dt,intervene_matrix_filtered_merge,1,perturbations_filter);
else,
    [otu_scale,X,F_hat_prime,growth_mean_p,growth_std_p,self_reg_mean_p,interact_std_p,deriv_std,perturb_std_p] = estMatricesHyperparametersFDiff(f,df_dt,intervene_matrix_filtered_merge,1,perturbations_filter);
end;

numPerturb = 0;
if ~isempty(perturbations_filter),
    for s=1:length(perturbations_filter),
        numPerturb = max(size(perturbations_filter{s},2));
    end;
end;

% infer dynamical system using adaptive Bayesian lasso algorithm
disp('Inferring dynamical system with adaptive Bayesian lasso');
%[Theta_samples_lasso] = MDSINELassoMCMC(X,F_hat_prime,growth_mean_p,growth_std_p,self_reg_mean_p,interact_std_p,perturb_std_p,numPerturb,1,MDSINELassoMCMCParamsFileName);
[Theta_samples_lasso] = MDSINELassoMCMC(X,F_hat_prime,growth_mean_p,growth_std_p,self_reg_mean_p,interact_std_p,perturb_std_p,numPerturb,1,cfg.bayesianLasso);
[Theta_lasso_low,Theta_lasso,Theta_lasso_high,Theta_lasso_indicator] = estInteractMatrixFromSamples(Theta_samples_lasso,0.95);
Theta_lasso_indicator = Theta_lasso_indicator(:,2:size(Theta_lasso,2));

Theta_global = Theta_lasso;
Theta_samples_global = Theta_samples_lasso;
Theta_bayes_factors = [];

mkdir(outputDir);
save(outputFileName, '-v7.3');
