function [h] = plotIndividualOTUSignature(subjNum,otuNum,signatures,hparams,times,timeIntervals,D_filtered,psi,medianCounts,origOTUIdx)
% plot individual signature for a refOTU, including median & 95% CI +
% scaled counts from original data
%
% inputs:
% subjNum = subject number
% otuNum = # of refOTU in subject
% signatures = cell array of signatureObjects, calculated by
% individualOTUSignatures
% hparams = hyperParameters object, calculated by estimateHyperParameters
% times = cell array of time-points for each subject
% timeIntervals = array specifying time-intervals for each subject
% D_filtered = filtered data, calculated by filterData
% psi = subject & time-point specific offsets, calculated by filterData
% medianCounts = median counts from data across all subjects & time-points,
% calculated by filterData
% origOTUIdx = indices of refOTUs so that we can map back to original #'s,
% calculated by filterData
% 
% output:
% h = figure handle
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% plotIndividualOTUSignature.m (version 1.00)

% vector of times w/ daily sampling
itimes = (1:max(times{subjNum}))';

h = figure;
hold on;

% create trajectory & 95% CI, and add subject & refOTU specific offset
Y = [signatures{subjNum}{otuNum}.q025_trajectory ; signatures{subjNum}{otuNum}.med_trajectory ; signatures{subjNum}{otuNum}.q975_trajectory];
Y = Y + repmat(hparams.gamma{subjNum}(otuNum),3,length(itimes));

% create scaling factor, so total # of reads per sample is scaled to
% 10K
scaleF = medianCounts/10000;

% transform trajectory into count space
Y = exp(Y)/scaleF;

% set up colored patches to render 95% CI
col1 = [1.0 0 0];
edgeColor1 = col1;
patchSaturation1=0.40;
faceAlpha1 = 1.0;
patchColor1=col1+(1-col1)*(1-patchSaturation1);
uE1=Y(3,:);
lE1=Y(1,:);
yP1=[lE1,fliplr(uE1)];
xP1=[itimes',fliplr(itimes')];
patch(xP1,yP1,1,'facecolor',patchColor1,...
              'edgecolor','none',...
              'facealpha',faceAlpha1);
hold on;
plot(itimes,lE1,'-','color',edgeColor1,'LineWidth',2);
plot(itimes,uE1,'-','color',edgeColor1,'LineWidth',2);

% get max of 95% CI and data, so can set Y-axis limits
maxL = max(Y(3,:));
maxL2 = max(D_filtered{subjNum}(otuNum,:)./(exp(psi{subjNum})*scaleF));
maxL = max(maxL2,maxL)+10;

% plot vertical lines to show antibiotic treatment intervals
for intv=1:length(timeIntervals{subjNum}),
    ti = timeIntervals{subjNum}{intv};
    if ti.treat == 1,
        plot([times{subjNum}(ti.startTimeIdx) times{subjNum}(ti.startTimeIdx)],[0 maxL],'-b','LineWidth',2);
        plot([times{subjNum}(ti.endTimeIdx) times{subjNum}(ti.endTimeIdx)],[0 maxL],'-b','LineWidth',2);
    end;
end;

% plot median trajectory
plot(itimes,Y(2,:),'--r','LineWidth',4);

titl = sprintf('subj#%i refOTU#%i',subjNum,origOTUIdx{subjNum}(otuNum));
title(titl,'FontName', 'Arial', 'FontSize', 14);
ylabel('normalized counts/10K reads','FontName', 'Arial', 'FontSize', 20);
xlabel('time','FontName', 'Arial', 'FontSize', 20);

% plot data
plot(times{subjNum},D_filtered{subjNum}(otuNum,:)./(exp(psi{subjNum})*scaleF),'-k');
plot(times{subjNum},D_filtered{subjNum}(otuNum,:)./(exp(psi{subjNum})*scaleF),'ok','MarkerSize',8,'LineWidth',2);

axis tight;

set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 20);
set(gca, 'LineWidth',3);