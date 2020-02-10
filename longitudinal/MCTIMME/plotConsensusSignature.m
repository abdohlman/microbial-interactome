function [h] = plotConsensusSignature(consensusSignature,refSubj,times,timeIntervals)
% plot consensus signature, including median & 95% CI +
% scaled counts from original data
%
% inputs:
% consensusSignature = a signatureObject
% refSubj = the reference subject for common time-line
% times = cell array of time-points for each subject
% timeIntervals = array specifying time-intervals for each subject
% output:
% h = figure handle
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% plotCensensusSignature.m (version 1.00)

% vector of days w/ daily sampling
itimes = (1:max(times{refSubj}))';

h = figure;
hold on;

% create trajectory & 95% CI, and add subject & refOTU specific offset
Y = [consensusSignature.q025_trajectory ; consensusSignature.med_trajectory ; consensusSignature.q975_trajectory];

% transform & normalize trajectory
Y = exp(Y)./repmat(std(exp(Y(2,:))),3,length(itimes));

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

% get max of 95% CI, so can set Y-axis limits
maxL = max(Y(3,:))+0.5;

% plot vertical lines to show antibiotic treatment intervals
for intv=1:length(timeIntervals{refSubj}),
    ti = timeIntervals{refSubj}{intv};
    if ti.treat == 1,
        plot([times{refSubj}(ti.startTimeIdx) times{refSubj}(ti.startTimeIdx)],[0 maxL],'-b','LineWidth',2);
        plot([times{refSubj}(ti.endTimeIdx) times{refSubj}(ti.endTimeIdx)],[0 maxL],'-b','LineWidth',2);
    end;
end;

% plot median trajectory
plot(itimes,Y(2,:),'--r','LineWidth',4);

axis tight;

set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 20);
set(gca, 'LineWidth',3);

set(gca,'xtick',[],'ytick',[]);
