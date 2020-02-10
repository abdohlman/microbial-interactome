function [D_raw,times,timeIntervals,refOTU_taxonomy] = loadDethlefsenData(dataDir)
% load data for the three subjects from Dethlefsen et al. PNAS 2010 dataset
% 
% input:
% dataDir = directory w/ data files
%
% outputs:
% D_raw = cell array w/ refOTU counts for each subject
% times = cell array w/ times for each subject
% timeIntervals = cell array w/ timeIntervalObjects for each subject
% refOTU_taxonomy = array w/ taxonomy labels for each refOTU
%
% each subject has three data files:
% 1) subjX.txt - counts for refOTUs, tab-delimited; rows=refOTUs, cols=time-points
% 2) subjX_times.txt - days (starting w/ day 1) for each time-point
% 3) subjX_timeIntervals.txt - specifies treatment/non-treatment intervals
%   (see readTimeIntervalFile.m for details)
%
% refOTU_taxonomy.txt contains the SILVA taxonomy names for refOTUs,
% which are shared across all three subjects
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% loadDethlefsenData.m (version 1.00)

% create cell array to store refOTU counts
D_raw = cell(3,1);

% load refOTU counts
D_raw{1} = load([dataDir 'subjD_data.txt']);
D_raw{2} = load([dataDir 'subjE_data.txt']);
D_raw{3} = load([dataDir 'subjF_data.txt']);

% create cell array to store times
times = cell(3,1);

% load times
times{1} = load([dataDir 'subjD_times.txt']);
times{2} = load([dataDir 'subjE_times.txt']);
times{3} = load([dataDir 'subjF_times.txt']);

% load time intervals
timeIntervals = cell(3,1);
timeIntervals{1} = readTimeIntervalFile([dataDir 'subjD_timeIntervals.txt'],times{1});
timeIntervals{2} = readTimeIntervalFile([dataDir 'subjE_timeIntervals.txt'],times{2});
timeIntervals{3} = readTimeIntervalFile([dataDir 'subjF_timeIntervals.txt'],times{3});

% load refOTU taxonomy
fid = fopen([dataDir 'refOTU_taxonomy.txt'],'r');
refOTU_taxonomy = cell(size(D_raw{1},1),1);
i = 1;
tline = fgetl(fid);
while ischar(tline),
    refOTU_taxonomy{i} = tline;
    tline = fgetl(fid);
    i = i+1;
end;
fclose(fid);

