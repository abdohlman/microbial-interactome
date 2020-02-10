function [CSGs] = consensusSignatureGroups(hparams,origOTUIdx,baseDir)
% find consensus signature groups (CSGs)
%
% inputs:
% hparams = hyperparameters object
% origOTUIdx = indices for original refOTUs, calculated by filterData
% baseDir = directory containing MCMC sample output
%
% outputs:
% CSGs = cell array, CSGs{i} is a cell array of refOTUObjects belonging to
% consensus signature group i
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% consensusSignatureGroups.m (version 1.00)

numSubjects = length(hparams.gamma);

% calculate matrix of co-clustering probabilities, where the (i,j)th entry
% for coClusterM is the frequency w/ which refOTU i & j are assigned to the
% same signature (note that the same refOTU in different subjects is
% separate entry in the matrix)
% the expected # of CSGs is the median of the # of signatures across all
% samples
[coClusterM,expectedNumCSGs] = coClusterProbs(baseDir,numSubjects);
% perform agglomerative clustering to obtain CSG
[Cl] = agglomCluster(1-coClusterM,expectedNumCSGs);
% create cell array of CSG members
[CSGs] = convertCSGs(Cl,origOTUIdx);

function [coClusterM,expectedNumCSGs,ofs] = coClusterProbs(baseDir,numSubjects)
% calculate co-clustering probabilities for refOTUs

% load assignments of refOTU time-series to signatures
signatureAssignsFileN = cell(numSubjects,1);
signatureAssigns = cell(numSubjects,1);
numOTUs = zeros(numSubjects,1);
ofs = zeros(numSubjects,1);
for ds=1:numSubjects,
    signatureAssignsFileN{ds} = [baseDir '_signatureAssigns' int2str(ds) '.txt'];
    signatureAssigns{ds} = dlmread(signatureAssignsFileN{ds});
    numOTUs(ds) = size(signatureAssigns{ds},2);
    if ds > 1,
        ofs(ds) = ofs(ds-1) + numOTUs(ds-1);
    end;
end;
% # of MCMC samples
numSamples = size(signatureAssigns{1},1);
% initialize (sparse) matrix for computing pairwise co-clustering
% probabilities between refOTUs
coClusterM = sparse(sum(numOTUs),sum(numOTUs));
expectedNumCSGs = zeros(1,numSamples);

for s=1:numSamples,
    numSignatures = 0;
    % determine # of signatures across all subjects
    for ds=1:numSubjects,
        numSignatures = max([numSignatures signatureAssigns{ds}(s,:)]);
    end;
    expectedNumCSGs(s) = numSignatures;
    for k=1:numSignatures,
        f = [];
        for ds=1:numSubjects,
            % find all refOTUs assigned to signature k
            ft = find(signatureAssigns{ds}(s,:) == k);
            if ~isempty(ft),
                ft = ft + ofs(ds);
                f = [f ft];
            end;
        end;
        % increment matrix entry for refOTUs sharing the same signature in
        % sample s
        for j=1:length(f),
            ff = find(f > f(j));
            if ~isempty(ff),
                coClusterM(f(j),f(ff)) =  coClusterM(f(j),f(ff)) + 1;
            end;
        end;
    end;
end;
coClusterM = coClusterM/numSamples;

expectedNumCSGs = median(expectedNumCSGs);

function [Cl] = agglomCluster(D,expectedNumClusters)
% agglomerative clustering using average linkage
% input:
% D = matrix of pairwise distances between data-points
% expectedNumClusters = final # of clusters desired

n = size(D,1);
Cl = cell(n,1);
for i=1:n,
    Cl{i} = i;
end;

merge1 = 0;
merge2 = 0;
d = 0.0;
minD = 100000;

so = sprintf('expected # CSGs = %i',expectedNumClusters);
disp(so);

% merge clusters until the expected # of clusters is reached
while merge1 >= 0 && length(Cl) > expectedNumClusters,
    merge1 = -1;
    merge2 = -1;
    minD = 10000;
    nc = length(Cl);
    so = sprintf('Merge level %i',nc);
    disp(so);
    for i=1:(nc-1),
        for j=(i+1):nc,
            d = avgLinkage(D,i,j,Cl);
            if d < minD,
                merge1 = i;
                merge2 = j;
                minD = d;
            end;
        end;
    end;
    if merge1 >=0,
        Cl{merge1} = union(Cl{merge1},Cl{merge2});
        Cl(merge2) = [];
    end;
end;

function [avgD] = avgLinkage(D,cluster1,cluster2,Cl)
% compute average linkage between cluster1 and cluster2
% input:
% D = distance matrix for data-points
% cluster1 = index of cluster1
% cluster2 = index of cluster2
% Cl = cell array specifying clusters, where Cl{i} is an array of indices
% of data-points belonging to cluster i

avgD = 0.0;
nc = 0;
d = 0.0;
n1 = length(Cl{cluster1});
n2 = length(Cl{cluster2});

for i=1:n1,
    if length(Cl{cluster1}) == 1,
        id1 = Cl{cluster1};
    else,
        id1 = Cl{cluster1}(i);
    end;
    for j=1:n2,
        if length(Cl{cluster2}) == 1,
            id2 = Cl{cluster2};
        else,
            id2 = Cl{cluster2}(j);
        end;
        if id2 > id1,
            d = D(id1,id2);
        else
            d = D(id2,id1);
        end;
        avgD = avgD + d;
        nc = nc + 1;
    end;
end;
avgD = avgD/nc;

function [CSGs] = convertCSGs(Cl,origOTUIdx)
% utility function to create cell array of refOTUObjects
% Cl = cell array, such that Cl{i} is an array specifying members of consensus signature group i
% hparams = hyperparameters object
% origOTUIdx = indices used for mapping back to original refOTUs, generated
% by filterData
% output:
% CSGs = cell array, such that Cl{i} is an array specifying members of
% consensus signature group i as refOTUObjects
numSubjects = length(origOTUIdx);
numOTUs = zeros(numSubjects,1);
ofs = zeros(numSubjects,1);
bdidx = zeros(numSubjects,2);
for ds=1:numSubjects,
    numOTUs(ds) = length(origOTUIdx{ds});
    if ds > 1,
        ofs(ds) = ofs(ds-1) + numOTUs(ds-1);
        bdidx(ds,1) = 1+ofs(ds);
        bdidx(ds,2) = bdidx(ds,1) + numOTUs(ds)-1;
    else,
        bdidx(ds,1) = 1;
        bdidx(ds,2) = numOTUs(ds);
    end;
end;

CSGs = cell(length(Cl),1);

for ncc=1:length(Cl),
    numItems = length(Cl{ncc});
    nri = 1;
    CSGMembers = cell(numItems,1);
    for ds=1:numSubjects,
        % find refOTUs from each subject that belong to this CSG
        f = find(Cl{ncc} >= bdidx(ds,1) & Cl{ncc} <= bdidx(ds,2));
        for onn=1:length(f),
            otu = Cl{ncc}(f(onn)) - ofs(ds);
            roo = refOTUObject;
            roo.subjNum = ds;
            roo.OTUNum = otu;
            roo.origOTU = origOTUIdx{ds}(otu);
            CSGMembers{nri} = roo;
            nri = nri+1;
        end;
    end;
    CSGs{ncc} = CSGMembers;
end;