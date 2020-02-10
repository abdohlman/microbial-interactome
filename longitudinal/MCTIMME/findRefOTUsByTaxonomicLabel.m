function [otuList] = findRefOTUsByTaxonomicLabel(tlabel,tlevel,origOTUIdx,refOTU_taxonomy)
% helper function to find indices of refOTUs with a certain taxonomic label
%
% inputs:
% tlabel = text string containing taxonomic label, e.g., 'Bacteroidales')
% tlevel = integer for taxonomic level (order = 4, family = 5, genus = 6)
% origOTUIdx = indices of refOTUs so that we can map back to original #'s,
% calculated by filterData
% refOTU_taxonomy = cell array of taxonomic labels for each refOTU
%
% output:
% otuList = cell array of refOTUs that match label in each subject
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% findRefOTUsByTaxonomicLabel.m (version 1.00)

numSubjects = length(origOTUIdx);
taxonomy = parseTaxonomy(refOTU_taxonomy);
otuList = cell(numSubjects,1);

for ds=1:numSubjects,
    numOTUs = length(origOTUIdx{ds});
    ml = [];
    for otu=1:numOTUs,
        origOTU = origOTUIdx{ds}(otu);
        if strcmp(taxonomy{origOTU,tlevel},tlabel) == 1,
            ml = [ml otu];
        end;
    end;
    otuList{ds} = ml;
end;

function [taxonomy2] = parseTaxonomy(taxonomy)
% utility function to parse ';' delimited taxonomy into cell structure
taxonomy2 = cell(length(taxonomy),6);
for i=1:length(taxonomy),
    [s] = getLevels(taxonomy{i});
    for j=1:length(s),
        taxonomy2(i,j) = s(j);
    end;
end;

function [s2] = getLevels(tax)
% helper function for parse taxonomy
s =  regexp(tax, '[;_]', 'split');
s2 = cell(1,length(s));

filterC = [];
for i=1:length(s),
    s2{i} = s{i};
    [x, status] = str2num(s{i});
    unc = strmatch('uncultured',s{i});
    unc2 = strmatch('Incertae Sedis',s{i});
    if status==1 | length(unc)>0 | length(unc2)>0,,
        filterC = [filterC i];
    end;
end;
s2(filterC) = [];

