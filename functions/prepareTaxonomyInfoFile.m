function adaptedInfoFilePath = prepareTaxonomyInfoFile(taxpath)
% This function prepares the microbial taxonomy data for further processing
% by the prepareInputData function
%
% USAGE: adaptedInfoFilePath = prepareTaxonomyInfoFile(taxpath)
%
% INPUT
% taxpath               Path to taxonomies
%
% OUTPUT
% adaptedInfoFilePath   Path to adapted taxonomy info files
%
% AUTHOR: Tim Hensen, 12/2023.

% Load taxonomy path
taxonomyfile = readtable(taxpath);

% Change microbeID
taxonomyfile.GenomeID = replace(taxonomyfile.GenomeID,' ','_');
taxonomyfile.GenomeID = replace(taxonomyfile.GenomeID,'.','');


% Load AGORA2 input file path
global CBTDIR
AGORA2_infofolder = [CBTDIR filesep 'papers' filesep '2021_demeter' filesep 'input' filesep 'AGORA2_infoFile.xlsx'];
AGORA2_info = readtable(AGORA2_infofolder);

% Add gram status
idx = matches(taxonomyfile.Phylum,'Firmicutes') | matches(taxonomyfile.Phylum,'Actinobacteria');
taxonomyfile.GramStaining(idx) = "Gram+";
taxonomyfile.GramStaining(~idx) = "Gram-";

% Obtain AGORA2 presence
idx = matches(taxonomyfile.GenomeID,AGORA2_info.Strain);
taxonomyfile.PresentInAGORA2(idx) = "Yes";
taxonomyfile.PresentInAGORA2(~idx) = "No";

% Create table for demeter
inputTable=struct();
inputTable.MicrobeID = taxonomyfile.GenomeID;
inputTable.PubSeedID = repmat('N/A',height(taxonomyfile),1);
inputTable.("x2_CompleteComparativeGenomics_1_CertainComparativeGeno") = nan(height(taxonomyfile),1);
inputTable.Strain = taxonomyfile.GenomeID;
inputTable.Species = taxonomyfile.Species;
inputTable.Genus = taxonomyfile.Genus;
inputTable.Family = nan(height(taxonomyfile),1);
inputTable.Order = taxonomyfile.Order;
inputTable.Class = nan(height(taxonomyfile),1);
inputTable.Phylum = taxonomyfile.Phylum;
inputTable.NCBItaxID = nan(height(taxonomyfile),1);
inputTable.AnnotationVersionID = nan(height(taxonomyfile),1);
inputTable.Gram_Staining = taxonomyfile.GramStaining;
inputTable.("PresentInAGORA2") = taxonomyfile.PresentInAGORA2;

inputTable=struct2table(inputTable);

inputTable.('Oxygen Requirement') = repmat("",length(taxonomyfile.GenomeID),1);


adaptedInfoFilePath = 'InputData\adaptedInfoFile.txt';
writetable(inputTable,'InputData\adaptedInfoFile.xlsx','WriteVariableNames',true);
writetable(inputTable,adaptedInfoFilePath,'WriteVariableNames',true,'Delimiter','tab');
end

