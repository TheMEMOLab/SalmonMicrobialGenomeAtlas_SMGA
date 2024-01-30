function [consumedMets, secretedMets, averages] = getConsumedSecretedMets(resultPath)
% This function obtains the metabolite exchange fluxes obtained using
% parsimonuous FBA and collects all consumed metabolites, secreted metabolites, and their fluxes.
% In additin, it calculates the average fluxes for each exchange reaction.
%
% USAGE: [consumedMets, secretedMets, averages] = getConsumedSecretedMets(resultPath)
%
% INPUT
% subsysResultPath      Path to the processed pFBA results from the
%                       saveAndProcessFluxes function
% resultFolder          Path to folder where the outputs will be saved
%
% OUTPUT
% consumedMets          Top 25 most consumed metabolites and their average fluxes
%                       per genus
% secretedMets          Top 25 most secreted metabolites and their average fluxes
%                       per genus
% secretedMets          Total flux averages for each investigated
%                       metabolite
% averages
% AUTHOR: Tim Hensen, 12/2023.

% Get fluxes
fluxes=readtable(resultPath); 

% Remove subsystem info
fluxes.Subsystem=[];

% Filter on exchange reactions
fluxes(~contains(fluxes.Reaction,'EX_'),:)=[];
fluxes(contains(fluxes.Reaction,'biomass'),:)=[];

% Add average exchange flux for all microbes
fluxVals = table2array(fluxes(:,3:end));
fluxes.Mean = mean(fluxVals,2,'omitnan');
fluxes = sortrows(fluxes,'Mean','descend');

% Get metabolite average
averages = table(fluxes.Reaction, fluxes.Mean,'VariableNames',{'Reaction','Average flux'});

% Remove Mean
%fluxes.Mean = [];

% Get microbe taxonomy info
taxonomy = readtable('InputData\adaptedInfoFile.txt');
taxonomy = taxonomy(:,{'MicrobeID','Genus'});

% Change microbe ID
taxonomy.MicrobeID = replace(taxonomy.MicrobeID,' ','_');
taxonomy.MicrobeID = replace(taxonomy.MicrobeID,'.','');

% Get VMH databases
database= loadVMHDatabase;

% Filter on top 25 and bottom 25
top = fluxes(1:25,:);
bottom = fluxes(end-24:end,:);

% Obtain summarised fluxes per genus for most consumed and secreted
% metabolites
secretedMets = processExchangeFluxes(top,taxonomy, database);
consumedMets = processExchangeFluxes(bottom,taxonomy, database);
end

% Create function to process top and bottom fluxes
function genusFluxTable = processExchangeFluxes(fluxes,taxonomy, database)
% transpose fluxes 
fluxes_transposed = rows2vars(fluxes,'VariableNamesSource', 'Reaction','VariableNamingRule','preserve');
fluxes_transposed.Properties.VariableNames(1) = {'MicrobeID'};

% join tables
fluxes_Genus=outerjoin(fluxes_transposed,taxonomy,'Keys','MicrobeID','MergeKeys',true);
fluxes_Genus.Genus(matches(fluxes_Genus.MicrobeID,'Mean')) = {'Overall mean'};

% Find genus groups
[genus_groups,genus] = findgroups(fluxes_Genus.Genus);

% Get reaction names
rxnNames=string(fluxes_Genus.Properties.VariableNames(2:end-1));

% Preallocate table
genusFluxes = nan(length(rxnNames),length(genus));
for i=1:height(genusFluxes)
    meanf = @(x) mean(x,'omitnan');
    genusFluxes(i,:) =  splitapply(meanf, fluxes_Genus.(rxnNames(i)),genus_groups);
end

% Create table
genusFluxTable = array2table(genusFluxes,'VariableNames',genus,'RowNames',rxnNames);

% Replace reactions with VMH metabolite names
VMHID = extractBetween(genusFluxTable.Properties.RowNames,'EX_','(e)');

% Find metabolites
[~,~,ib]=intersect(VMHID,database.metabolites(:,1),'stable');
genusFluxTable.Properties.RowNames = database.metabolites(ib,2);

% move mean to last spot
genusFluxTable = movevars(genusFluxTable,'Overall mean','after','Weissella');
end