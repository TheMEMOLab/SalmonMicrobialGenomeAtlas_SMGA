function [resPath,stats] = testModelsAgainstMetabolomics(modelPath, metabolomicsPath,additionalTestPath)
% This functions performs qualitatively tests if the models can consume metabolites consumed in the exometabolomics
% data and secreted those secreted in the exometabolomics using flux variability analysis. 
%
% USAGE: [resPath,stats] = testModelsAgainstMetabolomics(modelPath, metabolomicsPath,additionalTestPath)
%
% INPUT
% modelPath                 Path to refined and gapfilled models
% metabolomicsPath          Path to metabolites to investigate
% ('additionalRxns.mat')
%
% OUTPUT
% resPath                   Path to test results    
% stats                     Output in matlab
%
% AUTHOR: Tim Hensen, 12/2023.

% Test models using FVA and compare against metabolomics

% Set resultPath
resPath = [additionalTestPath filesep 'FVA_metabolomicsCheck.xlsx'];

% Remove previous file
delete(resPath)
% Get reconstruction paths
folder = what(modelPath);
paths = string(append(folder.path, filesep, folder.mat));

% Get strain names
strains = string(erase(folder.mat,'.mat'));

% Get microbe names
microbes = extractBefore(strains,'_');

% check which reactions confirm with the metabolomics
additionalRxns = load(metabolomicsPath);
additionalRxns = additionalRxns.(string(fieldnames(additionalRxns)));

% Preallocate statistics table
stats = table('Size',[3,4], 'VariableTypes',{'string','double','double','double'},...
    'VariableNames',{'Strain','TP','FN','Sensitivity'});
stats.Strain = strains;

for i = 1:length(strains)

    % load models
    model = load(paths(i));
    model = model.(string(fieldnames(model)));
    
    % Get all exchanged metabolites
    metabolomicsExMets=additionalRxns.(strains(i));
    
    % Remove metabolites without exchanges
    metabolomicsExMets(metabolomicsExMets.Consumed == 0 & metabolomicsExMets.Secreted == 0,:)=[];
    
    % Remove metabolites from testing that are probably wrong
    metabolomicsExMets(contains(metabolomicsExMets.VMHID,{'3mhis','allop','M03062','n8aspmd'}),:)=[];
    metabolomicsExMets(contains(metabolomicsExMets.VMH_metabolite,{'Acetyl','N-acetyl'}),:)=[];
    
    % Set diet
    
    % Get dietary reactions
    idx = metabolomicsExMets.Consumed==1;
    DietRxns = cellstr(append('EX_',metabolomicsExMets.VMHID(idx),'(e)'));
    
    % Close all consumptions reactions
    model = changeRxnBounds(model, model.rxns(contains(model.rxns,'EX_')), 1000, 'u');
    model = changeRxnBounds(model, model.rxns(contains(model.rxns,'EX_')), 0, 'l');
    
    % Open diet reactions
    model = changeRxnBounds(model, model.rxns(matches(model.rxns,DietRxns)), -1000, 'l');
    model = changeRxnBounds(model, model.rxns(matches(model.rxns,DietRxns)), 0, 'u');
    
    % Test reactions
    rxnsToCheck = cellstr(append('EX_',metabolomicsExMets.VMHID,'(e)'));
    [minFlux, maxFlux] = guidedSim(model, rxnsToCheck);
    
    % Obtain correctly predicted compounds
    fvares = table(rxnsToCheck, metabolomicsExMets.VMHID, minFlux, maxFlux,'VariableNames',{'Rxn','VMHID','minFlux','maxFlux'});
    fvares = innerjoin(fvares,metabolomicsExMets);
    fvares.Correct = (fvares.minFlux<0 & fvares.Consumed==1) | (fvares.maxFlux>0 & fvares.Secreted==1);
    
    % Obtain summary statistics
    stats.TP(i) = sum(fvares.Correct);
    stats.FN(i) = sum(~fvares.Correct);
    stats.Sensitivity(i) = stats.TP(i)/(stats.TP(i)+stats.FN(i));
    
    % Annotate the fva results for consumption
    fvares.Consumption_in(fvares.Consumed & fvares.minFlux>=0) = {'metabolome'};
    fvares.Consumption_in(fvares.Consumed & fvares.minFlux<0) = {'metabolome and model'};
    fvares.Consumption_in(~fvares.Consumed & fvares.minFlux<0) = {'model'};
    fvares.Consumption_in(cellfun(@isempty,fvares.Consumption_in)) = {'not consumed'};
    
    % Annotate the fva results for secretion
    fvares.Secretion_in(fvares.Secreted & fvares.maxFlux<=0) = {'metabolome'};
    fvares.Secretion_in(fvares.Secreted & fvares.maxFlux>0) = {'metabolome and model'};
    fvares.Secretion_in(~fvares.Secreted & fvares.maxFlux>0) = {'model'};
    fvares.Secretion_in(cellfun(@isempty,fvares.Secretion_in)) = {'not secreted'};
    
    % Save test results
    writetable(fvares,resPath,"Sheet",microbes(i));
end

% Save statistics
writetable(stats,resPath,"Sheet",'Statistics');
end
