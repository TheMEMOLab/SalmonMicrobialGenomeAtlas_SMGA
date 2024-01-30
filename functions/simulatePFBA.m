function resPath = simulatePFBA(reconstructionFolder,resultFolder)
% This function performs parsimonious FBA to generate a flux vector for all
% reactions in the models given a given maximal growth constraint.
%
% USAGE: resPath = simulatePFBA(reconstructionFolder,resultFolder)
%
% INPUT
% reconstructionFolder    Folder with refined reconstructions
% resultFolder            Folder to save the pFBA results
%
% AUTHOR: Tim Hensen, 12/2023.

% function: Get pFBA results for all models 

% input
%refinedReconstructionFolder = [pwd filesep 'refinedReconstructionsMet'];
%resultFolder = [pwd filesep 'analysis'];

% Set solvers
changeCobraSolver('ibm_cplex','LP',-1);
changeCobraSolver('ibm_cplex','QP',-1);

% Get reconstruction paths
folder = what(reconstructionFolder);
modelPaths = string(append(folder.path, filesep, folder.mat));

% Get strain names
strainNames = string(erase(folder.mat,'_testmultiplegenomes'));
strainNames = string(erase(strainNames,'.mat'));

% Get complex medium
constraints = readtable('ComplexMedium.txt', 'Delimiter', 'tab');
constraints=table2cell(constraints);
constraints=cellstr(string(constraints));

pFBA = struct();

% Load model
warning('off','all')
for i=1:length(modelPaths)
    model = load(modelPaths(i));
    model = model.(string(fieldnames(model)));
    
    % apply complex medium
    warning('off','all')

    model = useDiet(model,constraints);

    % Initiate table
     resTable = table('Size',[length(model.rxns) 3],...
    'VariableTypes', {'string','string','double'},...
    'VariableNames',{'Reaction','Subsystem',char(strainNames(i))});

    % Add reaction info
    resTable.Reaction = model.rxns;
    resTable.Subsystem = string(model.subSystems);

    % Set model to operate in an anaerobe environment
    model = changeRxnBounds(model,'EX_o2(e)',0,'b');

    % Maximise for the biomass production
    model = changeObjective(model,'bio1');

    % Constrain biomass production
    model.osenseStr = 'max';
    FBA = optimizeCbModel(model);
    newBounds = FBA.f*0.75;
    model = changeRxnBounds(model,'bio1',newBounds,'b');

    % Perfrom pFBA
    model.osenseStr = 'min';
    QP = optimizeCbModel(model, 'min', 1e-6, 1);
    QP.v(abs(QP.v)<1e-6)=0;
    resTable.(strainNames(i)) = QP.v;

    % Save in structure
    pFBA.(strainNames(i)) = resTable;
end
warning('on','all')

resPath = [resultFolder filesep 'pFBA_results.mat'];
save(resPath,'pFBA')
end
