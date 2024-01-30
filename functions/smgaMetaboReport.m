%% % MetaboRePort workflow tutorial:
%cd('C:\Users\mspg\Documents\reconstruction')
% Set path to the cobratoolbox
global CBTDIR

% user defined path
folder = [pwd filesep 'translatedDraftReconstructions']; % Set path to folder with reconstructions

% Make folder to save the updated reconstructions
updatedReconstructPath = [pwd filesep 'updatedReconstructions'];
if exist(updatedReconstructPath,'dir')~=7
    mkdir(updatedReconstructPath)
end

Make folder to save the updated reconstructions as sbml files
annotatedSBMLreconstructions = [pwd filesep 'annotatedSBMLreconstructions'];
if exist(annotatedSBMLreconstructions,'dir')~=7
    mkdir(annotatedSBMLreconstructions)
end

% Make folder to where the reports are saved 
reportDir = [pwd filesep 'draftReports'];
if exist(reportDir,'dir')~=7
    mkdir(reportDir)
end

% Set metabolite structure path
metstructPath = [CBTDIR filesep 'tutorials' filesep 'dataIntegration' filesep...
 'metaboAnnotator' filesep 'data' filesep 'met_strc_rBioNet_new.mat'];

% Load rBioNet metabolite structure information
load(metstructPath);

% Get reconstruction filenames
[modelPaths, modelList] = getModelPaths(folder);

% Preallocate ScoresOverall table for speed
ScoresOverall = cell(length(modelList),2);
modelProperties = struct();
%%
tic;
for i = 1 : length(modelList)
    disp(i)
    % Load model 
    model = load(modelPaths(i));
    model = model.(string(fieldnames(model)));
        
    % Populate and further annotate model with metabolite info
    [modelUpdated] = populateModelMetStr(model, metabolite_structure_rBioNet,1);
    [modelUpdated] = annotateSBOTerms(modelUpdated);
    
    if any(contains(fieldnames(modelUpdated), {'metInChIString'}))
        modelUpdated = rmfield(modelUpdated,'metInChIString'); % wrongly in microbe models
    end
    
    % Add reaction information
    [modelUpdated] = populateModelwithRxnIDs(modelUpdated);
    
    % Get metaboscore
    [modelProp2,ScoresOverall2] = generateMetaboScore(model);
    
    % Store scores in modelProperties and ScoresOverall
    modelProperties.(modelList(i)).ScoresOverall = ScoresOverall2;
    modelProperties.(modelList(i)).modelUpdated = model;
    modelProperties.(modelList(i)).modelProp2 = modelProp2;
    
    ScoresOverall{i,1} = modelList(i);
    ScoresOverall{i,2} = num2str(ScoresOverall2);
    
    if mod(i,10)
        save('MetaboRePorts.mat','modelProperties','ScoresOverall');
    end
    
    % save updated mat file
    model = modelUpdated;
    save(strcat(updatedReconstructPath, filesep, modelList(i), '.mat'),'model');
    
    generate sbml file
    remove description from model structure as this causes issues
    if any(contains(fieldnames(modelUpdated), {'description'}))
        modelUpdated = rmfield(modelUpdated,'description');
    end
    sbmlPath = char(strcat(annotatedSBMLreconstructions, filesep, 'Annotated_',modelList(i)));
try 
    outmodel = writeCbModel(modelUpdated, 'format','sbml', 'fileName', sbmlPath);    
catch ME
    % Remove nested cells
    nestedCells = find(cellfun(@iscell,modelUpdated.subSystems));
    for k = 1:length(nestedCells)
        modelUpdated.subSystems{nestedCells(k)} = char(modelUpdated.subSystems{nestedCells(k)});
    end
    outmodel = writeCbModel(modelUpdated, 'format','sbml', 'fileName', sbmlPath);    
end

end
toc;

% Create MetaboReport
generateMetaboReport(modelProperties,reportDir);
save('MetaboRePorts.mat','modelProperties','ScoresOverall')
%%
% Save summary of scores
fields = string(fieldnames(modelProperties));
totalScores = cell(length(fields),6);
for i=1:length(fields)
    score = modelProperties.(fields(i)).modelProp2.Scores;
    reconScores = [score.Consistency,...
        score.AnnotationMetabolites,...
        score.AnnotationReactions,...
        score.AnnotationGenes,...
        score.AnnotationSBO,...
        score.Overall];
    totalScores(i,:) = num2cell(reconScores);
end
headers = {'Consistency', 'Metabolite','Reaction','Genes','SBO','Overall'};
totalScores = [headers; totalScores];
writecell(totalScores,'metaboReport_drafts.xlsx')




%%
function [paths, names] = getModelPaths(directory)
% This function finds the paths and names of COBRA models in a directory.
% Tim Hensen, November 2023

% Get reconstruction names
modelDir = dir(directory);
nameList = {modelDir.name}';

% Filter on .mat files
nameList = string(nameList(contains(nameList, '.mat')));

% Obtain paths
if matches(char(directory(end)),filesep) % Check if the path contains a fileseparator at the end
    paths = append(directory, nameList);
else
    paths = append(directory, filesep, nameList);
end

% Get reconstruction names
names = string(erase(nameList,'.mat'));
end