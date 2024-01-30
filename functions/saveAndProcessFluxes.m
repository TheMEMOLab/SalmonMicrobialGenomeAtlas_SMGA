function resultPath = saveAndProcessFluxes(pFBA_path,additionalTestPath)
% This function loads the pFBA results and puts them in a single table with
% rows for the union of all reactions in the models.
%
% USAGE: resultPath = saveAndProcessFluxes(pFBA_path,additionalTestPath)
%
% INPUT
% pFBA_path             Path to pFBA results
% additionalTestPath    Path to folder where the outputs will be saved
%
% OUTPUT
% resultPath            Path to file containing the function outputs
%
% AUTHOR: Tim Hensen, 12/2023.

% Load pFBA results
pFBA = load(pFBA_path);
pFBA = pFBA.(string(fieldnames(pFBA)));

% Step 1: Iteratively find all unique reactions 

% Get fieldnames
strains = string(fieldnames(pFBA));

% Get all unique reactions and subsystems
totalRxns = repmat("",1,2);
for i=1:length(strains)
    modelRes = pFBA.(strains(i));
    totalRxns = [totalRxns ; [modelRes.Reaction modelRes.Subsystem]];
    [~,idx] = unique(totalRxns(:,1));
    totalRxns = totalRxns(idx,:);
end
totalRxns(1,:)=[];

% Step 2: Create table for all results
rxnTable = table(totalRxns(:,1),totalRxns(:,2),'VariableNames',{'Reaction','Subsystem'});
strainTable = array2table(nan(height(totalRxns),length(strains)),'VariableNames',strains);
pFBAtable = [rxnTable strainTable];

% Step 3: Iteratively populate result table
for i=1:length(strains)
    % Get strain
    modelRes = pFBA.(strains(i));
    
    % Find indices of reactions in strain
    idRxn = matches(pFBAtable.Reaction,modelRes.Reaction);
    
    % Add results to indices in the strains
    pFBAtable.(strains(i))(idRxn) = modelRes.(strains(i));
end

% Remove rows without flux in none of the strains
checkFluxes = table2array(pFBAtable(:,3:end));
pFBAtable(all(isnan(checkFluxes)| checkFluxes==0,2),:)=[];

% Save fluxes
resultPath = [additionalTestPath filesep 'pFBA_fluxes.xlsx'];
writetable(pFBAtable,resultPath)
end

