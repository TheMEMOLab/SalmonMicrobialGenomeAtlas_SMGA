function [uptakeTable, secretionProductTable] = extendInputTables(metabolomeRxnPath, inputDataFolder)
% This function extends the input tables in DEMETER with the metabolites to
% be added in the reconstructions.
%
% USAGE: [uptakeTable, secretionProductTable] = extendInputTables(metabolomeRxnPath, inputDataFolder)
%
% INPUT
% metabolomeRxnPath     Path to .mat file with the metabolites to be added 
%                       (additionalRxns.mat)
% inputDataFolder       Path to folder with all the input data for DEMETER  
%
% OUTPUT
% uptakeTable           Table with all metabolites to be added by uptake
%                       exchange reactions.
% secretionProductTable Table with all metabolites to be added by secretion
%                       exchange reactions.
% AUTHOR: Tim Hensen, 12/2023.

% Load metabolites
additionalRxns = load(metabolomeRxnPath);
additionalRxns = additionalRxns.(string(fieldnames(additionalRxns)));

% Get microbes
microbes = string(fieldnames(additionalRxns));

% Get tables
uptakeTable=readInputTableForPipeline([inputDataFolder filesep 'uptakeTable.txt']);
uptakeTable(:,strncmp(uptakeTable(1,:),'Ref',3))=[];

secretionProductTable=readInputTableForPipeline([inputDataFolder filesep 'secretionProductTable.txt']);
secretionProductTable(:,strncmp(secretionProductTable(1,:),'Ref',3))=[];

for i=1:length(microbes)
    % Get consumed and secreted metabolites
    metabolites = additionalRxns.(microbes(i));
    
    % Extend uptake table
    exchanged = metabolites(metabolites.Consumed==1,:);
    metsToAdd = exchanged.VMH_metabolite;
    microbeUptakeIDX = find(matches(uptakeTable(:, 1),microbes(i)));
    uptakeTable = extendSecretionUptakeTable(uptakeTable, metsToAdd, microbeUptakeIDX);
    
    % Extend secretion table
    exchanged = metabolites(metabolites.Secreted==1,:);
    metsToAdd = exchanged.VMH_metabolite;
    microbeSecretionIDX = find(matches(secretionProductTable(:, 1),microbes(i)));
    secretionProductTable = extendSecretionUptakeTable(secretionProductTable, metsToAdd, microbeSecretionIDX);
end
uptakeTable(:,end+1) = [{'Ref'}; repmat({''},height(uptakeTable)-1,1)];
secretionProductTable(:,end+1) = [{'Ref'}; repmat({''},height(secretionProductTable)-1,1)];


% Save tables
writetable(cell2table(uptakeTable),[inputDataFolder filesep 'uptakeTable'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
writetable(cell2table(secretionProductTable),[inputDataFolder filesep 'secretionProductTable'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
end


function data = extendSecretionUptakeTable(data, metsToAdd, microbeIndex)
% Check if metabolites are present in uptakeTable
dataMets = data(1,2:end);

% Add metabolites to table 
data(microbeIndex,matches(dataMets,metsToAdd)) = {1};

% if there are metabolites not yet in the table, they will be added
if any(~matches(metsToAdd,dataMets))
    
    % Extend cell array
    metsToExtend = metsToAdd(~matches(metsToAdd,dataMets));
    dataExtension = num2cell(zeros(height(data)-1,length(metsToExtend)));
    dataExtension = [metsToExtend' ; dataExtension];

    % Add metabolites
    dataExtension(microbeIndex,:) = {1};

    % Add extension to data
    data = [data dataExtension];
end
end