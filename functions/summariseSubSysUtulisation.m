function resultPath = summariseSubSysUtulisation(subsysResultPath, resultFolder)
% This function summarises the relative flux abundances per subsystem
% calculated in the getRelativeSubsystemFluxes function and calculates the
% mean relative flux abundances per genus.
%
% USAGE: resultPath = summariseSubSysUtulisation(subsysResultPath, resultFolder)
%
% INPUT
% subsysResultPath      Path to the output file from the
%                       getRelativeSubsystemFluxes function
% resultFolder          Path to folder where the outputs will be saved
%
% OUTPUT
% resultPath            Path to file containing the function outputs
%
% AUTHOR: Tim Hensen, 12/2023.

% Concatinate the subsystemAbundances by genus

% Load subsystem abundances
subUtul = readtable(subsysResultPath);

% Remove duplicate subsystem 'Wood-Ljungdahl Pathway'
subUtul(end,:) = [];


% Transpose results
subUtul_transposed = rows2vars(subUtul,'VariableNamesSource', 'Subsystem','VariableNamingRule','preserve');
subUtul_transposed.Properties.VariableNames(1) = {'MicrobeID'};

% Change microbeIDs to align with taxonomy info
subUtul_transposed.MicrobeID = replace(subUtul_transposed.MicrobeID,'.','');

% Get microbe taxonomy info
taxonomy = readtable('InputData\adaptedInfoFile.txt');
taxonomy = taxonomy(:,{'MicrobeID','Genus'});

% Change microbe ID
taxonomy.MicrobeID = replace(taxonomy.MicrobeID,' ','_');
taxonomy.MicrobeID = replace(taxonomy.MicrobeID,'.','');

% join tables
subUtul_Genus=outerjoin(subUtul_transposed,taxonomy,'Keys','MicrobeID','MergeKeys',true);

% Find genus groups
[genus_groups,genus] = findgroups(subUtul_Genus.Genus);

% Get subsystem names
subsysNames=string(subUtul_Genus.Properties.VariableNames(2:end-1));

% Preallocate table
genusSubsys = nan(length(subsysNames),length(genus));
for i=1:height(genusSubsys)
meanf = @(x) mean(x,'omitnan');
genusSubsys(i,:) =  splitapply(meanf, subUtul_Genus.(subsysNames(i)),genus_groups);
end

% Convert nan to zero
genusSubsys(isnan(genusSubsys))=0;

% Create table
genusSubsysTable = array2table(genusSubsys,'VariableNames',genus,'RowNames',subsysNames);

% Sort average subsystem utulisation over all genera
[~,idx] = sort(mean(genusSubsys,2),'descend');
genusSubsysTable = genusSubsysTable(idx,:);

% Add rownames as first column
genusSubsysTable.Subsystem = genusSubsysTable.Properties.RowNames;
% 
genusSubsysTable = movevars(genusSubsysTable,'Subsystem','Before',genusSubsysTable.Properties.VariableNames(1));

% Save all results
resultPath = [resultFolder filesep 'genus_RelativeSubsystemUtulisation.xlsx'];
writetable(genusSubsysTable,resultPath)
end



