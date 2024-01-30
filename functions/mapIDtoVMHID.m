function [path, mappingData]= mapIDtoVMHID(mappingPath)
% This function maps metabolite names, HMDB IDs, InchiKey IDs, or SMILES to
% VMH IDs. 
%
% USAGE: [path, mappingData]= mapIDtoVMHID(mappingPath)
%
% INPUT
% mappingPath            Directory to file with metabolite names to map
%                        
%
% OUTPUT
% path                   Path to file with mapping results
% mappingData            Mapped metabolites
%
% AUTHOR: Tim Hensen, 12/2023.

% Load mediaData
metAnnotations = readtable(mappingPath,'ReadVariableNames',false);
metAnnotations.Properties.VariableNames = {'metabolite','alternative name','HMDB','InchiKey','SMILES'};

% Map metabolomics compounds onto VMH IDs
targetIDtypes = {'HMDB','InchiKey','SMILES'};
sourceIDtypes = {'HMDB','InchiKey','SMILES'};

% Load mapping database
metaboliteDatabase=readtable('metabolon_crossmatch_IT_withUpdatedInchiKey.xlsx','PreserveVariableNames',1);

% Remove duplicate sugar c6 metabolite
metAnnotations(115,:) = [];

% Add VMH metabolite names
VMH=loadVMHDatabase;
metabolites = VMH.metabolites;

% Create new dataset
mappingData = metAnnotations;
for i = 1:length(targetIDtypes)  
    % input
    targetIDs = metAnnotations.(targetIDtypes{i});
    % Set empty cells as NA
    targetIDs(cellfun(@isempty,targetIDs)) = {'N/A'};
    IDtype = sourceIDtypes{i};

    % Get IDtype
    IDtype = upper(IDtype);
    sourceIDs = metaboliteDatabase.(IDtype);

    % Map target IDs against source IDs
    [~,ia,ib] = intersect(targetIDs,sourceIDs,'stable');

    % Obtain ID mapping
    mappingData.VMHID(ia) = metaboliteDatabase.('VMH ID')(ib);

    % Map target against VMH database
    switch targetIDtypes{i}
        case 'InchiKey'
            [~,ia,ib] = intersect(targetIDs,metabolites(:,9),'stable');
        case 'SMILES'
            [~,ia,ib] = intersect(targetIDs,metabolites(:,10),'stable');
        case 'HMDB'
            [~,ia,ib] = intersect(targetIDs,metabolites(:,11),'stable');
    end
    % Add mapped metabolites
    mappingData.VMHID(ia) = metabolites(ib,1);
end

% Make empty cells into strings
id=cellfun(@isempty,mappingData.VMHID);
mappingData.VMHID(id) = {''};
mappingData.VMHID = string(mappingData.VMHID);

%% Manually add VMH metabolites
%mappingData.metabolite(mappingData.VMHID=="");
mappingData.VMHID(matches(mappingData.metabolite,{'4-Hydroxybenzaldehyde'})) = "4hbald";
%mappingData.VMHID(matches(mappingData.metabolite,{'7-Methylguanine'})) =
%"M03062"; % mismapping
mappingData.VMHID(matches(mappingData.metabolite,{'AMP'})) = "3amp";
mappingData.VMHID(matches(mappingData.metabolite,{'Lysine'})) = "lys_L";
mappingData.VMHID(matches(mappingData.metabolite,{'Methionine sulfoxide'})) = "metsox_S_L";
mappingData.VMHID(matches(mappingData.metabolite,{'Nicotinic acid'})) = "nac";
mappingData.VMHID(matches(mappingData.metabolite,{'O-Acetylserine'})) = "acser";
mappingData.VMHID(matches(mappingData.metabolite,{'Threonine'})) = "thr_L";
mappingData.VMHID(matches(mappingData.metabolite,{'Urocanate'})) = "urcan";
mappingData.VMHID(matches(mappingData.metabolite,{'? CMP'})) = "3cmp";
mappingData.VMHID(matches(mappingData.metabolite,{'4-Coumaric acid'})) = "T4hcinnm";
mappingData.VMHID(matches(mappingData.metabolite,{'5-Methylcytosine'})) = "5mtc";
mappingData.VMHID(matches(mappingData.metabolite,{'4-Coumaric acid'})) = "cbasp";
mappingData.VMHID(matches(mappingData.metabolite,{'Diaminopimelic acid'})) = "26dap_LL";
mappingData.VMHID(matches(mappingData.metabolite,{'Glucosamine'})) = "gam";
mappingData.VMHID(matches(mappingData.metabolite,{'Glutathione oxidized'})) = "gthox";
mappingData.VMHID(matches(mappingData.metabolite,{'Malic acid'})) = "mal_L";
mappingData.VMHID(matches(mappingData.metabolite,{'N8-Acetylspermidine'})) = "n8aspmd";
mappingData.VMHID(matches(mappingData.metabolite,{'Norepinephrine'})) = "nrpphr";
mappingData.VMHID(matches(mappingData.metabolite,{'Putrescine'})) = "ptrc";
mappingData.VMHID(matches(mappingData.metabolite,{'Pyridoxamine'})) = "pydam";
mappingData.VMHID(matches(mappingData.metabolite,{'Sugar (C6)'})) = "glc_D";
mappingData.VMHID(matches(mappingData.metabolite,{'UMP'})) = "3ump";
%%
% Map against VMH database to obtain metabolite names
[~,ia,ib] = intersect(mappingData.VMHID,metabolites(:,1),'stable');
mappingData.VMH_metabolite(ia) = metabolites(ib,2);

% Make empty cells into strings
id=cellfun(@isempty,mappingData.VMH_metabolite);
mappingData.VMH_metabolite(id) = {''};
mappingData.VMH_metabolite = string(mappingData.VMH_metabolite);

% Add mapping info
mappingData.Mapped = ~matches(mappingData.VMHID,"");

% Save mapping
path = [pwd filesep 'metaboliteMapping.xlsx'];
writetable(mappingData, path,'WriteVariableNames',true);
end