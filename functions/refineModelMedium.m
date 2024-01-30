function [addedReactions,BlockedReactions, modelsRefined] = refineModelMedium(reconstructionFolder, metabolomicsPath)
% The correct exchange and transport reactions defined from the exometabolomics are not always correctly
% added in DEMETER. This function adds these missing reactions by finding
% which metabolites in the in vivo medium are blocked and then adding
% transport and exchange reactions for them. 
%
% USAGE: [addedReactions,BlockedReactions, modelsRefined] = refineModelMedium(reconstructionFolder, metabolomicsPath)
%
% INPUT
% reconstructionFolder  Folder with refined reconstructions
% metabolomicsPath      Path to the additionalRxns.mat file 
%
% OUTPUT
% addedReactions    List of reactions added to each model
% BlockedReactions  List of reactions that transport or exchange
%                   metabolites that are exchanged in vivo, but cannot be exchanged in the
%                   models.
% modelsRefined     Updated models
%
% AUTHOR: Tim Hensen, 12/2023.


% Get reconstruction paths
folder = reconstructionFolder;
paths = string(append(folder.path, filesep, folder.mat));

% Get strain names
strains = string(erase(folder.mat,'.mat'));

% check which reactions confirm with the metabolomics
additionalRxns = load(metabolomicsPath);
additionalRxns = additionalRxns.(string(fieldnames(additionalRxns)));

addedReactions = struct();
BlockedReactions = struct();
modelsRefined = struct();

for i = 1:length(paths)
    % load models
    model = load(paths(i));
    model = model.(string(fieldnames(model)));

    % Get all exchanged metabolites
    metabolomicsExMets=additionalRxns.(strains(i));

    % Remove uridine duplicate
    metabolomicsExMets = sortrows(metabolomicsExMets,'VMHID');
    metabolomicsExMets(90,:)=[];

    % Get reactions
    metabolomicsExMets(metabolomicsExMets.Consumed==0 & metabolomicsExMets.Secreted==0,:)=[];

    %% Add exchange reactions to conform with metabolomics
    rxn = append('EX_',metabolomicsExMets.VMHID,'(e)');
    rxn_formula = append(metabolomicsExMets.VMHID, '[e]  <=> ');
    for j=1:length(rxn)
        [model,rxnIDexists] = addReaction(model,rxn{j}, rxn_formula{j});    

        % Get added reactions
        exchangeAdded = setdiff(rxn{j},model.rxns(rxnIDexists));
        if ~isempty(exchangeAdded)
            if j==1
                addedReactions.(strains(i))= rxn(j);
                addedReactions.(strains(i))(1,2) = rxn_formula(j);
            else
                addedReactions.(strains(i))(end+1,1) = rxn(j);
                addedReactions.(strains(i))(end,2) = rxn_formula(j);
            end
        end
    end

    %% Add transport reactions where needed

    % Find blocked exchange reactions
    metsToCheck = append('EX_',metabolomicsExMets.VMHID,'(e)');
    [BlockedRxns] = identifyFastBlockedRxns(model,metsToCheck);

    % Add transport reactions for these blocked reactions
    metList = erase(BlockedRxns,'EX_');
    metList = erase(metList,'(e)');

    % Check if there are already transport reactions for the metabolites 
    database=loadVMHDatabase;


    inDatabase=0;
    for j = 1:length(metList)

        % Check if transport reaction exists for this metabolite
        databaseTransport = matches(database.reactions(:,11),'Transport, extracellular') & ...
        contains(database.reactions(:,3),append(metList{j },'[c]')) & ...
        contains(database.reactions(:,3),append(metList{j },'[e]'));

        % Get reaction name and formula
        if 0%sum(databaseTransport)>0
            transportRxn = database.reactions(databaseTransport,1);
            rxnForm = database.reactions(databaseTransport, 3);
            inDatabase = inDatabase+1;
        else
            transportRxn = append(upper(metList(j )),'t');
            rxnForm = append(metList{j },'[c]  <=> ',metList{j },'[e]');
        end

        rxnForm = cellstr(rxnForm);

        for k = 1:length(transportRxn)
            % Get reaction formula
            reaction = transportRxn(k);
            formula = rxnForm(k);
            % Add reaction
            [model,rxnIDexists] = addReaction(model, char(reaction), 'reactionFormula', char(formula), 'geneRule', 'ManualMetaboliteGapfill');
            % Add to added reactions
            transportRxnAdded=setdiff(reaction,model.rxns(rxnIDexists));
            if ~isempty(transportRxnAdded)
                addedReactions.(strains(i))(end+1,1) = reaction;
                addedReactions.(strains(i))(end,2) = formula;
            end
        end
    end
    
    % Add medium info
    % Get dietary reactions
    idx = metabolomicsExMets.Consumed==1;
    DietRxns = cellstr(append('EX_',metabolomicsExMets.VMHID(idx),'(e)'));

    % Close all consumptions reactions
    model = changeRxnBounds(model, model.rxns(contains(model.rxns,'EX_')), 1000, 'u');
    model = changeRxnBounds(model, model.rxns(contains(model.rxns,'EX_')), 0, 'l');

    % Open diet reactions
    model = changeRxnBounds(model, model.rxns(matches(model.rxns,DietRxns)), -1000, 'l');
    model = changeRxnBounds(model, model.rxns(matches(model.rxns,DietRxns)), 0, 'u');

    % Check again for blocked reactions. Did this unblock anything?
    BlockedReactions.(strains(i)).Blocked = identifyFastBlockedRxns(model,metsToCheck);
    BlockedReactions.(strains(i)).Unblocked = setdiff(BlockedRxns, BlockedReactions.(strains(i)).Blocked);
    
    % Save models
    modelsRefined.(strains(i)) = model;
end
end