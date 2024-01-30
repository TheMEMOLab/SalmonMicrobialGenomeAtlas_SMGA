function [modelsManGF,blockedRxnManStruct] = manualGapfilling(modelsGapfilled,blockedRxnGapFillStruct)
% This function performs manual gapfilling for remaining blocked reactions
%
% USAGE: [modelsManGF,blockedRxnManStruct] = manualGapfilling(modelsGapfilled,blockedRxnGapFillStruct)
%
% INPUT
% modelsGapfilled           Models gapfilled by automatic gapfilling
% blockedRxnGapFillStruct   Structure with remaining blocked reactions
% OUTPUT
% modelsManGF               Manually gapfilled models
% blockedRxnManStruct       Remaining blocked reactions
%
% AUTHOR: Tim Hensen, 12/2023.

% Manually gapfill blocked exchange reactions in the metabolome data

% Get strains
strains = string(fieldnames(modelsGapfilled));

% Get all unique blocked reactions
blockedRxns = {};
for i = 1:length(fieldnames(blockedRxnGapFillStruct))
    blocked = blockedRxnGapFillStruct.(strains(i)).Blocked;
    blockedRxns = [blockedRxns; blocked];
end
blockedRxns = unique(blockedRxns);

% Find blocked metabolites
blockedRxns(:,2) = extractBetween(blockedRxns,'EX_','(e)');

% Manually find reactions in the VMH database
blockedRxns{matches(blockedRxns(:,2),'btn'),3} = {'BTS4','BACCLi'};
blockedRxns{matches(blockedRxns(:,2),'glutar'),3} = {'OXPTNDH'};
blockedRxns{matches(blockedRxns(:,2),'ind3ac'),3} = {'L_TRPCOO','INHINDH'};
blockedRxns{matches(blockedRxns(:,2),'trp_L'),3} = {'TRPAS2','LTDCL','TRPS2r','TRPS1'}; % check
blockedRxns{matches(blockedRxns(:,2),'leu_L'),3} = {'LEUTA','LEUO','LLEUDr','ALALEU1c'};
blockedRxns{matches(blockedRxns(:,2),'xyl_D'),3} = {'XNOX'};

% Remove all other empty cells
blockedRxns(cellfun(@isempty, blockedRxns(:,3)),:)=[];

% Get VMH database
database = loadVMHDatabase;

% Preallocate outputs
modelsManGF = struct();
blockedRxnManStruct = struct();

for i=1:length(strains)
    % Get model
    model = modelsGapfilled.(strains(i));
    
    % Add sink for thm
    model = addSinkReactions(model, {'thm[c]'});

    % Find reactions to be gapfilled
    Blocked_before = blockedRxnGapFillStruct.(strains(i)).Blocked;
    [~,id_blocked] = setdiff(blockedRxns(:,1),Blocked_before);
    blockedRxns_2 = blockedRxns(id_blocked,:);
    
    % Check for each blocked reaction if manual gapfilling helps
    for j = 1:length(blockedRxns_2(:,1))
        gapfillRxns = blockedRxns_2{j,3};

        % Add gapfill reactions
        for k = 1:length(gapfillRxns)
            % Get reaction formula
            reaction = gapfillRxns(k);
            formula = database.reactions{ismember(database.reactions(:, 1), reaction), 3};
            % Add reaction
            model = addReaction(model, char(reaction), 'reactionFormula', formula, 'geneRule', 'ManualMetaboliteGapfill');
        end

        % Check if gapfilling worked
        %gapfillWorked = identifyFastBlockedRxns(model,blockedRxns_2(j,1));
    end
    
    % Check unblocked reactions
    Blocked_after = identifyFastBlockedRxns(model,Blocked_before);
    blockedRxnManStruct.(strains(i)).Blocked = Blocked_after;
    blockedRxnManStruct.(strains(i)).Unblocked = setdiff(Blocked_before,Blocked_after);
    
    % Save updated model
    modelsManGF.(strains(i)) = model;
end

end