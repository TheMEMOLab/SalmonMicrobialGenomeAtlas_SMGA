function [modelsGapfilled,gapfillRxns,blockedRxnStruct] = automaticGapfilling(modelsRefined,blockedRxnStruct)
% This function performs automic gapfilling on the refined reconstructions using the fastGapFill function.
%
% USAGE: [modelsGapfilled,gapfillRxns,blockedRxnStruct] = automaticGapfilling(modelsRefined,blockedRxnStruct)
%
% INPUT
% modelsRefined         Structure with refined reconstructions
% blockedRxnStruct      Structure with blocked reactions for each
%                       reconstruction for which exomometabolomics 
%                       data was available
%
% OUTPUT
% modelsGapfilled       Structure with gapfilled refined reconstructions
% gapfillRxns           Structure with previously blocked gapfilled reactions
% blockedRxnStruct      Structure with remaining blocked reactions
%
% AUTHOR: Tim Hensen, 12/2023.

strains = string(fieldnames(modelsRefined));

% Preallocate outputs
modelsGapfilled = struct();
gapfillRxns = struct();

for i=1:length(strains)
    % Get model
    model = modelsRefined.(strains(i));

    % Find blocked reactions
    metsToCheck = blockedRxnStruct.(strains(i)).Blocked;
    [Blocked_before] = identifyFastBlockedRxns(model,metsToCheck);


    % Gapfill model
    disp('start prepareFastGapFill')

    % Obtain consistMatricesSUX
    [consistModel,consistMatricesSUX,BlockedRxns] = prepareFastGapFill(model);

    % Set parameters
    epsilon = 1e-4;
    weights.MetabolicRxns = 0.1; % Kegg metabolic reactions
    weights.ExchangeRxns = 0.5; % Exchange reactions
    weights.TransportRxns = 10; % Transport reactions

    disp('start fastGapFill')
    % Perform gapfill
    [AddedRxns] = fastGapFill(consistMatricesSUX,epsilon, weights);

    % Annotate reactions
    disp('start postProcessGapFillSolutions')
    [AddedRxnsExtended] = postProcessGapFillSolutions(AddedRxns,consistModel,BlockedRxns);

    % Add reactions to model
    disp('Add gapfilled reactions')
    for j=1:length(AddedRxnsExtended.rxns)
        model = addReaction(model,AddedRxnsExtended.rxns{j}, AddedRxnsExtended.rxnFormula{j});
    end

    % Find unblocked reactions and blocked reactions
    Blocked_after = identifyFastBlockedRxns(model,metsToCheck);
    blockedRxnStruct.(strains(i)).Blocked = Blocked_after;
    blockedRxnStruct.(strains(i)).Unblocked = setdiff(Blocked_before,Blocked_after);

    % Save added reactions
    gapfillRxns.(strains(i)) = AddedRxnsExtended;
    % Save model
    modelsGapfilled.(strains(i)) = model;
end
end
