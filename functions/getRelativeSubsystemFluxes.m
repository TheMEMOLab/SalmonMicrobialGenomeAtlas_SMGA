function resultPath = getRelativeSubsystemFluxes(pFBA_path, resultFolder)
% This function calculates the relative abundances of the total fluxes per
% subsystem.Transport and exchange fluxes are excluded.
%
% USAGE: resultPath = getRelativeSubsystemFluxes(resPath, resultFolder)
%
% INPUT
% resPath               Path to the output file from the
%                       saveAndProcessFluxes function 
% resultFolder          Path to folder where the outputs will be saved
%
% OUTPUT
% resultPath            Path to file containing the function outputs
%
% AUTHOR: Tim Hensen, 12/2023.

% Load pFBA results
pFBA = load(pFBA_path);
pFBA = pFBA.(string(fieldnames(pFBA)));

% Get fields
fields = string(fieldnames(pFBA));

% Calculate the subsystem utulisation by 
% sum(Sub_flux)/sum(Total_flux) * length(Sub_rxns)/length(Total_rxns)
pFBA_norm = struct();
for i=1:length(fields)

    % Remove transport and demand reactions
    idrem = matches(pFBA.(fields(i)).Subsystem,...
    {'','Exchange/demand reaction','Transport, periplasmatic','Transport, extracellular','Exchange'});
    pFBA.(fields(i)) = pFBA.(fields(i))(~idrem,:);

    % Calculate the total flux in strain
    %t_flux = sum(abs(pFBA.(fields(i)).(fields(i))),'omitnan');

    % Get the subsystem groups
    [subsysGroups,subs] = findgroups(pFBA.(fields(i)).Subsystem);

    % Calculate the total flux per subsystem
    sumf = @(x) sum(abs(x),'omitnan');
    s_flux=splitapply(sumf,pFBA.(fields(i)).(fields(i)),subsysGroups);

    % Calculate the number of reactions in total and per subsystem
    %t_rxns = length(pFBA.(fields(i)).Reaction);

    % Find number of reactions per subsystem
    s_rxns=splitapply(@length,pFBA.(fields(i)).Reaction,subsysGroups);

    pFBA_table_norm = table(subs,'VariableNames',{'Subsystem'});

    % Calculate the total flux per reaction for each subsystem
    pFBA_table_norm.(fields(i)) = s_flux./s_rxns;

    % Calculate the relative total flux per reaction for each subsystem
    pFBA_table_norm.(fields(i)) = pFBA_table_norm.(fields(i)) ./ sum(pFBA_table_norm.(fields(i)));
    pFBA_table_norm.(fields(i)) = pFBA_table_norm.(fields(i)) * 100;
    pFBA_norm.(fields(i)) = pFBA_table_norm;
end

% Combine all results
fields = string(fieldnames(pFBA_norm));
for i=1:length(fields)
    if i == 1
        pFBA_table = pFBA_norm.(fields(i));
    else
        pFBA_table = outerjoin(pFBA_table,pFBA_norm.(fields(i)),'Keys',{'Subsystem'},'MergeKeys',true);
    end
end

% Save all results
resultPath = [resultFolder filesep 'relativeSubsystemUtulisation.xlsx'];
writetable(pFBA_table,resultPath)
end

