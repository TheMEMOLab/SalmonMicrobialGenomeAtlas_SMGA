function fileName = getMetabolomeUptakeSecretionInfo(mediapath,mappingDataPath)
% This function loads the exometabolomics data and obtains for each mapped
% metabolite if it is secreted or consumed. 
%
% USAGE: fileName = getMetabolomeUptakeSecretionInfo(mediapath,mappingDataPath)
%
% INPUT
% mediapath             Path to exometabolomic data file
% mappingDataPath       Path to file with mapped metabolites. This file must be the output from the
%                       mapIDtoVMHID function      
%
% OUTPUT
% fileName              Path to file with additional reactions to add to
%                       the reconstructions
%
% AUTHOR: Tim Hensen, 12/2023.

% Obtain microbe names
mapping = readtable(mediapath,'ReadVariableNames',true,'Sheet','Mapping');

% Get mapped data
mappingData = readtable(mappingDataPath);
mappingData(~mappingData.Mapped,:)=[];

additionalRxns = struct();
for i=1:length(mapping.Microbe)
    % For each microbe and medium, obtain which metabolites are secreted and consumed. 
    mediadata = readtable(mediapath,'ReadVariableNames',true,'Sheet',mapping.Microbe{i});
    
    % debug
    mappingData = mappingData(string(mappingData.metabolite) ~= "Urocanate",:);
    mediadata = mediadata(string(mediadata.Name) ~= "Urocanic acid",:);
    
    % Filter on metabolites which were mapped onto VMH IDs
    mediadata = innerjoin(mediadata,mappingData);
    mediadata_mapped = mediadata(~matches(mediadata.VMHID,""),:);
    
    % Add which metabolites are consumed or secreted
    mediadata_mapped.significant = mediadata_mapped.p<mediadata_mapped.FDR; % Check if p-value is below fdr threshold
    mediadata_mapped.Consumed = mediadata_mapped.log2_Control_bacterium>0 & mediadata_mapped.significant;
    mediadata_mapped.Secreted = mediadata_mapped.log2_Control_bacterium<0 & mediadata_mapped.significant;
    
    % Store results
    microbename = mapping.Strain{i};
    additionalRxns.(microbename) = mediadata_mapped(:,{'VMHID','VMH_metabolite','Consumed','Secreted'});
end
fileName = 'additionalRxns.mat';
save(fileName,'additionalRxns')