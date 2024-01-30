function visualiseModelTests(modelPath, fvaMetabolomeTestPath)
% This function loads the refined reconstructions and performs automic
% gapfilling using the fastGapFill function.
%
% USAGE: visualiseModelTests(modelPath, fvaMetabolomeTestPath)
%
% INPUT
% modelPath             Path to models to visualise
% fvaMetabolomeTestPath Path to qualitative comparisons between FVA
%                       exchange fluxes and exomometabolomics
%
% AUTHOR: Tim Hensen, 12/2023.


% INPUT: modelPath: Path to refined models
% fvaMetabolomeTestPath: Path to file with model test results

folder = what(modelPath);

% Get strains and microbes
strains = string(erase(folder.mat,'.mat'));
microbes = extractBefore(strains,'_');

% Get fva results 
resPath = [fvaMetabolomeTestPath filesep 'FVA_metabolomicsCheck.xlsx'];

% Obtain FVA results 
for i=1:length(microbes)
    result = readtable(resPath,"Sheet",microbes(i));
    result = result(:,{'VMH_metabolite','Consumption_in','Secretion_in'});
    result.Properties.VariableNames = cellstr(["metabolite",append(strains(i),'_consumption'),append(strains(i),'_secretion')]);

    if i == 1
        consumption = result(:,[1 2]);
        secretion = result(:,[1 3]);
    else
        consumption = outerjoin(consumption,result(:,[1 2]),'Keys',{'metabolite'},'MergeKeys',true);
        secretion = outerjoin(secretion,result(:,[1 3]),'Keys',{'metabolite'},'MergeKeys',true);
    end
end

% Remove appendices on names
consumption.Properties.VariableNames = erase(consumption.Properties.VariableNames,'_consumption');
secretion.Properties.VariableNames = erase(secretion.Properties.VariableNames,'_secretion');

FVA = {consumption, secretion};
% Remove metabolites not exchanged in the metabolome
for i = 1:length(FVA)
    data=table2array(FVA{i}(:,2:end));

    % Set empty cells to not secreted or consumed
    for j = 1:width(data)
        if i ==1
            data(cellfun(@isempty,data(:,j)),j) = {'not consumed'};
        elseif i==2
            data(cellfun(@isempty,data(:,j)),j) = {'not secreted'};
        end

        FVA{i}{:,2:end} = data;
    end

    % Remove metabolites not secreted or consumed
    rmRows = ~any(contains(data,{'metabolome and model','metabolome'}),2);
    FVA{i}(rmRows,:)=[];
end

% Encode data in numbers for creating the heatmaps
for i = 1:length(FVA)
    for j = 1:length(strains)
        col = FVA{i}.(strains(j));
        col(matches(col,'not consumed')) = {'1'};
        col(matches(col,'not secreted')) = {'1'};
        col(matches(col,'metabolome')) = {'2'};
        %col(matches(col,'model')) = {'3'};
        col(matches(col,'metabolome and model')) = {'3'};
        col = str2double(string(cell2mat(col)));
        FVA{i}.(strains(j)) = col;
    end
end

% Setup encodings
direction = {'Metabolites consumed in vitro','Metabolites secreted in vitro'};

for i=1:length(FVA)
    % Get data
    data_table = FVA{i};
    data_table = sortrows(data_table,'metabolite','ascend');
    %data_table = sortrows(data_table,'Photobacterium_phosphoreum_S39_bc34','descend');
    data=table2array(data_table(:,2:end));
    
    % Setup mesh
    s = [height(data_table) width(data_table)-1];
    xrange = [-1 1]; % imagesc only needs the endpoints
    yrange = [-1 1];
    dx = diff(xrange)/(s(2)-1);
    dy = diff(yrange)/(s(1)-1);
    xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,s(2)+1);
    yg = linspace(yrange(1)-dy/2,yrange(2)+dy/2,s(1)+1);
    
    % Visualise data
    fig = figure('Position',[28,-3,800,(height(data_table)*16)+40],'Visible',"on");
    h = imagesc(xrange, yrange, data);
    
    % Add mesh
    hold on
    hm = mesh(xg,yg,zeros(s+1));
    hm.FaceColor = 'none';
    hm.EdgeColor = 'k';
    hold off
    
    % Set title
    title(direction{i});
    
    % Set font sizes
    ax = ancestor(h, 'axes');
    ax.FontSize=8;
    ax.TitleFontSizeMultiplier=1.7;
    
    % Set y-axis
    yrule = ax.YAxis;
    yrule.TickValues = linspace(-1,1,height(data_table));
    yrule.TickLabels = data_table.metabolite;
    yrule.TickLabelInterpreter = 'none';
    
    % Set x-axis
    xrule = ax.XAxis;
    xrule.TickValues = linspace(-1,1,width(data_table)-1);
    xrule.TickLabels = data_table.Properties.VariableNames(2:end);
    xrule.TickLabelInterpreter = 'none';
    
    % Add colormap
    c = flipud(parula(3));
    colormap(c)
    caxis([0.5 3.5])
    
    % Set legend
    a=colorbar;
    a.Location = 'southoutside';
    a.Ticks = linspace(1,3,3);
    if i == 1
        a.TickLabels = {'Not consumed','metabolome','metabolome & model'};
        t = a.Title;
        t.String = "Metabolite consumed in:";
        t.FontWeight = 'bold';
        t.Position = [221.2500 -30 0];
    else
        a.TickLabels = {'not secreted','metabolome','metabolome & model'};
        t = a.Title;
        t.String = "Metabolite secreted in:";
        t.FontWeight = 'bold';
        t.Position = [221.2500 -30 0];
    end
    a.FontSize = 9;
    
    % Save figures
    exportgraphics(fig,[fvaMetabolomeTestPath filesep strcat(direction{i},'_metabolites.pdf')],'BackgroundColor','none','ContentType','vector')
    exportgraphics(fig,[fvaMetabolomeTestPath filesep strcat(direction{i},'_metabolites.png')],'Resolution',400)
end
end