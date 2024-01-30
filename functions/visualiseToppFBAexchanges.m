function visualiseToppFBAexchanges(consumedMets, secretedMets)
% This function creates heatmaps visualising the top 25 most consumed and
% secreted metabolites in the pFBA flux predictions.
%
% USAGE: visualisepFBAexchanges(res, figTitle, y_title,figNumber)
%
% INPUT
% consumedMets          Top 25 most consumed metabolites and their average fluxes
%                       per genus
% secretedMets          Top 25 most secreted metabolites and their average fluxes
%                       per genus
% AUTHOR: Tim Hensen, 12/2023.

additionalTestPath = pwd;
fig = figure('Position',[28,-3,1200,800]);

tiledlayout(2,1);
nexttile
visualisepFBAexchanges(consumedMets, 'Predicted top consumed metabolites in a complex medium', 'Consumption in mmol/g dry weight/hour','A');
nexttile
visualisepFBAexchanges(secretedMets, 'Predicted top secreted metabolites in a complex medium', 'Secretion in mmol/g dry weight/hour','B');

% save figure
exportgraphics(fig,[additionalTestPath filesep 'genus_exchange_fluxes.pdf'],'BackgroundColor','none','ContentType','vector')
exportgraphics(fig,[additionalTestPath filesep 'genus_exchange_fluxes.png'],'Resolution',400)

% Create heatmaps
function visualisepFBAexchanges(res, figTitle, y_title,figNumber)
mat = table2array(res(:,2:end));
idzero = all(mat==0,2);
res = res(~idzero,:);
mat = mat(~idzero,:);

if figNumber=='A'
    mat = -1*mat;
end

% Set dimensions
y = height(res);
x = width(res)-1;
s = [y x];

xrange = [-1 1]; % imagesc only needs the endpoints
yrange = [-1 1];

dx = diff(xrange)/(s(2)-1);
dy = diff(yrange)/(s(1)-1);
xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,s(2)+1);
yg = linspace(yrange(1)-dy/2,yrange(2)+dy/2,s(1)+1);

%close all
%fig = figure('Position',[28,-3,1200,400]);
imagesc(xrange, yrange,mat);
hold on
hm = mesh(xg,yg,zeros(s+1));
hm.FaceColor = 'none';
hm.EdgeColor = 'k';

ylabel('Metabolite','FontWeight',"bold")
yticks(linspace(-1,1,y))
yticklabels(res.Properties.RowNames)
ax = gca;
ax.YAxis.TickLabelInterpreter = 'none';

xlabel('Genus','FontWeight',"bold")
xticks(linspace(-1,1,x))
xticklabels(res.Properties.VariableNames(2:end))

t = title(figTitle);
t.Position(2) = -1.07;
ax=gca;
ax.XTickLabelRotation=-25;
ax.TitleFontSizeMultiplier=1.9;
ax.FontSize=8;
%ax.FontWeight = 'bold';
ax.XAxis.TickLabelInterpreter = 'none';

colormap('parula')
if figNumber=='A'
    caxis([0 round(max(max(mat)),0)])
end
a=colorbar;
ylabel(a,y_title,'FontSize',12,'Rotation',270);
a.Label.Position(1) = 4;
a.Location = 'manual';
a.Position(1) = 0.92;

text(-1.1,-1.25,figNumber,'FontSize',20)
end
end