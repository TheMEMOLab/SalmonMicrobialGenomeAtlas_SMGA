function visualiseSubsystemFluxes(genusResultPath)
% This function visualises the summarised pFBA results
%
% USAGE: visualiseModelTests(modelPath, fvaMetabolomeTestPath)
%
% INPUT
% genusResultPath       Path to summariseSubSysUtulisation function
%                       output
%
% AUTHOR: Tim Hensen, 12/2023.

% Visualise pFBA subsystem utulisation results

%genusResultPath = 'genus_RelativeSubsystemUtulisation.xlsx';
% load results
res = readtable(genusResultPath);

% Remove all subsystems with zero for all genera
mat = table2array(res(:,2:end));
idzero = all(mat==0,2);
res = res(~idzero,:);
mat = mat(~idzero,:);

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

close all
fig = figure('Position',[28,-3,1200,1400]);
imagesc(xrange, yrange,mat);
hold on
hm = mesh(xg,yg,zeros(s+1));
hm.FaceColor = 'none';
hm.EdgeColor = 'k';

yticks(linspace(-1,1,y))
yticklabels(res.Subsystem)
set(groot,'defaultAxesTickLabelInterpreter','none');  

xticks(linspace(-1,1,x))
xticklabels(res.Properties.VariableNames(2:end))

t = title('Relative flux activity by subsystem');
t.Position(2) = -1.03;
ax=gca;
ax.XTickLabelRotation=-25;
ax.TitleFontSizeMultiplier=1.9;
ax.FontSize=8;
ax.FontWeight = 'bold';

a=colorbar;
ylabel(a,'Percentage (%) of total normalised predicted flux in subsystem','FontSize',12,'Rotation',270);
a.Label.Position(1) = 4;
a.Location = 'manual';
a.Position(1) = 0.92;

% save figure
exportgraphics(fig,[pwd filesep 'analysis' filesep 'genus_total_fluxes.pdf'],'BackgroundColor','none','ContentType','vector')
exportgraphics(fig,[pwd filesep 'analysis' filesep 'genus_total_fluxes.png'],'Resolution',400)


