clear
close all

%% Paths and calibration inputs
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);

addpath(fullfile(repoRoot, 'ode'))
addpath(fullfile(repoRoot, 'data_rich'))

loadedResults = load(fullfile(repoRoot, 'data_rich', 'results_Danio_rerio.mat'), 'par');
par = loadedResults.par;
if isfield(par, 'free')
    par = rmfield(par, 'free');
end

[data, auxData] = mydata_Danio_rerio;

lawrHighData = data.tL_LawrEber2008_high;
lawrLowData = data.tL_LawrEber2008_low;

TC_LawrEber2008 = tempcorr(auxData.temp.tL_LawrEber2008_high, par.T_ref, par.T_A);
init_cond = [par.V_0; getE0(par.f, par, parscomp_st(par)); 0; 0; 1; 0];

tEnd = max([lawrHighData(:,1); lawrLowData(:,1)]) + 5;
tPlot = linspace(0, tEnd, 500)';

%% tL_LawrEber2008 at low and high
predLow = predictLengthCurve(tPlot, init_cond, par, par.f_LawrEber2008_low, TC_LawrEber2008);
predHigh = predictLengthCurve(tPlot, init_cond, par, par.f_LawrEber2008_high, TC_LawrEber2008);

%% Figure
fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [100 100 980 620]);
ax = axes(fig);
hold(ax, 'on')
set(ax, 'FontSize', 14, 'Position', [0.11 0.13 0.84 0.78])
ax.Toolbar.Visible = 'off';
ax.Interactions = [];

highColor = [0.1216 0.4667 0.7059];
lowColor = [0.8392 0.1529 0.1569];

plot(ax, tPlot, predHigh, 'LineWidth', 2.5, ...
    'Color', highColor, 'DisplayName', 'High food prediction')
plot(ax, tPlot, predLow, 'LineWidth', 2.5, ...
    'Color', lowColor, 'DisplayName', 'Low food prediction')

scatter(ax, lawrHighData(:,1), lawrHighData(:,2), 42, ...
    'o', 'MarkerEdgeColor', highColor, 'MarkerFaceColor', highColor, ...
    'DisplayName', 'High food data')
scatter(ax, lawrLowData(:,1), lawrLowData(:,2), 42, ...
    's', 'MarkerEdgeColor', lowColor, 'MarkerFaceColor', lowColor, ...
    'DisplayName', 'Low food data')

xlabel(ax, 'Time since fertilization (d)')
ylabel(ax, 'Total length (cm)')
title(ax, 'Lawrence et al. 2008 growth trajectories')
legend(ax, 'Location', 'northwest')
box(ax, 'on')
grid(ax, 'off')
xlim(ax, [0, tEnd])
yMin = min([lawrHighData(:,2); lawrLowData(:,2); predHigh(:); predLow(:)]);
yMax = max([lawrHighData(:,2); lawrLowData(:,2); predHigh(:); predLow(:)]);
yPad = 0.08 * (yMax - yMin);
ylim(ax, [max(0, yMin - yPad), yMax + yPad])

%% Save figure
figOutDir = fullfile(scriptDir, 'figures');
if ~exist(figOutDir, 'dir')
    mkdir(figOutDir);
end

exportgraphics(fig, fullfile(figOutDir, 'lawrence_et_al_2008_validation.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(figOutDir, 'lawrence_et_al_2008_validation.pdf'), 'ContentType', 'vector');

%% Prediction function
function pred = predictLengthCurve(curveData, init_cond, par, F, TC)
a = [0; curveData(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order

[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1);
L = V.^(1/3);

L = L(it_sort); L = L(2:end); % reconstuct L
pred = L / par.del_Mt;
end