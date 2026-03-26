clear
close all

%% Figure options
transparentBackground = true;
figureWidth = 1220;
figureHeight = 680;
axisFontSize = 18;
labelFontSize = 21;
legendFontSize = 17;
titleFontSize = 21;
exportPadding = 60;

if transparentBackground
    displayBackgroundColor = 'none';
else
    displayBackgroundColor = 'w';
end

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

yMin = min([lawrHighData(:,2); lawrLowData(:,2)]);
yMax = max([lawrHighData(:,2); lawrLowData(:,2)]);

%% tL_LawrEber2008 at low and high
predLow = predictLengthCurve(tPlot, init_cond, par, par.f_LawrEber2008_low, TC_LawrEber2008);
predHigh = predictLengthCurve(tPlot, init_cond, par, par.f_LawrEber2008_high, TC_LawrEber2008);
yMin = min([yMin; predHigh(:); predLow(:)]);
yMax = max([yMax; predHigh(:); predLow(:)]);
yPad = 0.08 * (yMax - yMin);

highColor = [0.1216 0.4667 0.7059];
lowColor = [0.8392 0.1529 0.1569];

%% Figure with predictions
fig = figure('Color', displayBackgroundColor, 'Units', 'pixels', 'Position', [100 100 figureWidth figureHeight]);
tl = tiledlayout(fig, 1, 1, 'Padding', 'loose', 'TileSpacing', 'loose');
ax = nexttile(tl);
hold(ax, 'on')
set(ax, 'FontSize', axisFontSize, 'Color', displayBackgroundColor)
ax.Toolbar.Visible = 'off';
ax.Interactions = [];

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

xlabel(ax, 'Time since fertilization (d)', 'FontSize', labelFontSize)
ylabel(ax, 'Total length (cm)', 'FontSize', labelFontSize)
sgtitle(tl, 'Growth trajectories from Lawrence et al. 2008', 'FontSize', titleFontSize)
legend(ax, 'Location', 'northwest', 'FontSize', legendFontSize)
box(ax, 'off')
grid(ax, 'off')
xlim(ax, [0, tEnd])
ylim(ax, [max(0, yMin - yPad), yMax + yPad])

%% Figure with data only
dataFig = figure('Color', displayBackgroundColor, 'Units', 'pixels', 'Position', [100 100 figureWidth figureHeight]);
dataTl = tiledlayout(dataFig, 1, 1, 'Padding', 'loose', 'TileSpacing', 'loose');
dataAx = nexttile(dataTl);
hold(dataAx, 'on')
set(dataAx, 'FontSize', axisFontSize, 'Color', displayBackgroundColor)
dataAx.Toolbar.Visible = 'off';
dataAx.Interactions = [];

scatter(dataAx, lawrHighData(:,1), lawrHighData(:,2), 42, ...
    'o', 'MarkerEdgeColor', highColor, 'MarkerFaceColor', highColor, ...
    'DisplayName', 'High food data')
scatter(dataAx, lawrLowData(:,1), lawrLowData(:,2), 42, ...
    's', 'MarkerEdgeColor', lowColor, 'MarkerFaceColor', lowColor, ...
    'DisplayName', 'Low food data')

xlabel(dataAx, 'Time since fertilization (d)', 'FontSize', labelFontSize)
ylabel(dataAx, 'Total length (cm)', 'FontSize', labelFontSize)
sgtitle(dataTl, 'Growth data from Lawrence et al. 2008', 'FontSize', titleFontSize)
legend(dataAx, 'Location', 'northwest', 'FontSize', legendFontSize)
box(dataAx, 'off')
grid(dataAx, 'off')
xlim(dataAx, [0, tEnd])
ylim(dataAx, [max(0, yMin - yPad), yMax + yPad])

%% Save figures
figOutDir = fullfile(scriptDir, 'figures');
if ~exist(figOutDir, 'dir')
    mkdir(figOutDir);
end

saveFigurePair(fig, tl, ax, fullfile(figOutDir, 'lawrence_et_al_2008_predictions.png'), fullfile(figOutDir, 'lawrence_et_al_2008_predictions.pdf'), transparentBackground, exportPadding);
saveFigurePair(dataFig, dataTl, dataAx, fullfile(figOutDir, 'lawrence_et_al_2008_data_only.png'), fullfile(figOutDir, 'lawrence_et_al_2008_data_only.pdf'), transparentBackground, exportPadding);

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

%% Save helper
function saveFigurePair(figHandle, layoutHandle, axesHandle, pngPath, pdfPath, transparentBackground, exportPadding)
if transparentBackground
    saveTransparentPng(figHandle, layoutHandle, axesHandle, pngPath, exportPadding);
    exportgraphics(layoutHandle, pdfPath, 'ContentType', 'vector', 'Padding', exportPadding, 'BackgroundColor', 'w');
else
    exportgraphics(layoutHandle, pngPath, 'Resolution', 300, 'Padding', exportPadding, 'BackgroundColor', 'w');
    exportgraphics(layoutHandle, pdfPath, 'ContentType', 'vector', 'Padding', exportPadding, 'BackgroundColor', 'w');
end
end

%% Transparent PNG helper
function saveTransparentPng(figHandle, layoutHandle, axesHandle, pngPath, exportPadding)
keyColor = [1 0 1];
tmpPngPath = [pngPath '.tmp.png'];

originalFigureColor = figHandle.Color;
originalAxesColor = axesHandle.Color;
figHandle.Color = keyColor;
axesHandle.Color = keyColor;

cleanupObj = onCleanup(@() restoreColors(figHandle, axesHandle, originalFigureColor, originalAxesColor));
exportgraphics(layoutHandle, tmpPngPath, 'Resolution', 300, 'Padding', exportPadding, 'BackgroundColor', keyColor);
makePngBackgroundTransparent(tmpPngPath, pngPath, keyColor);
delete(tmpPngPath);
clear cleanupObj
end

%% Restore helper
function restoreColors(figHandle, axesHandle, figureColor, axesColor)
figHandle.Color = figureColor;
axesHandle.Color = axesColor;
end

%% Alpha helper
function makePngBackgroundTransparent(inputPath, outputPath, keyColor)
rgbImage = im2double(imread(inputPath));
if size(rgbImage, 3) == 1
    rgbImage = repmat(rgbImage, [1 1 3]);
end

[height, width, ~] = size(rgbImage);
backgroundImage = repmat(reshape(keyColor, 1, 1, 3), height, width);
distanceToBackground = sqrt(sum((rgbImage - backgroundImage) .^ 2, 3));
alphaThreshold = 0.35;
alphaChannel = min(1, distanceToBackground / alphaThreshold);
alphaChannel(distanceToBackground < 1e-8) = 0;

alphaImage = repmat(max(alphaChannel, eps), 1, 1, 3);
foregroundImage = (rgbImage - (1 - alphaImage) .* backgroundImage) ./ alphaImage;
foregroundImage = min(max(foregroundImage, 0), 1);
foregroundImage(repmat(alphaChannel == 0, 1, 1, 3)) = 0;

imwrite(foregroundImage, outputPath, 'Alpha', alphaChannel);
end

