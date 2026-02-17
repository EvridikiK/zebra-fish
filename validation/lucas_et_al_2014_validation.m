clear
addpath('../ode/')

%% ---------- PLOT OPTIONS ----------
% Toggle grid for ALL plots
opts.gridOn = false;  

% Define colors for each estimation (edit to taste)
% Must be Nx3 RGB in [0,1], same order as `estimations`
modelColors = [
    0.1216 0.4667 0.7059;   % data_rich (seaborn blue)
    1.0000 0.4980 0.0549;   % data_moderate (seaborn orange)
    0.1725 0.6275 0.1725;   % data_limited (seaborn green)
];

% Data styling
legendLabelData = 'Lucas et. al 2014 data'; 
dataColor = [1 0 0];  % red
markerAlpha = 0.5;   % transparency for overlapping points

%% ---------- SETTINGS ---------
% Estimation folders to compare
estimations = {'data_rich', 'data_moderate', 'data_limited'};
% Food levels to compare
F = 0.7;

% Save options
saveFigs = false;
saveFormats = {'png', 'pdf'};    
figOutDir = fullfile(pwd, 'figures');

%% Load oxygen vs weight data
dataFolder = fullfile('..', 'data', 'Lucas et al. 2014');

larvae_data   = readtable(fullfile(dataFolder, '5-day_larvae.csv'));
juvenile_data = readtable(fullfile(dataFolder, '2-month_juveniles.csv'));
adult_data    = readtable(fullfile(dataFolder, '6-month_adults.csv'));

o2_vs_weight = [ ...
    larvae_data.weight,   larvae_data.oxygen_consumption;
    juvenile_data.weight, juvenile_data.oxygen_consumption;
    adult_data.weight,    adult_data.oxygen_consumption
];

o2_vs_age = [
    5   mean(larvae_data.oxygen_consumption);
    60  mean(juvenile_data.oxygen_consumption);
    180 mean(adult_data.oxygen_consumption)
];

age_length_weight_data = readtable(fullfile(dataFolder, 'age_length_weight.csv'));
weight_vs_age = [age_length_weight_data.age, age_length_weight_data.weight];
length_vs_age = [age_length_weight_data.age, age_length_weight_data.length];

% Standard errors for error bars (from table columns)
se_weight_vs_age = age_length_weight_data.se_weight;
se_length_vs_age = age_length_weight_data.se_length;

%% ---------- Initialize figures once (plot DATA once) ----------

% Figure 1: Oxygen vs weight (log-log)
figure(1); clf; hold on
set(gca, 'FontSize', 14, 'XScale', 'log', 'YScale', 'log')
h1 = scatter(o2_vs_weight(:,1), o2_vs_weight(:,2), 45, ...
    'MarkerFaceColor', dataColor, 'MarkerEdgeColor', dataColor, ...
    'MarkerFaceAlpha', markerAlpha, 'MarkerEdgeAlpha', markerAlpha, ...
    'DisplayName', legendLabelData);
xlabel('Weight (g)')
ylabel('Oxygen consumption (mg O_2/h)')

% Figure 2: Oxygen vs age
figure(2); clf; hold on
set(gca, 'FontSize', 14)
h2 = scatter(o2_vs_age(:,1), o2_vs_age(:,2), 55, ...
    'MarkerFaceColor', dataColor, 'MarkerEdgeColor', dataColor, ...
    'DisplayName', legendLabelData);
xlabel('Age (d)')
ylabel('Oxygen consumption (mg O_2/h)')

% Figure 3: Weight vs age
figure(3); clf; hold on
set(gca, 'FontSize', 14)

% Data with error bars (standard error)
h3 = errorbar(weight_vs_age(:,1), weight_vs_age(:,2), se_weight_vs_age, 'o', ...
    'Color', dataColor, 'MarkerFaceColor', dataColor, 'MarkerEdgeColor', dataColor, ...
    'LineStyle', 'none', 'CapSize', 6, ...
    'DisplayName', legendLabelData);

xlabel('Age (d)')
ylabel('Weight (g)')

% Figure 4: Length vs age
figure(4); clf; hold on
set(gca, 'FontSize', 14)

% Data with error bars (standard error)
h4 = errorbar(length_vs_age(:,1), length_vs_age(:,2), se_length_vs_age, 'o', ...
    'Color', dataColor, 'MarkerFaceColor', dataColor, 'MarkerEdgeColor', dataColor, ...
    'LineStyle', 'none', 'CapSize', 6, ...
    'DisplayName', legendLabelData);

xlabel('Age (d)')
ylabel('Length (mm)')

%% ---------- Loop over estimations: compute + plot predictions ----------
for i = 1:numel(estimations)
    estimation = estimations{i};

    % Load parameters for this estimation
    par = loadEstimationParameters(estimation);
    cpar = parscomp_st(par);
    vars_pull(par); vars_pull(cpar);

    % Simulate from birth
    [t, L, W, J_O] = getPredictions(par, cpar, F);

    % Pretty legend label: "Data rich" / "Data moderate" / "Data limited"
    prettyName = formatEstimationName(estimation);

    % Add prediction line to each plot
    figure(1);
    plot(W, J_O, 'LineWidth', 2, 'Color', modelColors(i,:), ...
        'DisplayName', prettyName)

    figure(2);
    plot(t, J_O, 'LineWidth', 2, 'Color', modelColors(i,:), ...
        'DisplayName', prettyName)

    figure(3);
    plot(t, W, 'LineWidth', 2, 'Color', modelColors(i,:), ...
        'DisplayName', prettyName)

    figure(4);
    plot(t, L, 'LineWidth', 2, 'Color', modelColors(i,:), ...
        'DisplayName', prettyName)
end

%% ---------- Finalize legends + grid + save figures ----------

% Output folder for figures (relative to this script)
if ~exist(figOutDir, 'dir')
    mkdir(figOutDir);
end

for f = 1:4
    figure(f);
    legend('Location','NorthWest');
    if opts.gridOn
        grid on
    else
        grid off
    end

    if saveFigs
        % Keep filenames stable and meaningful
        switch f
            case 1
                baseName = 'o2_vs_weight';
            case 2
                baseName = 'o2_vs_age';
            case 3
                baseName = 'weight_vs_age';
            case 4
                baseName = 'length_vs_age';
        end
        baseName = ['lucas_et_al_2014_' baseName];

        for k = 1:numel(saveFormats)
            ext = lower(saveFormats{k});
            outPath = fullfile(figOutDir, [baseName '.' ext]);

            switch ext
                case 'png'
                    exportgraphics(gcf, outPath, 'Resolution', 300);
                case 'pdf'
                    exportgraphics(gcf, outPath, 'ContentType', 'vector');
                case 'fig'
                    savefig(gcf, outPath);
                otherwise
                    error('Unsupported save format: %s', ext);
            end
        end
    end
end

%% ---------- Helper functions ----------

function prettyName = formatEstimationName(estimation)
% Converts:
%   data_rich     -> Data rich estimation
%   data_moderate -> Data moderate estimation
%   data_limited  -> Data limited estimation
s = strrep(estimation, '_', ' ');
core = strrep(s, 'data ', ''); % remove leading "data " when present
prettyName = ['Data ' core ' estimation'];
end

function par = loadEstimationParameters(estimation)
currentFolder = pwd;
cleanup = onCleanup(@() cd(currentFolder));


load(fullfile('..', estimation, 'results_Danio_rerio.mat'), 'par')
% cd(fullfile('..', estimation))
% [~, ~, metaData, ~, ~] = mydata_Danio_rerio;
% [par, ~, ~] = pars_init_Danio_rerio(metaData);
end

function [t_out, L, W, J_O] = getPredictions(par, cpar, F)
% Simulate from fertilization
TC = tempcorr(C2K(28), par.T_ref, par.T_A);
E_0 = getE0(F, par, cpar);
init_cond = [1e-10; E_0; 0; 0; 1; 0];

t = linspace(0, 200, 201);
[~, VEHRsMG] = ode45(@ode_VEHRsMG, t, init_cond, [], par, F, TC);

% Keep consistency with your original indexing (drop t=0 row)
t_out = t(2:end);

V   = VEHRsMG(2:end, 1);
E   = VEHRsMG(2:end, 2);
E_H = VEHRsMG(2:end, 3);
E_R = VEHRsMG(2:end, 4);
s_M = VEHRsMG(2:end, 5);
% E_egg = VEHRsMG(2:end, 6); % not used here

[p_A, ~, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M, TC, F, par);
p_D = compute_dissipation_power(p_S, p_J, p_R, p_C2, E_H, par.E_Hp, par.kap_R);

eta_M = -inv(cpar.n_M) * cpar.n_O * cpar.eta_O;
del_Ms = par.del_Mt * 1.25; % shape coefficient for standard length
L = V.^(1/3) / del_Ms * 10; % mm

J_O = [p_A, p_D, p_G] * eta_M(3, :)'; % mol/d
J_O = -J_O * 32 * 1e3 / 24;          % mg O2/h

W = 1 * V + (E + E_R) * cpar.w_E / par.mu_E / par.d_E; % g
end
