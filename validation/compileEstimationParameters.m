function paramTable = compileEstimationParameters(estimations, parameters)
%COMPILEESTIMATIONPARAMETERS Compile estimated parameters across estimations.
%   paramTable = compileEstimationParameters(estimations)
%   paramTable = compileEstimationParameters(estimations, parameters)
%
%   Inputs
%   estimations : cell array or string array with estimation folder names.
%   parameters  : cell array or string array with parameter names.
%
%   Output
%   paramTable  : table with parameters as rows and estimations as columns.

if nargin < 1 || isempty(estimations)
    error('Input ''estimations'' must be provided and cannot be empty.');
end

if nargin < 2 || isempty(parameters)
    parameters = {'z', 'kap_X', 'kap_P', 'v', 'kap', 'p_M', ...
        'E_G', 'k_J', 'E_Hb', 'E_Hj', 'E_Hp', 'kap_R'};
end

if isstring(estimations)
    estimations = cellstr(estimations);
end
if isstring(parameters)
    parameters = cellstr(parameters);
end

nEstimations = numel(estimations);
nParameters = numel(parameters);
values = nan(nParameters, nEstimations);

for i = 1:nEstimations
    estimation = estimations{i};

    resultPath = fullfile('..', estimation, 'results_Danio_rerio.mat');
    if ~exist(resultPath, 'file')
        error('No results file found for estimation ''%s'' at ''%s''.', ...
            estimation, resultPath);
    end

    loaded = load(resultPath, 'par');
    if ~isfield(loaded, 'par')
        error('File ''%s'' does not contain a ''par'' struct.', resultPath);
    end
    par = loaded.par;

    for j = 1:nParameters
        paramName = parameters{j};
        if isfield(par, paramName)
            values(j, i) = par.(paramName);
        else
            warning('Parameter ''%s'' not found in estimation ''%s''.', ...
                paramName, estimation);
        end
    end
end

varNames = matlab.lang.makeValidName(estimations);
varNames = matlab.lang.makeUniqueStrings(varNames, {}, namelengthmax);

paramTable = array2table(values, ...
    'VariableNames', varNames, ...
    'RowNames', parameters);
end
