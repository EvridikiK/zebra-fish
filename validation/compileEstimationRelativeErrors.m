function reTable = compileEstimationRelativeErrors(estimations, datasets)
%COMPILEESTIMATIONRELATIVEERRORS Compile relative errors across estimations.
%   reTable = compileEstimationRelativeErrors(estimations)
%   reTable = compileEstimationRelativeErrors(estimations, datasets)
%
%   Inputs
%   estimations : cell array or string array with estimation folder names.
%   datasets    : optional cell array or string array with dataset names.
%
%   Output
%   reTable     : table with datasets as rows and estimations as columns.
%                 Includes a final row named 'MRE'.

if nargin < 1 || isempty(estimations)
    error('Input ''estimations'' must be provided and cannot be empty.');
end

if isstring(estimations)
    estimations = cellstr(estimations);
end

useAllDatasets = (nargin < 2 || isempty(datasets));
if ~useAllDatasets && isstring(datasets)
    datasets = cellstr(datasets);
end
if ~useAllDatasets
    datasets = datasets(:)';
    datasets(strcmp(datasets, 'psd')) = [];
end

nEstimations = numel(estimations);
estimationData = cell(nEstimations, 1);
allDatasets = {};

for i = 1:nEstimations
    estimation = estimations{i};

    resultPath = fullfile('..', estimation, 'results_Danio_rerio.mat');
    if ~exist(resultPath, 'file')
        error('No results file found for estimation ''%s'' at ''%s''.', ...
            estimation, resultPath);
    end

    loaded = load(resultPath, 'metaPar', 'data');
    if ~isfield(loaded, 'metaPar') || ~isfield(loaded, 'data')
        error('File ''%s'' must contain both ''metaPar'' and ''data''.', resultPath);
    end

    if ~isfield(loaded.metaPar, 'RE') || ~isfield(loaded.metaPar, 'MRE')
        error('File ''%s'' missing ''metaPar.RE'' or ''metaPar.MRE''.', resultPath);
    end

    datasetNames = fieldnames(loaded.data);
    datasetNames(strcmp(datasetNames, 'psd')) = [];
    RE = loaded.metaPar.RE;
    if size(RE, 2) < 2
        error('metaPar.RE in ''%s'' must have at least two columns.', resultPath);
    end

    nDatasets = numel(datasetNames);
    if size(RE, 1) < nDatasets
        error(['metaPar.RE in ''%s'' has fewer rows (%d) than non-psd ', ...
            'datasets in data (%d).'], resultPath, size(RE, 1), nDatasets);
    end

    reValues = RE(1:nDatasets, 1);
    weights = RE(1:nDatasets, 2);

    keep = weights ~= 0 & ~strcmp(datasetNames, 'psd');
    datasetNames = datasetNames(keep);
    reValues = reValues(keep);

    if ~useAllDatasets
        wanted = ismember(datasetNames, datasets);
        datasetNames = datasetNames(wanted);
        reValues = reValues(wanted);
    else
        for k = 1:numel(datasetNames)
            if ~ismember(datasetNames{k}, allDatasets)
                allDatasets{end + 1} = datasetNames{k}; %#ok<AGROW>
            end
        end
    end

    estimationData{i}.datasetNames = datasetNames;
    estimationData{i}.reValues = reValues;
    estimationData{i}.MRE = loaded.metaPar.MRE;
end

if useAllDatasets
    finalDatasets = allDatasets;
else
    finalDatasets = datasets;
end

finalRows = [finalDatasets, {'MRE'}];
values = nan(numel(finalRows), nEstimations);

for i = 1:nEstimations
    names = estimationData{i}.datasetNames;
    vals = estimationData{i}.reValues;

    for j = 1:numel(finalDatasets)
        idx = find(strcmp(names, finalDatasets{j}), 1, 'first');
        if ~isempty(idx)
            values(j, i) = vals(idx);
        end
    end

    values(end, i) = estimationData{i}.MRE;
end

varNames = matlab.lang.makeValidName(estimations);
varNames = matlab.lang.makeUniqueStrings(varNames, {}, namelengthmax);

reTable = array2table(values, ...
    'VariableNames', varNames, ...
    'RowNames', finalRows);
end
