function outTable = compileEstimationDataPredictionRelativeErrors(estimation, datasets, skipMultivariate, addUnits, addSources)
%COMPILEESTIMATIONDATAPREDICTIONRELATIVEERRORS Compile data/prediction/RE.
%   outTable = compileEstimationDataPredictionRelativeErrors(estimation)
%   outTable = compileEstimationDataPredictionRelativeErrors(estimation, datasets)
%   outTable = compileEstimationDataPredictionRelativeErrors(estimation, datasets, skipMultivariate)
%   outTable = compileEstimationDataPredictionRelativeErrors(estimation, datasets, skipMultivariate, addUnits, addSources)
%
%   Inputs
%   estimation      : estimation folder name (char or string).
%   datasets        : optional dataset names (cell array or string array).
%                     Defaults to all datasets in data (excluding 'psd').
%   skipMultivariate: optional logical flag. If true, multivariate datasets
%                     are removed from the output table. Default is false.
%   addUnits        : optional logical flag. If true, add Units column.
%                     Default is true.
%   addSources      : optional logical flag. If true, add Source column.
%                     Default is true.
%
%   Output
%   outTable        : table with row names as datasets and columns:
%                     Data, Prediction, RelativeError and optional
%                     Units/Source.

if nargin < 1 || isempty(estimation)
    error('Input ''estimation'' must be provided and cannot be empty.');
end
if nargin < 2
    datasets = [];
end
if nargin < 3 || isempty(skipMultivariate)
    skipMultivariate = false;
end
if nargin < 4 || isempty(addUnits)
    addUnits = true;
end
if nargin < 5 || isempty(addSources)
    addSources = true;
end

if ~(ischar(estimation) || (isstring(estimation) && isscalar(estimation)))
    error('Input ''estimation'' must be a character vector or scalar string.');
end
estimation = char(estimation);

resultPath = fullfile('..', estimation, 'results_Danio_rerio.mat');
if ~exist(resultPath, 'file')
    error('No results file found for estimation ''%s'' at ''%s''.', ...
        estimation, resultPath);
end

loaded = load(resultPath, 'data', 'prdData', 'txtData');
if ~isfield(loaded, 'data') || ~isfield(loaded, 'prdData')
    error('File ''%s'' must contain both ''data'' and ''prdData''.', resultPath);
end

allDatasets = fieldnames(loaded.data);
allDatasets(strcmp(allDatasets, 'psd')) = [];

if isempty(datasets)
    selectedDatasets = allDatasets;
else
    if isstring(datasets)
        datasets = cellstr(datasets);
    end
    selectedDatasets = datasets(:)';
    selectedDatasets(strcmp(selectedDatasets, 'psd')) = [];
end

reTable = compileEstimationRelativeErrors({estimation}, selectedDatasets);

nRows = numel(selectedDatasets);
dataVals = nan(nRows, 1);
predVals = nan(nRows, 1);
reVals = nan(nRows, 1);
unitsVals = strings(nRows, 1);
sourceVals = strings(nRows, 1);
keep = true(nRows, 1);

for i = 1:nRows
    dName = selectedDatasets{i};

    hasData = isfield(loaded.data, dName);
    hasPred = isfield(loaded.prdData, dName);

    if ~hasData || ~hasPred
        warning('Dataset ''%s'' missing in data or prdData for ''%s''.', ...
            dName, estimation);
        if skipMultivariate
            keep(i) = false;
        end
        continue
    end

    dVal = loaded.data.(dName);
    pVal = loaded.prdData.(dName);

    isZeroVariate = isnumeric(dVal) && isscalar(dVal) && ...
        isnumeric(pVal) && isscalar(pVal);

    if isZeroVariate
        dataVals(i) = dVal;
        predVals(i) = pVal;
    elseif skipMultivariate
        keep(i) = false;
    end

    if any(strcmp(reTable.Properties.RowNames, dName))
        reVals(i) = reTable{dName, 1};
    end

    if addUnits
        unitsVals(i) = extractTxtDataField(loaded, 'units', dName);
    end
    if addSources
        sourceVals(i) = extractTxtDataField(loaded, 'bibkey', dName);
    end
end

selectedDatasets = selectedDatasets(keep);
dataVals = dataVals(keep);
predVals = predVals(keep);
reVals = reVals(keep);
unitsVals = unitsVals(keep);
sourceVals = sourceVals(keep);

outTable = table('Size', [numel(selectedDatasets), 0], 'RowNames', selectedDatasets);
outTable.Data = dataVals;
outTable.Prediction = predVals;
outTable.RelativeError = reVals;
if addUnits
    outTable.Units = unitsVals;
end
if addSources
    outTable.Source = sourceVals;
end

end

function value = extractTxtDataField(loaded, containerName, datasetName)
value = "";
if ~isfield(loaded, 'txtData') || ~isstruct(loaded.txtData)
    return
end
if ~isfield(loaded.txtData, containerName)
    return
end

container = loaded.txtData.(containerName);
if ~isstruct(container) || ~isfield(container, datasetName)
    return
end

raw = container.(datasetName);
if ischar(raw) || (isstring(raw) && isscalar(raw))
    value = string(raw);
elseif iscell(raw)
    value = strjoin(string(raw(:))', ', ');
else
    value = string(raw);
end
end
