function [data, auxData, txtData] = load_calibration_data(datasetIds)
%LOAD_CALIBRATION_DATA Assemble calibration data structs from external files.
%   datasetIds can be a char vector, string array, cell array of names,
%   or a legacy struct array with an "id" field.

repoRoot = get_repo_root();
registry = read_json_file(fullfile(repoRoot, 'data', 'dataset_registry.json'));
lookup = build_dataset_lookup(repoRoot, registry);
datasetIds = normalize_dataset_ids(datasetIds);

data = struct();
temp = struct();
init = struct();
treat = struct();
units = struct();
label = struct();
comment = struct();
bibkey = struct();

for i = 1:numel(datasetIds)
    datasetId = datasetIds{i};
    if ~isfield(lookup, datasetId)
        error('data_pipeline:UnknownDataset', 'Dataset "%s" is not registered.', datasetId);
    end

    entry = lookup.(datasetId);

    tablePath = fullfile(entry.study_root, entry.derived_file);
    tbl = readtable(tablePath);
    tbl = apply_row_filter(tbl, entry);

    values = tbl{:, entry.data_columns};
    if strcmp(entry.type, 'scalar')
        if numel(values) ~= 1
            error('data_pipeline:ScalarShape', ...
                'Dataset "%s" must resolve to exactly one scalar value.', datasetId);
        end
        data.(datasetId) = values(1);
        unitValues = normalize_to_cellstr(entry.units);
        labelValues = normalize_to_cellstr(entry.labels);
        units.(datasetId) = unitValues{1};
        label.(datasetId) = labelValues{1};
    else
        data.(datasetId) = values;
        units.(datasetId) = normalize_to_cellstr(entry.units);
        label.(datasetId) = normalize_to_cellstr(entry.labels);
    end

    if isfield(entry, 'temperature_celsius') && ~isempty(entry.temperature_celsius)
        temp.(datasetId) = convert_temperature(entry.temperature_celsius);
        units.temp.(datasetId) = 'K';
        label.temp.(datasetId) = 'temperature';
    end

    if isfield(entry, 'init') && ~isempty(entry.init)
        init.(datasetId) = entry.init.value;
        units.init.(datasetId) = entry.init.units;
        label.init.(datasetId) = entry.init.label;
    end

    if isfield(entry, 'treat') && ~isempty(entry.treat)
        treat.(datasetId) = normalize_aux_value(entry.treat);
        units.treat.(datasetId) = entry.treat.units;
        label.treat.(datasetId) = entry.treat.label;
    end

    if isfield(entry, 'comment') && ~isempty(entry.comment)
        comment.(datasetId) = entry.comment;
    end

    if isfield(entry, 'bibkey') && ~isempty(entry.bibkey)
        bibkey.(datasetId) = normalize_text_metadata(entry.bibkey);
    else
        error('data_pipeline:MissingBibkey', ...
            'Dataset "%s" is missing a bibkey in its manifest.', datasetId);
    end
end

auxData.temp = temp;
if ~isempty(fieldnames(init))
    auxData.init = init;
end
if ~isempty(fieldnames(treat))
    auxData.treat = treat;
end

txtData.units = units;
txtData.label = label;
txtData.comment = comment;
txtData.bibkey = bibkey;
end

function lookup = build_dataset_lookup(repoRoot, registry)
lookup = struct();

for i = 1:numel(registry.studies)
    studyEntry = get_collection_item(registry.studies, i);
    manifestPath = fullfile(repoRoot, studyEntry.manifest);
    manifest = read_json_file(manifestPath);
    studyRoot = fileparts(manifestPath);

    for j = 1:numel(manifest.datasets)
        dataset = get_collection_item(manifest.datasets, j);
        dataset.study_root = studyRoot;
        if isfield(lookup, dataset.id)
            error('data_pipeline:DuplicateDatasetId', ...
                'Dataset "%s" is declared more than once across manifests.', dataset.id);
        else
            lookup.(dataset.id) = dataset;
        end
    end
end
end

function tbl = apply_row_filter(tbl, entry)
if ~isfield(entry, 'row_filter') || isempty(entry.row_filter)
    return
end

columnValues = tbl.(entry.row_filter.column);
expected = entry.row_filter.equals;

if iscell(columnValues)
    mask = strcmp(columnValues, expected);
elseif isstring(columnValues)
    mask = columnValues == string(expected);
elseif iscategorical(columnValues)
    mask = strcmp(cellstr(columnValues), expected);
else
    mask = columnValues == expected;
end

tbl = tbl(mask, :);
end

function temperatureKelvin = convert_temperature(temperatureCelsius)
if numel(temperatureCelsius) == 1
    temperatureKelvin = C2K(temperatureCelsius);
else
    temperatureKelvin = reshape(arrayfun(@C2K, temperatureCelsius(:).'), 1, []);
end
end

function item = get_collection_item(collection, idx)
if iscell(collection)
    item = collection{idx};
else
    item = collection(idx);
end
end

function values = normalize_to_cellstr(value)
if iscell(value)
    values = reshape(value, 1, []);
elseif isstring(value)
    values = reshape(cellstr(value), 1, []);
else
    values = {value};
end
end

function datasetIds = normalize_dataset_ids(datasetIds)
if ischar(datasetIds)
    datasetIds = {datasetIds};
elseif isstring(datasetIds)
    datasetIds = reshape(cellstr(datasetIds), 1, []);
elseif isstruct(datasetIds)
    ids = cell(1, numel(datasetIds));
    for i = 1:numel(datasetIds)
        if ~isfield(datasetIds(i), 'id')
            error('data_pipeline:InvalidDatasetId', ...
                'Legacy dataset structs must define an "id" field.');
        end
        ids{i} = datasetIds(i).id;
    end
    datasetIds = ids;
elseif iscell(datasetIds)
    for i = 1:numel(datasetIds)
        item = datasetIds{i};
        if ischar(item)
            datasetIds{i} = item;
        elseif isstring(item) && isscalar(item)
            datasetIds{i} = char(item);
        elseif isstruct(item) && isfield(item, 'id')
            datasetIds{i} = item.id;
        else
            error('data_pipeline:InvalidDatasetId', ...
                'Dataset ids must be names or legacy structs with an "id" field.');
        end
    end
else
    error('data_pipeline:InvalidDatasetIds', ...
        'Dataset ids must be a name, a cell array of names, a string array, or a legacy struct array.');
end
end

function value = normalize_aux_value(entry)
value = entry.value;
if isfield(entry, 'wrap_in_cell') && entry.wrap_in_cell
    value = {value};
end
end

function value = normalize_text_metadata(value)
if isstring(value)
    if isscalar(value)
        value = char(value);
    else
        value = reshape(cellstr(value), 1, []);
    end
elseif iscell(value)
    value = reshape(value, 1, []);
end
end
