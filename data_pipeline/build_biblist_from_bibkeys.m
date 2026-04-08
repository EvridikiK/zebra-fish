function biblist = build_biblist_from_bibkeys(bibkeyStruct)
%BUILD_BIBLIST_FROM_BIBKEYS Build metaData.biblist from data/references.bib.
%   biblist = BUILD_BIBLIST_FROM_BIBKEYS(bibkeyStruct) collects the unique
%   BibTeX keys referenced in the txtData.bibkey-style struct bibkeyStruct,
%   loads data/references.bib, and returns a struct whose fields are those
%   keys and whose values are the raw BibTeX entries.

if nargin < 1 || ~isstruct(bibkeyStruct)
    error('data_pipeline:InvalidBibkeyStruct', ...
        'Input must be a txtData.bibkey-style struct.');
end

requestedKeys = collect_unique_bibkeys(bibkeyStruct);
entryLookup = parse_bibtex_entries(fullfile(get_repo_root(), 'data', 'references.bib'));

biblist = struct();
for i = 1:numel(requestedKeys)
    bibkey = requestedKeys{i};
    if ~isKey(entryLookup, bibkey)
        error('data_pipeline:MissingBibEntry', ...
            'BibTeX key "%s" was requested but is missing from data/references.bib.', ...
            bibkey);
    end
    biblist.(bibkey) = entryLookup(bibkey);
end
end

function requestedKeys = collect_unique_bibkeys(bibkeyStruct)
requestedKeys = {};
seenKeys = containers.Map('KeyType', 'char', 'ValueType', 'logical');

datasetNames = fieldnames(bibkeyStruct);
for i = 1:numel(datasetNames)
    datasetName = datasetNames{i};
    rawKeys = normalize_bibkey_value(bibkeyStruct.(datasetName), datasetName);

    for j = 1:numel(rawKeys)
        bibkey = rawKeys{j};
        if ~isKey(seenKeys, bibkey)
            seenKeys(bibkey) = true;
            requestedKeys{end + 1} = bibkey; %#ok<AGROW>
        end
    end
end
end

function entryLookup = parse_bibtex_entries(bibPath)
if ~exist(bibPath, 'file')
    error('data_pipeline:BibFileMissing', ...
        'Unable to find BibTeX file: %s', bibPath);
end

rawBib = fileread(bibPath);
entryLookup = containers.Map('KeyType', 'char', 'ValueType', 'char');

cursor = 1;
textLength = numel(rawBib);
while cursor <= textLength
    atPos = find(rawBib(cursor:end) == '@', 1, 'first');
    if isempty(atPos)
        break
    end
    atPos = atPos + cursor - 1;

    [entryKey, entryText, nextCursor] = parse_single_entry(rawBib, atPos);
    if isempty(entryKey)
        cursor = atPos + 1;
        continue
    end

    if isKey(entryLookup, entryKey)
        error('data_pipeline:DuplicateBibEntry', ...
            'BibTeX key "%s" appears more than once in data/references.bib.', ...
            entryKey);
    end

    entryLookup(entryKey) = entryText;
    cursor = nextCursor;
end
end

function [entryKey, entryText, nextCursor] = parse_single_entry(rawBib, atPos)
entryKey = '';
entryText = '';
nextCursor = atPos + 1;

cursor = atPos + 1;
textLength = numel(rawBib);

while cursor <= textLength && isspace(rawBib(cursor))
    cursor = cursor + 1;
end

while cursor <= textLength && isstrprop(rawBib(cursor), 'alpha')
    cursor = cursor + 1;
end

while cursor <= textLength && isspace(rawBib(cursor))
    cursor = cursor + 1;
end

if cursor > textLength || ~any(rawBib(cursor) == ['{', '('])
    return
end

openingDelimiter = rawBib(cursor);
if openingDelimiter == '{'
    closingDelimiter = '}';
else
    closingDelimiter = ')';
end

commaOffset = find(rawBib(cursor + 1:end) == ',', 1, 'first');
if isempty(commaOffset)
    error('data_pipeline:InvalidBibEntry', ...
        'Malformed BibTeX entry starting at character %d in data/references.bib.', ...
        atPos);
end

commaPos = cursor + commaOffset;
entryKey = strtrim(rawBib(cursor + 1:commaPos - 1));
if isempty(entryKey)
    error('data_pipeline:InvalidBibEntry', ...
        'Encountered a BibTeX entry without a key in data/references.bib.');
end

depth = 0;
entryEnd = [];
for i = cursor:textLength
    currentChar = rawBib(i);
    if currentChar == openingDelimiter
        depth = depth + 1;
    elseif currentChar == closingDelimiter
        depth = depth - 1;
        if depth == 0
            entryEnd = i;
            break
        end
    end
end

if isempty(entryEnd)
    error('data_pipeline:InvalidBibEntry', ...
        'BibTeX entry "%s" is missing a closing delimiter in data/references.bib.', ...
        entryKey);
end

entryText = rawBib(atPos:entryEnd);
nextCursor = entryEnd + 1;
end

function normalizedKeys = normalize_bibkey_value(rawValue, datasetName)
if ischar(rawValue)
    normalizedKeys = {rawValue};
elseif isstring(rawValue)
    normalizedKeys = reshape(cellstr(rawValue), 1, []);
elseif iscell(rawValue)
    normalizedKeys = cell(1, numel(rawValue));
    for i = 1:numel(rawValue)
        item = rawValue{i};
        if ischar(item)
            normalizedKeys{i} = item;
        elseif isstring(item) && isscalar(item)
            normalizedKeys{i} = char(item);
        else
            error('data_pipeline:InvalidBibkeyValue', ...
                'Dataset "%s" has a bibkey entry with unsupported type.', ...
                datasetName);
        end
    end
else
    error('data_pipeline:InvalidBibkeyValue', ...
        'Dataset "%s" has a bibkey entry with unsupported type.', ...
        datasetName);
end

for i = 1:numel(normalizedKeys)
    normalizedKeys{i} = strtrim(normalizedKeys{i});
    if isempty(normalizedKeys{i})
        error('data_pipeline:InvalidBibkeyValue', ...
            'Dataset "%s" has an empty bibkey entry.', datasetName);
    end
end
end
