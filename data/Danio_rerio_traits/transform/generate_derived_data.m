function generate_derived_data()
%GENERATE_DERIVED_DATA Regenerate calibration-ready canonical scalar traits.

studyRoot = fileparts(fileparts(mfilename('fullpath')));
rawPath = fullfile(studyRoot, 'raw', 'traits_raw.csv');
derivedPath = fullfile(studyRoot, 'derived', 'traits.csv');

rawTable = readtable(rawPath);
derivedTable = rawTable(:, {'dataset_id', 'value'});

writetable(derivedTable, derivedPath);
end
