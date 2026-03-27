function generate_derived_data()
%GENERATE_DERIVED_DATA Convert Eaton and Farley length data from mm to cm.

studyRoot = fileparts(fileparts(mfilename('fullpath')));
rawPath = fullfile(studyRoot, 'raw', 'tL_EatoFarl1974_raw.csv');
derivedDir = fullfile(studyRoot, 'derived');
derivedPath = fullfile(derivedDir, 'tL_EatoFarl1974.csv');

if ~exist(derivedDir, 'dir')
    mkdir(derivedDir);
end

rawTable = readtable(rawPath);
derivedTable = table();
derivedTable.age_d = rawTable.age_d;
derivedTable.standard_length_cm = rawTable.standard_length_mm / 10;

writetable(derivedTable, derivedPath);
end
