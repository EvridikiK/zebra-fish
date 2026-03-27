function generate_derived_data()
%GENERATE_DERIVED_DATA Convert the Schi2002 growth curve from mm to cm.

studyRoot = fileparts(fileparts(mfilename('fullpath')));
rawPath = fullfile(studyRoot, 'raw', 'tL_Schi2002_raw.csv');
derivedDir = fullfile(studyRoot, 'derived');
derivedPath = fullfile(derivedDir, 'tL_Schi2002.csv');

if ~exist(derivedDir, 'dir')
    mkdir(derivedDir);
end

rawTable = readtable(rawPath);
derivedTable = table();
derivedTable.time_d = rawTable.time_d;
derivedTable.total_length_cm = rawTable.total_length_mm / 10;

writetable(derivedTable, derivedPath);
end
