function generate_derived_data()
%GENERATE_DERIVED_DATA Convert BagaPels2001 raw measurements to calibration units.

studyRoot = fileparts(fileparts(mfilename('fullpath')));
rawPath = fullfile(studyRoot, 'raw', 'tLWWY_raw.csv');
derivedPath = fullfile(studyRoot, 'derived', 'tLWWY_converted.csv');

rawTable = readtable(rawPath);
derivedTable = table();
derivedTable.age_d = rawTable.age_d;
derivedTable.length_cm = rawTable.length_mm / 10;
derivedTable.wet_weight_g = rawTable.wet_weight_mg / 1000;
derivedTable.dry_weight_g = rawTable.dry_weight_mg / 1000;
derivedTable.yolk_cm3 = rawTable.yolk_mm3 / 1000;

writetable(derivedTable, derivedPath);
end
