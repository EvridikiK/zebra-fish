function generate_derived_data()
%GENERATE_DERIVED_DATA Regenerate validation-ready Lucas2014 tables.

studyRoot = fileparts(fileparts(mfilename('fullpath')));
rawDir = fullfile(studyRoot, 'raw');
derivedDir = fullfile(studyRoot, 'derived');

if ~exist(derivedDir, 'dir')
    mkdir(derivedDir);
end

larvae = readtable(fullfile(rawDir, '5_day_larvae.csv'));
juveniles = readtable(fullfile(rawDir, '2_month_juveniles.csv'));
adults = readtable(fullfile(rawDir, '6_month_adults.csv'));

oxygenByStage = [ ...
    repack_stage_table(larvae, "larvae", 5); ...
    repack_stage_table(juveniles, "juvenile", 60); ...
    repack_stage_table(adults, "adult", 180)];
writetable(oxygenByStage, fullfile(derivedDir, 'oxygen_by_stage.csv'));

summaryRaw = readtable(fullfile(rawDir, 'age_length_weight_raw.csv'));
summaryTable = table(summaryRaw.age, summaryRaw.weight, summaryRaw.se_weight, ...
    summaryRaw.length, summaryRaw.se_length, ...
    'VariableNames', {'age_d', 'weight_g', 'se_weight_g', 'length_mm', 'se_length_mm'});
writetable(summaryTable, fullfile(derivedDir, 'age_length_weight_summary.csv'));
end

function stageTable = repack_stage_table(tbl, stageName, ageDays)
stageTable = table( ...
    repmat(string(stageName), height(tbl), 1), ...
    repmat(ageDays, height(tbl), 1), ...
    tbl.weight, ...
    tbl.oxygen_consumption, ...
    'VariableNames', {'stage', 'age_d', 'weight_g', 'oxygen_consumption_mg_o2_h'});
end
