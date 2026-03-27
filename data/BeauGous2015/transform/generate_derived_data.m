function generate_derived_data()
%GENERATE_DERIVED_DATA Regenerate BeauGous2015 derived calibration inputs.

studyRoot = fileparts(fileparts(mfilename('fullpath')));
derivedDir = fullfile(studyRoot, 'derived');
if ~exist(derivedDir, 'dir')
    mkdir(derivedDir);
end

tLf1 = readtable(fullfile(studyRoot, 'raw', 'tLf1_raw.csv'));
tLf2 = readtable(fullfile(studyRoot, 'raw', 'tLf2_raw.csv'));
tLf3 = readtable(fullfile(studyRoot, 'raw', 'tLf3_raw.csv'));
writetable(tLf1, fullfile(derivedDir, 'tLf1_BeauGous2015.csv'));
writetable(tLf2, fullfile(derivedDir, 'tLf2_BeauGous2015.csv'));
writetable(tLf3, fullfile(derivedDir, 'tLf3_BeauGous2015.csv'));

lengthTable = readtable(fullfile(studyRoot, 'raw', 'L_BeauGous2015_raw.csv'));
tL1 = table([0; 19], ...
    [mean(lengthTable.initial_length_mm); mean(lengthTable.final_length_mm)] / 10, ...
    'VariableNames', {'time_d', 'standard_length_cm'});
writetable(tL1, fullfile(derivedDir, 'tL1.csv'));

weightTable = readtable(fullfile(studyRoot, 'raw', 'Ww_raw.csv'));
Wwt = table(mean(weightTable.wet_weight_mg) * 1e-3, 'VariableNames', {'value_g'});
writetable(Wwt, fullfile(derivedDir, 'Wwt.csv'));

tNRaw = readtable(fullfile(studyRoot, 'raw', 'tN_raw.csv'));
eggMatrix = tNRaw{:, 2:end};
tN = table(tNRaw.time_d, mean(cumsum(eggMatrix, 1), 2), ...
    'VariableNames', {'time_d', 'cumulated_eggs_n'});
writetable(tN, fullfile(derivedDir, 'tN.csv'));
end
