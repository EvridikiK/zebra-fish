function data = read_json_file(jsonPath)
%READ_JSON_FILE Decode a JSON file into a MATLAB struct.

fid = fopen(jsonPath, 'r');
if fid == -1
    error('data_pipeline:FileOpenFailed', 'Unable to open JSON file: %s', jsonPath);
end

cleanup = onCleanup(@() fclose(fid));
raw = fread(fid, '*char')';
data = jsondecode(raw);
end
