# data_pipeline

This folder contains the shared MATLAB infrastructure for loading external
calibration and validation data from `data/`.

Current responsibilities:
- resolve the repository root: `get_repo_root.m`
- read JSON manifests: `read_json_file.m`
- assemble DEBtool-compatible `data`, `auxData`, and `txtData` structs: `load_calibration_data.m`
