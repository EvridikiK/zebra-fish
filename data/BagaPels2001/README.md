# BagaPels2001

## What this study folder stores
This folder stores the larval zebrafish growth and weight data from
Bagatto, Pelster, and Burggren (2001).

The exposed dataset ids are:
- `tL_BagaPels2001`: age versus standard length
- `tWw_BagaPels2001`: age versus wet weight
- `tWd_BagaPels2001`: age versus dry weight

## Files and columns
`raw/tLWWY_raw.csv` stores the digitized data in the original units:
- `age_d`: age in days
- `length_mm`: length in mm
- `wet_weight_mg`: wet weight in mg
- `dry_weight_mg`: dry weight in mg
- `yolk_mm3`: yolk volume in mm^3

`derived/tLWWY_converted.csv` stores the calibration-ready version:
- `age_d`: age in days
- `length_cm`: length converted to cm
- `wet_weight_g`: wet weight converted to g
- `dry_weight_g`: dry weight converted to g
- `yolk_cm3`: yolk volume converted to cm^3

`dataset_manifest.json` maps the study table into dataset ids and defines
labels, units, temperature, and bibkey metadata.

## Transformation
The MATLAB script `transform/generate_derived_data.m` converts the raw
measurements into calibration units.

Transformation steps:
- divide `length_mm` by 10 to obtain `length_cm`
- divide `wet_weight_mg`, `dry_weight_mg`, and `yolk_mm3` by 1000 to obtain
  gram and cm^3 units
- expose separate time-length, time-wet-weight, and time-dry-weight
  datasets from the same derived table

## Notes and assumptions
- Dry weight and yolk are retained in the derived file because they are
  part of the same source table.
- The label text keeps the legacy spelling `fertilizatixon` to preserve
  current output behavior.
