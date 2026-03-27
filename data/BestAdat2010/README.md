# BestAdat2010

## What this study folder stores
This folder stores the larval zebrafish growth trajectory from Best et al.
(2010).

The exposed dataset id is:
- `tL_BestAdat2010`: age versus total length

## Files and columns
`raw/tL_BestAdat2010_raw.csv` contains:
- `age_d`: age in days post fertilization
- `total_length_cm`: total length in cm

`derived/tL_BestAdat2010.csv` contains the same columns and values in the
same units:
- `age_d`
- `total_length_cm`

`dataset_manifest.json` maps the derived file into the calibration loader
and provides the label, units, temperature, comment, and bibkey metadata.

## Transformation
No transformation is required for this study.

The raw and derived files are intentionally identical because:
- the digitized data are already in the units used by the calibration
- keeping a separate derived file preserves the same external data contract
  as transformed studies

There is therefore no `transform/` script for this study.

## Notes and assumptions
- The legacy `mydata` notes describe these points as digitized larval
  growth data from Best and Adatto 2010.
