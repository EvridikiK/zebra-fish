# BeauGous2015

## What this study folder stores
This folder stores the juvenile growth, final size, final wet weight, and
reproduction summaries derived from Beaugrand, Goussen, and co-authors.

The exposed dataset ids are:
- `tLf1_BeauGous2015`
- `tLf2_BeauGous2015`
- `tLf3_BeauGous2015`
- `tL1`
- `Wwt`
- `tN`

## Files and columns
The raw folder contains:
- `tLf1_raw.csv`, `tLf2_raw.csv`, `tLf3_raw.csv`: age in dpf and standard length in mm
- `L_BeauGous2015_raw.csv`: initial and final lengths in mm
- `Ww_raw.csv`: final wet weights in mg
- `tN_raw.csv`: daily egg counts per female

The derived folder contains calibration-ready files:
- `tLf1_BeauGous2015.csv`, `tLf2_BeauGous2015.csv`, `tLf3_BeauGous2015.csv`
- `tL1.csv`: two-point mean length trajectory in cm
- `Wwt.csv`: mean final wet weight in g
- `tN.csv`: mean cumulative egg production over time

## Transformation
The MATLAB transform script performs the legacy `mydata` transformations:
- copy the three juvenile growth tables as-is
- average the initial and final length observations and convert mm to cm for `tL1`
- average final wet weights and convert mg to g for `Wwt`
- compute the mean cumulative egg production curve for `tN`

## Notes and assumptions
- Comments and treatment metadata remain in the manifest to preserve the
  rich calibration outputs.
- `tN` carries an init value equal to the first length in `tL1`, matching
  the legacy implementation.
