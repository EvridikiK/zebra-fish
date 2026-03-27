# LawrEber2008

## What this study folder stores
This folder stores zebrafish growth data from Lawrence, Ebersole, and
co-authors under two food regimes.

The exposed dataset ids are:
- `tL_LawrEber2002_high`: age versus total length under high food
- `tL_LawrEber2002_low`: age versus total length under low food
- `tL_LawrEber2008_high`: age versus total length under high food
- `tL_LawrEber2008_low`: age versus total length under low food

## Files and columns
`raw/tL_conditions_raw.csv` contains both treatments in one table:
- `condition`: treatment label, `high` or `low`
- `time_d`: time in days
- `total_length_cm`: total length in cm

`derived/tL_conditions.csv` contains the same columns and values:
- `condition`
- `time_d`
- `total_length_cm`

`dataset_manifest.json` uses row filters on `condition` to expose the high-
and low-food series as separate dataset ids while keeping a single
underlying study table.

## Transformation
No transformation is required for this study.

The raw and derived files are identical because:
- the source table is already in calibration units
- the only dataset-specific logic is row selection by treatment, which is
  handled declaratively in the manifest rather than by a MATLAB script

There is therefore no `transform/` script for this study.

## Notes and assumptions
- The legacy dataset ids `tL_LawrEber2002_high` and
  `tL_LawrEber2002_low` are retained for compatibility.
- Labels and temperatures are defined in the manifest to preserve current
  `mydata` behavior.
