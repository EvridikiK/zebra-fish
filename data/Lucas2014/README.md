# Lucas2014

## What this study folder stores
This folder stores the zebrafish oxygen-consumption, weight, and length
data digitized from Lucas et al. (2014).

The exposed dataset ids are:
- `oxygen_by_stage`: individual oxygen-consumption observations with body mass and age
- `age_length_weight_summary`: age-level summary values for weight and length, with standard errors

## Files and columns
The raw folder contains:
- `5_day_larvae.csv`: individual body mass and oxygen-consumption points for 5 dpf larvae
- `2_month_juveniles.csv`: individual body mass and oxygen-consumption points for 60 dpf juveniles
- `6_month_adults.csv`: individual body mass and oxygen-consumption points for 180 dpf adults
- `age_length_weight_raw.csv`: age-level mean weight and length values with standard errors

The derived folder contains:
- `oxygen_by_stage.csv`: combined oxygen-consumption table with `stage`, `age_d`, `weight_g`, and `oxygen_consumption_mg_o2_h`
- `age_length_weight_summary.csv`: renamed summary table with `age_d`, `weight_g`, `se_weight_g`, `length_mm`, and `se_length_mm`

`dataset_manifest.json` documents the derived files, column meanings, units,
and source notes for the validation script.

## Transformation
The MATLAB script `transform/generate_derived_data.m` performs two steps:
- combine the three oxygen-consumption raw tables into one derived table with explicit `stage` and `age_d` columns
- rename the summary-table columns into a consistent loader-friendly format

## Notes and assumptions
- These data were digitized from Figure 1(b) and the accompanying summary
  values reported in Lucas et al. (2014).
- The oxygen-consumption values are stored in `mg O_2/h`, matching the
  validation plotting script.
