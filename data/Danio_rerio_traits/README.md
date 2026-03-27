# Danio_rerio_traits

## What this study folder stores
This folder stores the canonical zero-variate scalar traits for zebrafish.
These are the scalar values that were previously embedded directly in
`mydata_Danio_rerio.m` files.

The exposed dataset ids are:
- `ab`: age at birth
- `aj`: age at metamorphosis
- `ap`: age at puberty
- `am`: life span
- `L0`: egg diameter
- `Lb`: total length at birth
- `Lj`: total length at metamorphosis
- `Lp`: standard length at puberty
- `Li`: ultimate total length
- `Wd0`: egg dry weight
- `Wwi`: ultimate wet weight
- `Ri`: maximum reproduction rate
- `GSI`: gonado-somatic index

## Files and columns
`raw/traits_raw.csv` contains one row per scalar trait with these columns:
- `dataset_id`: stable id used by the calibration loader
- `value`: observed scalar value
- `unit`: human-readable data unit
- `label`: human-readable meaning of the trait
- `temp_celsius`: measurement temperature in Celsius when relevant
- `comment`: interpretation note carried over from the legacy `mydata`

`derived/traits.csv` contains the calibration-ready values with:
- `dataset_id`
- `value`

`dataset_manifest.json` maps each trait id to:
- the derived file and row selector
- units and labels used to populate `txtData`
- bibkeys used to populate `txtData.bibkey`
- temperatures used to populate `auxData.temp`
- comments used to populate `txtData.comment`
- optional `init` metadata, currently used for `GSI`

## Transformation
The MATLAB script `transform/generate_derived_data.m` produces the derived
table from the richer raw table.

Transformation steps:
- keep one row per scalar trait
- retain only the columns needed at load time in the derived file
- move units, labels, comments, and temperatures to the JSON manifest

## Notes and assumptions
- This folder is the single source of truth for zebrafish scalar traits.
- The values and labels follow the current canonical dataset definition.
