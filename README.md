# Zebrafish DEB Parameter Estimation


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19481412.svg)](https://doi.org/10.5281/zenodo.19481412)


This repository develops Dynamic Energy Budget (DEB) models for zebrafish
(`Danio rerio`) under different levels of data availability. The main goal is
to compare sparse and rich calibrations while keeping a shared mechanistic DEB
formulation and a shared external-data pipeline.

The two main scientific products are:

- `data_2p5/`: a sparse calibration built around a low-completeness data set
  (`metaData.COMPLETE = 2.5`)
- `data_rich/`: the richest calibration in the repository, with multiple
  growth, feeding, respiration, and reproduction data sets
  (`metaData.COMPLETE = 5`)

The repository also contains:

- a centralized external data layer in `data/`
- shared MATLAB data-loading helpers in `data_pipeline/`
- shared DEB model equations in `ode/`
- validation and comparison scripts in `validation/`
- the original AmP calibration in `AmPentry/`, kept for reference only

## Citation

If you use this repository, please cite the version archived on Zenodo:
```bibtex
@software{klagkou2026zebrafish,
  title = {Zebrafish DEB Parameter Estimation},
  author = {Klagkou, Evridiki and Tan, Tjui Yeuw and Lagunes, María-José and Oliveira, Diogo F.},
  year = 2026,
  version = {1.0.0},
  doi = {10.5281/zenodo.19481412},
  url = {https://zenodo.org/records/19481412}
}
```

If you use this code, please also cite the associated manuscript:

Oliveira, D. F., Sulc, A., Le Moan, E., Donati, E., Klagkou, E., Lagunes, M.-J., Dajčman, U., and Tan, T. Y. (2026).
*How to build Dynamic Energy Budget models: a step-by-step guide*.  
Submitted to *Ecological Modelling*.

> Status: submitted, under review. This citation will be updated with full bibliographic details once available.
  

## Repository overview

### Top-level structure

| Path | Role |
| --- | --- |
| `data/` | Central source of truth for calibration and validation input data |
| `data_pipeline/` | Shared MATLAB loader infrastructure for reading data from `data/` |
| `data_2p5/` | Sparse zebrafish calibration |
| `data_rich/` | Rich zebrafish calibration |
| `ode/` | Shared mechanistic DEB equations and helper functions |
| `validation/` | Post-processing, figure generation, and comparison utilities |
| `AmPentry/` | Original AmP zebrafish entry kept as a read-only reference |
| `calibration_data_summary.md` | Human-readable summary of calibration data coverage |
| `CITATION.cff` | Repository citation metadata used by GitHub and release archiving workflows |
| `AGENTS.md` | Repository-specific working notes and conventions for coding agents |

### Core naming convention

The standard calibration entry points always keep the species suffix:

- `mydata_Danio_rerio.m`
- `pars_init_Danio_rerio.m`
- `predict_Danio_rerio.m`
- `run_Danio_rerio.m`

Even when discussed informally as `mydata`, `pars_init`, `predict`, or `run`,
the actual files in this repository keep the full `Danio_rerio` suffix.

## Folder-by-folder guide

### `data/`

`data/` is the centralized external-data layer used by the project
calibrations and validation scripts.

The layout is study-centric:

- each study has its own folder, such as `data/BestAdat2010/` or
  `data/Lucas2014/`
- each study folder typically contains:
  - `raw/`: raw or digitized source tables
  - `derived/`: calibration-ready or validation-ready tables actually consumed
    by code
  - `dataset_manifest.json`: machine-readable dataset definitions
  - `README.md`: human-readable notes on the study and transformations
  - `transform/`: optional MATLAB scripts used to regenerate derived tables
- `data/dataset_registry.json` is the shared entry point used by the loader to
  discover all study manifests
- `data/references.bib` stores the BibTeX entries referenced by dataset
  manifests and calibration metadata

Important contract:

- dataset manifests describe datasets themselves, not which calibration uses
  them
- calibration usage is declared only in each calibration's
  `mydata_Danio_rerio.m`
- bibkeys are stored in manifests and later assembled into `metaData.biblist`

### `data_pipeline/`

This folder contains the shared MATLAB infrastructure for reading external data
from `data/` and converting it into DEBtool-compatible structs.

Key files:

- `load_calibration_data.m`
  - main shared entry point for calibration data loading
  - reads `data/dataset_registry.json`
  - loads study manifests
  - reads derived CSV files
  - assembles `data`, `auxData`, and the externally sourced portions of
    `txtData`
  - attaches units, labels, comments, temperatures, optional `init`, optional
    `treat`, and bibkeys from manifest metadata
- `build_biblist_from_bibkeys.m`
  - collects the bibkeys referenced by a calibration and builds
    `metaData.biblist` from `data/references.bib`
- `get_repo_root.m`
  - resolves the repository root
- `read_json_file.m`
  - reads dataset manifests and the dataset registry
- `README.md`
  - brief local summary of the data-pipeline responsibilities

### `data_2p5/`

This is the sparse calibration representing limited data availability.

Key files:

- `mydata_Danio_rerio.m`
  - declares the zero-variate and univariate dataset IDs used by the sparse
    calibration
  - calls `data_pipeline/load_calibration_data.m`
  - builds metadata, weights, pseudodata, plot titles, and bibliography
- `pars_init_Danio_rerio.m`
  - defines the `abj` DEB parameterization
  - provides initial parameter values
  - marks which parameters are free during estimation
  - includes one dataset-specific food level parameter
    (`f_BestAdat2010`)
- `predict_Danio_rerio.m`
  - computes zero-variate predictions such as life-history traits, ultimate
    size, and reproduction rate
  - simulates the `BestAdat2010` larval growth trajectory using `ode45` and
    the shared `ode_VEHRsMG` model
- `run_Danio_rerio.m`
  - configures DEBtool/AmP estimation options
  - runs repeated Nelder-Mead calibrations with alternating simplex sign until
    improvement falls below tolerance or the run limit is reached
  - writes results, figures, and HTML output
  - reloads the fitted parameters, recomputes predictions, and saves
    `prdData` back into `results_Danio_rerio.mat`

Typical outputs inside this folder:

- `results_Danio_rerio.mat`
- `Danio_rerio_res.html`
- `results_Danio_rerio_*.png`

### `data_rich/`

This is the main rich calibration used for the broadest data coverage in the
repository.

Key files:

- `mydata_Danio_rerio.m`
  - declares the full set of study-specific dataset IDs used in the rich
    calibration
  - loads external data through `data_pipeline/load_calibration_data.m`
  - assembles `auxData.temp`, `auxData.treat`, `auxData.init`, weights,
    titles, and bibliography metadata
- `pars_init_Danio_rerio.m`
  - defines the full `abj` parameter set
  - frees a broader set of core DEB parameters
  - includes multiple dataset-specific feeding-level parameters for the
    different experiments
- `predict_Danio_rerio.m`
  - serves as the main forward model for the rich calibration
  - applies a custom parameter filter
  - computes zero-variate predictions
  - simulates multiple growth, feeding, oxygen-consumption, and reproduction
    data sets
  - derives secondary outputs such as oxygen consumption and gonado-somatic
    index
- `run_Danio_rerio.m`
  - follows the same overall estimation loop as `data_2p5/`
  - uses `estim_options('filter', 1)`

Typical outputs inside this folder:

- `results_Danio_rerio.mat`
- `Danio_rerio_res.html`
- `results_Danio_rerio_*.png`

### `ode/`

This folder contains the shared mechanistic core used by the project
calibrations and several validation scripts.

Key files:

- `ode_VEHRsMG.m`
  - core ODE for the state vector
    `[V, E, E_H, E_R, s_M, egg_buffer]`
- `compute_powers.m`
  - computes DEB power fluxes such as assimilation, mobilization, maintenance,
    growth, and reproduction-related terms
- `compute_dissipation_power.m`
  - computes dissipative power, used for oxygen-consumption calculations
- `getAgeAndLengthAtTransitions.m`
  - solves the ODE with event detection to obtain ages and lengths at birth,
    metamorphosis, and puberty
- `getE0.m`
  - computes initial egg reserve from DEB compound parameters
- `get_max_E_R.m`
  - estimates the maximum reproduction buffer for a given state and food level

### `validation/`

This folder contains post-processing scripts used to compare calibrations,
summarize results, and produce manuscript-ready figures.

Key files:

- `compileEstimationParameters.m`
  - loads `results_Danio_rerio.mat` from one or more calibration folders
  - compiles selected fitted parameters into a comparison table
- `compileEstimationRelativeErrors.m`
  - compiles dataset-level relative errors and mean relative error (`MRE`)
    across calibrations
- `compileEstimationDataPredictionRelativeErrors.m`
  - builds a table of observed values, predicted values, and relative errors
    for one calibration
  - can also include units and bibkey-based source information
- `lawrence_et_al_2008_visualization.m`
  - loads `data_rich/results_Danio_rerio.mat`
  - reads observed data from `data/LawrEber2008/derived/`
  - generates high-food and low-food growth comparison figures
- `lucas_et_al_2014_validation.m`
  - compares selected calibration outputs against external Lucas et al. 2014
    oxygen, weight, and length data
  - reads observed data from `data/Lucas2014/derived/`
  - produces four figure families:
    oxygen vs weight, oxygen vs age, weight vs age, and length vs age
- `table2latex.m`
  - converts MATLAB tables into simple LaTeX tabular output
- `validate_biblist_refactor.m`
  - checks that dataset manifests and calibration biblists are consistent with
    `data/references.bib`

The `validation/figures/` subfolder stores exported `.png` and `.pdf` figures.

### `AmPentry/`

This folder contains the original AmP zebrafish calibration kept for
reference.

It includes its own:

- `mydata_Danio_rerio.m`
- `pars_init_Danio_rerio.m`
- `predict_Danio_rerio.m`
- `run_Danio_rerio.m`

This folder is not the main project product and should be treated as
read-only unless there is explicit reason to modify it.
The code in `AmPentry/` follows the original AmP style and does not use the
shared data pipeline or the shared ODE layer developed in this repository.
It is only kept here for reference and comparison purposes.

Citation for the original AmP entry:

```bibtex
@misc{augustine_amp_2018,
  title = {{{AmP Danio}} Rerio, Version 2018/08/09},
  shorttitle = {{{AmP Danio}} Rerio},
  author = {Augustine, Starrlight},
  year = 2018,
  month = aug,
  journal = {Add-my-Pet (AmP) collection},
  urldate = {2025-07-17},
  howpublished = {https://www.bio.vu.nl/thb/deb/deblab/add\_my\_pet/entries\_web/Danio\_rerio/Danio\_rerio\_res.html}
}
```

## Calibration workflow

The standard calibration flow in `data_2p5/` and `data_rich/` is:

1. `run_Danio_rerio.m` configures AmP/DEBtool estimation options and runs
   `estim_pars`.
2. `estim_pars` calls:
   - `mydata_Danio_rerio.m`
   - `pars_init_Danio_rerio.m`
   - `predict_Danio_rerio.m`
3. `mydata_Danio_rerio.m` selects the dataset IDs used by that calibration and
   calls `data_pipeline/load_calibration_data.m`.
4. The shared loader reads the registry and manifests in `data/`, loads the
   appropriate derived tables, and assembles `data`, `auxData`, and
   calibration text metadata.
5. `pars_init_Danio_rerio.m` defines parameter starting values and the
   free/fixed status of parameters.
6. `predict_Danio_rerio.m` evaluates the forward DEB model, calling shared
   mechanistic functions in `ode/`.
7. The calibration writes outputs such as:
   - `results_Danio_rerio.mat`
   - `Danio_rerio_res.html`
   - result figures
8. The run script then reloads the fitted result, recomputes `prdData`, and
   saves the full working bundle back into `results_Danio_rerio.mat`.

That `.mat` file is the main downstream bundle used by the repository. It
typically stores at least:

- fitted parameters
- metadata
- observed data
- weights
- predictions (`prdData`)

## Validation and comparison flow

Validation scripts use two main information sources:

- saved calibration outputs from `results_Danio_rerio.mat`
- centralized external data in `data/`

In practice, the flow is:

1. Run one or more calibrations in `data_2p5/` or `data_rich/`.
2. Use scripts in `validation/` to load those saved `.mat` outputs.
3. Optionally combine them with external validation tables from `data/`.
4. Generate comparison tables, relative-error summaries, and manuscript-ready
   figures under `validation/figures/`.

## Data architecture

The current repository direction is to keep all externally sourced calibration
and validation data under `data/` and keep the calibration folders focused on:

- selecting which datasets to use
- defining parameter sets
- evaluating the model
- storing results

This means:

- new studies should be added under `data/`, not embedded directly into
  calibration scripts
- shared CSV/JSON parsing logic should live in `data_pipeline/`
- shared DEB equations should live in `ode/`
- the calibration `mydata_Danio_rerio.m` files should stay readable and
  study-oriented, even though the numeric payload now lives externally

## Common outputs

Within each calibration folder, expect:

- `results_Danio_rerio.mat`: saved parameter, data, and prediction bundle
- `Danio_rerio_res.html`: HTML summary produced by the AmP estimation workflow
- `results_Danio_rerio_*.png`: calibration figures

Within `validation/figures/`, expect:

- figure exports in `.png` and `.pdf` format used for validation and
  manuscript support

## External dependencies

These MATLAB scripts assume the AmP/DEBtool ecosystem is already available on
the MATLAB path. The repository depends on functions such as:

- `check_my_pet`
- `estim_options`
- `estim_pars`
- `parscomp_st`
- `addchem`
- `addpseudodata`
- `setweights`
- `tempcorr`
- `get_tm_mod`
- `vars_pull`
- `initial_scaled_reserve`

Those functions are not defined in this repository.

## Typical usage

From within a calibration folder such as `data_rich/` or `data_2p5/`, the
usual entry point is:

```matlab
run_Danio_rerio
```

From within `validation/`, typical entry points are scripts such as:

```matlab
lucas_et_al_2014_validation
lawrence_et_al_2008_visualization
```

or table-building utilities such as:

```matlab
compileEstimationParameters({'data_rich','data_2p5'})
compileEstimationRelativeErrors({'data_rich','data_2p5'})
```
