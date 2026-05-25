# DRC malaria ITN simulation: resistance, PBO nets and WHO-calibrated burden

This repository script runs the DRC malaria intervention analysis used to compare:

1. no ITNs,
2. ITNs under observed pyrethroid resistance,
3. ITNs without resistance,
4. pyrethroid-only ITNs versus the observed pyrethroid-PBO ITN mix.

The workflow produces province-level and national summaries for cases, deaths, lives saved, PBO impact, and resistance-attributable burden.

## Main script

```r
Rscript malaria_drc_three_scenario_analysis_github.R
```

## Recommended repository structure

```text
.
├── malaria_drc_three_scenario_analysis_github.R
├── data/
│   └── raw/
│       ├── DRC_Data.xlsx
│       ├── bednet_distribution_by_province_2008_2022.xlsx
│       ├── Fully_Corrected_ITN_Usage_Coverage_Data.csv
│       ├── gamma_p_first100_samples.csv
│       ├── itn_efficacy_posterior_draws_100samples.csv
│       ├── itn_efficacy_pbo_posterior_draws_100samples.csv
│       ├── type_itn_pbo.xlsx
│       └── who_new_data_2024.xlsx
└── outputs/
    └── tables/
```

The input data are usually large and may be confidential, so they should normally not be committed to GitHub. Use `data/raw/` locally, or set environment variables to point to your files.

## Required R packages

The script checks for these packages before running:

```r
malariasimulation
malariaEquilibrium
ggplot2
readxl
readr
dplyr
tidyr
future
furrr
janitor
stringr
purrr
tibble
stringi
scales
```

## Configuring paths

The script first checks environment variables, then project-relative paths under `data/raw/`, then the original Windows paths used during development.

Example:

```r
Sys.setenv(
  DRC_DATA_XLSX = "data/raw/DRC_Data.xlsx",
  WHO_BURDEN_XLSX = "data/raw/who_new_data_2024.xlsx",
  N_DRAWS = "100",
  N_WORKERS = "8"
)
source("malaria_drc_three_scenario_analysis_github.R")
```

## WHO calibration convention

The observed ITN-with-resistance scenario is treated as the status-quo scenario. WHO matching is applied only to resistance-attributable case/death outputs. Other ITN and PBO comparisons are kept on the raw population-weighted model scale.

## Outputs

Most summary tables are written to:

```text
outputs/tables/
```

Key outputs include:

- national yearly ITN impact tables,
- province-level ITN impact tables,
- PBO additional cases/deaths averted tables,
- resistance-attributable burden summaries,
- manuscript-ready printed sentences for replacing placeholder values.
