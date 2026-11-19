# Influenza vaccination across districts of South Korea, 2011-2024

Analysis code for a manuscript on the district geography of age-sex-standardized adult influenza vaccination coverage from the Korea Community Health Survey and its district-level correlates (physical activity, sports facilities, and social composition).

Archived version: Zenodo, https://doi.org/10.5281/zenodo.24360577 (v1.0.0).

This code is released under the MIT license.

## Contents

This repository contains code only. The processed district-year data (analytic panel, result tables, and boundary file) are not included; they are available from the author on reasonable request. All source statistics are public (see Data sources).

| Path | What it is |
|---|---|
| `01_스크립트/` | Analysis scripts (below) |
| `04_분석/FLU/` | Model settings (`sets_FLU_std3.json`, `config_FLU_std3.json`), the model chain (`run_fh_chain.sh`), and the district adjacency graph (`graph_593.graph`) |
| `04_분석/FLU/MS/build/` | Scripts that read the result tables and write every number, table, and figure in the manuscript and supplementary material |

Folder and some column names are in Korean because the source tables are Korean; each script has a header comment describing its inputs and outputs.

## Data sources

All inputs are public, aggregated statistics. No individual-level data are included.

- Influenza vaccination rate in the past 12 months, walking, and moderate-or-higher physical activity (age-sex-standardized rates, standard errors, and respondent counts by district): Korea Community Health Survey, Korea Disease Control and Prevention Agency (KOSIS).
- Public sports facilities and village sports facilities: Ministry of Culture, Sports and Tourism.
- University-educated share, basic livelihood security recipients, single-person households, structural controls, and negative controls: Statistics Korea (KOSIS) and local indicators.

## Pipeline

Set `ST_DATA_ROOT` and `ST_GPT_DATASET` to the local copies of the public source files and the read-only district database, then run from the repository root.

| Step | Script | Output |
|---|---|---|
| 1 | `99c_extra_canonical_pool_flu.py` | District covariate pool (welfare recipient correction, education) |
| 2 | `99FL_flu_fh_prep.py` | Logit outcome and known sampling variance (Fay-Herriot), survey-operation indicators |
| 3 | `99FL0_flu_sqlite_crosscheck.py` | Cell-by-cell check of inputs against the read-only district database |
| 4 | `80_build_disease_panel_v.py FLU` (`ST_VER=std3`) | Analytic panel (imputation, standardization, VIF) |
| 5 | `81fh_flu.R FLU` via `04_분석/FLU/run_fh_chain.sh` | Main spatiotemporal Fay-Herriot model, set comparison, multiverse of all 127 covariate subsets, Knorr-Held models (`kh_constraints.R`), decomposition, windows, exclusion of 2020-2021, lags, leave-one-province-out, likelihoods, graph and prior, RSR, Spatial+, survey-based covariate and survey-operation checks, post hoc models |
| 6 | `99FL4_flu_maps.py` | Maps |
| 7 | `04_분석/FLU/MS/build/` | Manuscript numbers, tables, figures, and supplementary material |

The boundary file is found at `data/boundary/final.shp` by default; set `ST_SHP` to use another path. R jobs are wrapped in `rslot.sh` and run with one INLA thread.

## Software

R 4.6.0 with R-INLA, spdep, and sf; Python 3.9 with pandas, geopandas, numpy, matplotlib, and python-docx.

## License

MIT (see `LICENSE`).
