# Early‑Life‑Antibodies

[![phiper on GitHub](https://img.shields.io/badge/phiper-GitHub-black?logo=github)](https://github.com/Polymerase3/phiper)
[![styled with styler](https://img.shields.io/badge/styled%20with-styler-black?logo=R)](https://github.com/r-lib/styler)

**R analyses of Lifelines mother-baby PhIP‑Seq cohorts**

Analysis pipelines for Lifelines infants and mothers (serum & breast‑milk). Alpha & beta diversity, t‑SNEs, visualisations, GLMMs and reproducible workflows driven by the development version of the `phiper` package .

`phiper` (development): https://github.com/Polymerase3/phiper

---

## Status / short warning

This repository and `phiper` are under active development. Some scripts assume the current (development) `phiper` on my GitHub and may break if versions differ. The code is intended to be run interactively in R/RStudio - it *can* be adapted for HPC but that is not the primary target here.

**Important:** the analyses can be memory-intensive. I recommend at least **64 GB RAM** for full runs. Some scripts work on 32 GB but heatmap/distance steps may fail.

---

## Quick start

1. Clone the repository (SSH):
```bash
mkdir -p ~/projects
cd ~/projects
git clone git@github.com:Polymerase3/Early-Life-Antibodies.git
# or https://github.com/Polymerase3/Early-Life-Antibodies.git
cd Early-Life-Antibodies

# optional sanity check
git remote -v
```
If SSH fails: add your SSH public key to GitHub.

2. Open the project in RStudio:
- `File --> New Project --> Existing Directory` --> select `~/projects/Early-Life-Antibodies` --> **Create Project**.  

Confirm:
```r
getwd()   # should end with .../Early-Life-Antibodies
```

3. Run the mandatory first script **interactively** (recommended):
- Open `scripts/01_data-prep.R` and run chunk-by-chunk (Ctrl+Enter / Run) so you can watch installs and long steps.
- To run at once:
```r
# in R (project root)
source("scripts/01_data-prep.R")
```

`01_data-prep.R` will install/load packages (including dev `phiper` if needed), prepare the DuckDB backend, import `.parquet` reads & metadata and create the `phip_data` object used by downstream scripts. After it finishes you can run `02-tsnes.R`, `03-heatmaps.R`, etc. interactively.

---

## Installation / dependencies

* R (latest recommended). Packages used are listed at the top of `scripts/01_data-prep.R`; please run that script first - it installs or loads the required packages.
* Development `phiper` package (my package for efficient PhIP‑Seq handling). The first script attempts to install it automatically from my GitHub. If you want to install manually:

```r
# example (remotes must be installed)
remotes::install_github("Polymerase3/phiper")
```

* DuckDB is used as the on‑disk backend (handled inside `phiper`).

**Hardware**: minimum recommended **64 GB RAM**. Heavy distance/heatmap steps require more memory and temporary space.

---

## Repository structure

```
├── data/                 # input parquet files, precomputed distance matrices, metadata
├── results/              # generated outputs (created when you run scripts)
├── scripts/              # analysis scripts (run interactively)
├── README.md             # this file
```

### scripts (short descriptions)

* `01_data-prep.R` — **MANDATORY**. Builds DuckDB, imports parquet PhIP‑Seq reads, appends metadata, recodes time variables and creates a single `phip_data` long format object (~90M rows in our case). Run this before anything else.

* `02-tsnes.R` — t‑SNE explorations across metadata. Quick and dirty; written for rapid testing (not fully optimised). All t‑SNE outputs are saved to `results/`.

* `03-heatmaps.R` — heatmaps and stability boxplots (Figures 2/3 + supplements). Computes stability using Kulczynski distance (chosen for stability across large richness differences).

* `04-outliers.R` — detect richness outliers with density/histograms/QQ plots and Mahalanobis distances.

* `05-stability_meta.R` — correlate stability metrics with metadata for infants; outputs boxplots and FDR‑corrected tables.

All scripts expect `data/` and `results/` to be in the repository root. `results/` is created automatically when scripts run.

---

## Data

* Raw PhIP‑Seq reads: stored as `.parquet` files under `data/`. `phiper` handles efficient import into DuckDB.
* Precomputed distance matrices, metadata tables, and sample collection time files are present in `data/` so scripts can re‑use them without recomputing everything.

If you add or replace files, keep filenames and relative paths consistent with the scripts (they look for files in `data/`).

---

## Outliers (detected)

The following sample IDs were flagged as outliers during the richness / QC pipeline:

```
LLNEXT303926    R27P03_08_303926-M-W2-Melk-2_LLNext_P10_A_T_C2
LLNEXT003520    R27P02_23_003520-M-M3-Melk-8_LLNext_P9_A_T_C2
LLNEXT304815    R27P03_07_304815-M-W2-Melk-2_LLNext_P10_A_T_C2
LLNEXT008048    R26P03_53_008048-M-M12-Melk-1_LLNext_P8_A_T_C2
LLNEXT120669    R27P04_28_120669-M-M3-Melk-1_LLNext_P11_A_T_C2
LLNEXT202113    R27P03_64_202113-M-M3-Melk-2_LLNext_P10_A_T_C2
LLNEXT200736    R27P03_20_200736-M-M1-Melk-4_LLNext_P10_A_T_C2
LLNEXT011783    R27P03_49_011783-M-M1-Melk-6_LLNext_P10_A_T_C2
```

If you re‑run the pipeline with different parameters these may change — treat this list as the current run's flags, not a definitive blacklist.

---

## Reproducibility notes

* The workflow uses a DuckDB on‑disk backend to avoid loading massive tables into RAM.
* Many operations use `dplyr`/`dbplyr` via `phiper`; remember to `collect()` only *after* heavy filtering/summarisation and **avoid** collecting full large objects into memory.
* Version your `phiper` and R package environment (e.g., `renv` or `pak`) if you want strict reproducibility across machines.

---

## Contributing

This repo is experimental and maintained by me. If you want to contribute:

* Email: **[mateusz.kolek@meduniwien.ac.at](mailto:mateusz.kolek@meduniwien.ac.at)** — brief description of the change you propose and which script you want to touch. Or alternatively open an issue directly on GitHub.
* Small, focused pull requests are preferred. Add tests or a short example when changing core data‑handling code.

---

## Contact

Mateusz Kolek — [mateusz.kolek@meduniwien.ac.at](mailto:mateusz.kolek@meduniwien.ac.at)

If something is unclear or a script fails, send the error message and the output of `sessionInfo()`.
