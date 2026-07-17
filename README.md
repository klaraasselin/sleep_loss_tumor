This repository contains the raw dataset and the codes used to perform the analyses in the following publication:
Meliani J, Asselin K, Berriahi Z, Renard J, Tökölyi J, M. Nedelcu A, Hamede R, Ujvari B, Thomas F, M. Dujon A. 2026 Sleep loss modulates tumor susceptibility and growth in Hydra. AAAAAA (doi:)

To reproduce the analyses, scripts must be run in the following order:
`packages.R` $\rightarrow$ `script_data.R` $\rightarrow$ `functions.R` $\rightarrow$ any analysis script

Note on format: The analysis scripts are written as annotated R scripts designed to be compiled using `knitr::spin()`.

### Data File
*   **`data.xlsx`**: Raw data frame used in all analyses.

    *   `replicate`: Unique identifier of the experimental replicate.
    *   `id`: Unique identifier of the individual Hydra.
    *   `parent`: Identifier of the parent Hydra (for maternal/genetic lineage tracking).
    *   `date_birth`: Birth date of the Hydra.
    *   `treatment`: Experimental group ("Control" or "Sleep deprivation").
    *   `lineage`: Genetic lineage of the Hydra.
    *   `date_tumor`: First date when a tumor was observed (or `NA` if healthy).
    *   `date_death`: Death date of the Hydra (or `NA` if censored).
    *   `date_FB`: Date when the first bud was detached/produced.
    *   `buds_1` to `buds_6`: Weekly total number of buds produced asexually (weeks 1 to 6).
    *   `tenta_1` to `tenta_6`: Weekly tentacle count (weeks 1 to 6).

### Setup scripts

*   **`packages.R`**: Loads and installs all R libraries required for the project.
*   **`script_data.R`**: Imports the raw Excel dataset, handles data cleaning, formatting, and variable transformations.
*   **`functions.R`**: Custom utility functions (specifically `plot_cox_snell_survregme` to evaluate goodness-of-fit of parametric survival models).

### Analysis scripts

*   **`1_tumoral_dynamics.R`**: Survival analysis (AFT models via `SurvregME`) focusing on tumor onset and development rates.
*   **`2_supernumerary_tentacles.R`**: Survival regressions investigating the occurrence and timing of supernumerary tentacles.
*   **`3_first_bud.R`**: Survival models evaluating the delay before the first asexual budding event.
*   **`4_bud_rate.R`**: Generalized linear mixed-effects models (GLMMs) analyzing weekly budding rates.
*   **`5_mortality.R`**: Computes mortality counts and cumulative proportions of survival across treatments and lineages.
