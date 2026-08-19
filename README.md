This repository contains the dataset and the codes used to perform the analyses in the following publication:
Meliani J, Asselin K, Berriahi Z, Renard J, Tökölyi J, M. Nedelcu A, Hamede R, Ujvari B, Thomas F, M. Dujon A. 2026 Sleep loss modulates tumor susceptibility and growth in Hydra. AAAAAA (doi:)

It also contains the output of the codes in PDF format. 

To reproduce the analyses, scripts must be run in the following order:
`packages.R` $\rightarrow$ `script_data.R` $\rightarrow$ `functions.R` $\rightarrow$ any analysis script

The scripts are written as annotated R scripts designed to be compiled using `knitr::spin()`.

*   **`data.xlsx`**: Raw data file used in all analyses.

    *   `replicate`: Id of the replicate.
    *   `id`: Id of the hydra.
    *   `parent`: Id of the parent
    *   `date_birth`: Birth date of the hydra.
    *   `treatment`: Experimental group ("Control" or "Sleep deprivation").
    *   `lineage`: Lineage of the hydra.
    *   `date_tumor`: First date when a tumor was observed (or NA if healthy).
    *   `date_death`: Death date of the hydra (or NA if alive during the experimental period).
    *   `date_FB`: Date when the first bud was produced.
    *   `buds_1` to `buds_6`: Weekly total number of buds produced asexually (weeks 1 to 6).
    *   `tenta_1` to `tenta_6`: Weekly tentacle count (weeks 1 to 6).

### Setup scripts

*   **`0a_packages.R`**: Loads all R libraries required for the analyses.
*   **`0b_script_data.R`**: Imports the raw Excel dataset, handles data cleaning, formatting, and variable transformations.
*   **`0c_functions.R`**: Custom functions to evaluate goodness-of-fit of survival models.

### Analysis scripts

*   **`1_tumoral_dynamics.R`**: Survival models focusing on tumoral dynamics (prevalence and time to tumor onset).
*   **`2_supernumerary_tentacles.R`**: Survival models focusing on the occurence of supernumerary tentacles.
*   **`3_first_bud.R`**: Survival models focusing on the first bud produced asexually.
*   **`4_bud_rate.R`**: Generalized linear mixed-effects models focusing on weekly budding rates.
*   **`5_mortality.R`**: Computes mortality counts and cumulative proportions across treatments and lineages.

### Code outputs

*   **`1_tumoral_dynamics.pdf`**
*   **`2_supernumerary_tentacles.pdf`**
*   **`3_first_bud.pdf`**
*   **`4_bud_rate.pdf`**
*   **`5_mortality.pdf`**
