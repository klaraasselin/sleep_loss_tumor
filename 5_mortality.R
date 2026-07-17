#' ---
#' title: "5. Mortality"
#' author: ""
#' date: ""
#' output:
#'   pdf_document:
#'     toc: true
#'     toc_depth: 6
#'     number_sections: true
#' header-includes:
#'   - |
#'     ```{=latex}
#'     \usepackage{fvextra}
#'     \DefineVerbatimEnvironment{Highlighting}{Verbatim}{
#'       breaksymbolleft={}, 
#'       showspaces=false, 
#'       showtabs=false, 
#'       breaklines, 
#'       commandchars=\\\{\}
#'     }
#'     ```
#' ---

#+ include=FALSE, warning=FALSE, message=FALSE, results='hide'

knitr::opts_knit$set(root.dir = "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/these/travail/3_annexes/ix_sommeil/260715/")

knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, 
                      fig.width = 17, fig.height = 10)

options(tibble.print_max = Inf, tibble.width = Inf)

#+ include=TRUE, warning=FALSE, message=FALSE, results='hide'

#' # Environment preparation

#' Run scripts with the needed packages and the raw data frame.

source("0a_packages.R")
source("0b_script_data.R")

#+ include = TRUE

#' # Counts and cumulative proportions

#' ## By lineage, tumoral status, and sleep condition

df %>%
  group_by(lineage, tum_state, sleep_condition) %>%
  summarise(
    n_total = n(),
    n_dead = sum(death_state, na.rm = TRUE),
    prop = round(n_dead / n_total * 100, 2),
    .groups = "drop")

#' ## By lineage and tumoral status

df %>%
  group_by(lineage, tum_state) %>%
  summarise(
    n_total = n(),
    n_dead = sum(death_state, na.rm = TRUE),
    prop = round(n_dead / n_total * 100, 2),
    .groups = "drop")

#' ## By lineage and sleep condition

df %>%
  group_by(lineage, sleep_condition) %>%
  summarise(
    n_total = n(),
    n_dead = sum(death_state, na.rm = TRUE),
    prop = round(n_dead / n_total * 100, 2),
    .groups = "drop")

#' ## By tumoral status and sleep condition

df %>%
  group_by(tum_state, sleep_condition) %>%
  summarise(
    n_total = n(),
    n_dead = sum(death_state, na.rm = TRUE),
    prop = round(n_dead / n_total * 100, 2),
    .groups = "drop")

#' ## By lineage

df %>%
  group_by(lineage) %>%
  summarise(
    n_total = n(),
    n_dead = sum(death_state, na.rm = TRUE),
    prop = round(n_dead / n_total * 100, 2),
    .groups = "drop")

#' ## By tumoral status

df %>%
  group_by(tum_state) %>%
  summarise(
    n_total = n(),
    n_dead = sum(death_state, na.rm = TRUE),
    prop = round(n_dead / n_total * 100, 2),
    .groups = "drop")

#' ## By sleep condition

df %>%
  group_by(sleep_condition) %>%
  summarise(
    n_total = n(),
    n_dead = sum(death_state, na.rm = TRUE),
    prop = round(n_dead / n_total * 100, 2),
    .groups = "drop")