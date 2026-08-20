#' ---
#' title: "3. First bud production"
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

knitr::opts_knit$set(root.dir = "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/these/travail/3_annexes/ix_sommeil/260820/")

knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, 
                      fig.width = 17, fig.height = 10)

options(tibble.print_max = Inf, tibble.width = Inf)

#+ include=TRUE, warning=FALSE, message=FALSE, results='hide'

#' # Environment preparation

#' Run scripts with the needed packages, functions, and the raw data frame.

source("0a_packages.R")
source("0b_script_data.R")
source("0c_functions.R")

#+ include=TRUE, warning=FALSE, message=FALSE, results='show'

#' # Data formatting

#' Create a new data frame adapted for survival analysis.

FB_df <- 
  
  #' Load the raw data frame.
  df %>%
  
  #' Define the follow-up time for each individual:
  #' if the hydra produced a first bud -> time of first bud production (event),
  #' if the hydra died -> time of death (censor),
  #' if the hydra did not produce a first bud and stayed alive
  #' -> 42 days (no event, censor at the end of the experiment).
  mutate(age_followup_days = case_when(
    FB_state == 1 ~ age_FB,    
    FB_state == 0 & death_state == 1 ~ age_death,
    TRUE ~ 42),
    
    age_followup = age_followup_days / 7,
    FB_event = FB_state) %>%
  
  #' Build a competing-risk variable for the cumulative incidence functions:
  #' if the hydra did not produce a first bud and stayed alive -> 0 (no event),
  #' if the hydra produced a first bud -> 1 (event),
  #' if the hydra died -> 2 (competing risk of death). 
  mutate(status_competing = factor(case_when(FB_event == 1 ~ 1, 
                                             death_state == 1 ~ 2, 
                                             TRUE ~ 0))) %>%
  
  #' Build the variable tum_state_early: an individual is considered
  #' tumoral only if the tumor appeared before or at the same time as the
  #' supernumerary tentacle event.
  mutate(age_tum_weeks = age_tum / 7,
         tum_state_early = case_when(tum_state == "Tumoral" & 
                                       !is.na(age_tum_weeks) & 
                                       age_tum_weeks <= age_followup ~ "Tumoral",
                                     TRUE ~ "Healthy"),
         tum_state_early = factor(tum_state_early, levels = c("Healthy", "Tumoral"))) %>%
  
  #' Remove missing data.
  filter(!is.na(age_followup)) %>%
  
  #' Keep only the variables needed for the analysis.
  select(age_followup, FB_event, status_competing, sleep_condition, lineage, tum_state_early, parent_unique, replicate_unique)

#' Set reference levels for the categorical variables.
FB_df$lineage <- relevel(factor(FB_df$lineage), ref = "HO_MT")
FB_df$sleep_condition <- relevel(factor(FB_df$sleep_condition), ref = "Control")

#' # Survival object

#' Create a survival object with follow-up time and event indicator.

FB_surv_object <- Surv(FB_df$age_followup, FB_df$FB_event)

#' # Proportional hazards assumption check : Schoenfeld residuals

#' Fit a preliminary full Cox model (fixed effects only).
FB_cox <- coxph(Surv(age_followup, FB_event) ~ sleep_condition * lineage * tum_state_early, data = FB_df)

#' Test the proportional hazards assumption.
print(cox.zph(FB_cox))

#' The proportional hazards assumption is violated for at least one term, 
#' an accelerated failure time mixed-effects model will be used. 
#' 
#' # AIC comparison of different distributions
#' 
#' Compare intercept-only accelerated failure time models fitted with different 
#' parametric distributions to identify the best-fitting one.

FB_dist_weib <- SurvregME(Surv(age_followup, FB_event) ~ 1, data = FB_df, 
                          dist = "weibull")

FB_dist_lnorm <- SurvregME(Surv(age_followup, FB_event) ~ 1, data = FB_df, 
                           dist = "lognormal")

FB_dist_llog <- SurvregME(Surv(age_followup, FB_event) ~ 1, data = FB_df, 
                          dist = "loglogistic")

FB_dist_exp <- SurvregME(Surv(age_followup, FB_event) ~ 1, data = FB_df, 
                         dist = "exponential")

AIC(FB_dist_weib, FB_dist_lnorm, FB_dist_llog, FB_dist_exp) %>% arrange(AIC)

#' The model with a log-normal distribution has the lowest AIC, we'll use it for the models. 
#'
#' # Model selection
#' 
#' ## Random effects selection
#' 
#' With the full fixed-effect structure held constant, compare different
#' random-effect structures to find the optimal one.

FB_R_0 <- SurvregME(Surv(age_followup, FB_event) ~ 
                      sleep_condition * lineage * tum_state_early, 
                    data = FB_df, dist = "lognormal")

FB_R_1 <- SurvregME(Surv(age_followup, FB_event) ~ 
                      sleep_condition * lineage * tum_state_early + 
                      (1 | parent_unique), 
                    data = FB_df, dist = "lognormal")

FB_R_2 <- SurvregME(Surv(age_followup, FB_event) ~ 
                      sleep_condition * lineage * tum_state_early + 
                      (1 | replicate_unique), 
                    data = FB_df, dist = "lognormal")

FB_R_3 <- SurvregME(Surv(age_followup, FB_event) ~ 
                      sleep_condition * lineage * tum_state_early + 
                      (1 | parent_unique) + (1 | replicate_unique), 
                    data = FB_df, dist = "lognormal") 

AIC(FB_R_0, FB_R_1, FB_R_2, FB_R_3) %>% arrange(AIC)

#' The model with the replicate id as a random effect (FB_R_2) has the lowest AIC, 
#' this random effect structure is therefore selected.
#' 
#' ## Fixed effects selection
#' 
#' With the selected random-effect structure fixed, compare different 
#' fixed-effect structures to find the optimal one.

FB_F_0 <- SurvregME(Surv(age_followup, FB_event) ~ 1 + (1 | replicate_unique), data = FB_df, dist = "lognormal")
FB_F_1 <- SurvregME(Surv(age_followup, FB_event) ~ sleep_condition + (1 | replicate_unique), data = FB_df, dist = "lognormal")
FB_F_2 <- SurvregME(Surv(age_followup, FB_event) ~ lineage + (1 | replicate_unique), data = FB_df, dist = "lognormal")
FB_F_3 <- SurvregME(Surv(age_followup, FB_event) ~ tum_state_early + (1 | replicate_unique), data = FB_df, dist = "lognormal")
FB_F_4 <- SurvregME(Surv(age_followup, FB_event) ~ sleep_condition + lineage + (1 | replicate_unique), data = FB_df, dist = "lognormal")
FB_F_5 <- SurvregME(Surv(age_followup, FB_event) ~ sleep_condition + tum_state_early + (1 | replicate_unique), data = FB_df, dist = "lognormal")
FB_F_6 <- SurvregME(Surv(age_followup, FB_event) ~ lineage + tum_state_early + (1 | replicate_unique), data = FB_df, dist = "lognormal")
FB_F_7 <- SurvregME(Surv(age_followup, FB_event) ~ sleep_condition * lineage + (1 | replicate_unique), data = FB_df, dist = "lognormal")
FB_F_8 <- SurvregME(Surv(age_followup, FB_event) ~ sleep_condition * tum_state_early + (1 | replicate_unique), data = FB_df, dist = "lognormal")
FB_F_9 <- SurvregME(Surv(age_followup, FB_event) ~ lineage * tum_state_early + (1 | replicate_unique), data = FB_df, dist = "lognormal")
FB_F_10 <- SurvregME(Surv(age_followup, FB_event) ~ sleep_condition + lineage + tum_state_early + (1 | replicate_unique), data = FB_df, dist = "lognormal")
FB_F_11 <- SurvregME(Surv(age_followup, FB_event) ~ sleep_condition * lineage + tum_state_early + (1 | replicate_unique), data = FB_df, dist = "lognormal")
FB_F_12 <- SurvregME(Surv(age_followup, FB_event) ~ sleep_condition * tum_state_early + lineage + (1 | replicate_unique), data = FB_df, dist = "lognormal")
FB_F_13 <- SurvregME(Surv(age_followup, FB_event) ~ lineage * tum_state_early + sleep_condition + (1 | replicate_unique), data = FB_df, dist = "lognormal")
FB_F_14 <- SurvregME(Surv(age_followup, FB_event) ~ sleep_condition * lineage * tum_state_early + (1 | replicate_unique), data = FB_df, dist = "lognormal")

AIC(FB_F_0, FB_F_1, FB_F_2, FB_F_3, FB_F_4, FB_F_5, FB_F_6, FB_F_7, 
    FB_F_8, FB_F_9, FB_F_10, FB_F_11, FB_F_12, FB_F_13, FB_F_14) %>% arrange(AIC)

#' The model with the interaction between the sleep condition and the lineage, and with
#' the additive effect of tumoral state has the lowest AIC and is therefore selected (FB_F_11).
#' 
#' ## Goodness-of-fit evaluation
#' 
#' Assess goodness-of-fit of the selected model using Cox-Snell residuals.

plot_cox_snell_survregme(model = FB_F_11, data = FB_df, time_col  = "age_followup", event_col = "FB_event")

#' The curve follows a 45 degree exponential, we validate the fit of the selected model. 
#' 
#' ## Selected model summary
#' 
#' #' Releveling of "lineage" to obtain all pairwise comparisons from the same model.

FB_df$lineage <- relevel(FB_df$lineage, ref = "HO_MT")
summary(SurvregME(Surv(age_followup, FB_event) ~ 
                    sleep_condition * lineage + tum_state_early + (1 | replicate_unique), 
                  data = FB_df, dist = "lognormal"))

FB_df$lineage <- relevel(FB_df$lineage, ref = "HO_SPC")
summary(SurvregME(Surv(age_followup, FB_event) ~ 
                    sleep_condition * lineage + tum_state_early + (1 | replicate_unique), 
                  data = FB_df, dist = "lognormal"))

FB_df$lineage <- relevel(FB_df$lineage, ref = "HO_SPT")
summary(SurvregME(Surv(age_followup, FB_event) ~ 
                    sleep_condition * lineage + tum_state_early + (1 | replicate_unique), 
                  data = FB_df, dist = "lognormal"))

FB_df$lineage <- relevel(FB_df$lineage, ref = "HO_VLN")
summary(SurvregME(Surv(age_followup, FB_event) ~ 
                    sleep_condition * lineage + tum_state_early + (1 | replicate_unique), 
                  data = FB_df, dist = "lognormal"))

FB_df$lineage <- relevel(FB_df$lineage, ref = "HC_MT")
summary(SurvregME(Surv(age_followup, FB_event) ~ 
                    sleep_condition * lineage + tum_state_early + (1 | replicate_unique), 
                  data = FB_df, dist = "lognormal"))

FB_df$lineage <- relevel(FB_df$lineage, ref = "HV_GAL")
summary(SurvregME(Surv(age_followup, FB_event) ~ 
                    sleep_condition * lineage + tum_state_early + (1 | replicate_unique), 
                  data = FB_df, dist = "lognormal"))

#' # Median times

#' ## By sleep condition

#' Estimate cumulative incidence functions for the first bud event,
#' accounting for the competing risk of death, stratified by sleep condition.

FB_cuminc_sleep <- cuminc(Surv(age_followup, status_competing) ~ sleep_condition, data = FB_df)

#' For each stratum, extract the median time to event and its 95CI.

FB_medians_sleep <- tidy(FB_cuminc_sleep) %>%
  dplyr::filter(outcome == "1") %>%
  dplyr::group_by(strata) %>%
  dplyr::arrange(time) %>%
  dplyr::summarise(
    median_time = { idx <- which(estimate >= 0.5)[1]
    if (is.na(idx)) NA_real_ else time[idx] },
    ci_low = { idx <- which(conf.high >= 0.5)[1]
    if (is.na(idx)) NA_real_ else time[idx] },
    ci_high = { idx <- which(conf.low >= 0.5)[1]
    if (is.na(idx)) NA_real_ else time[idx] },
    .groups = "drop")

print(FB_medians_sleep)

#' ## By lineage

#' Same method as above, but stratified by lineage.

FB_cuminc_lineage <- cuminc(Surv(age_followup, status_competing) ~ lineage, data = FB_df)

FB_medians_lineage <- tidy(FB_cuminc_lineage) %>%
  dplyr::filter(outcome == "1") %>%
  dplyr::group_by(strata) %>%
  dplyr::arrange(time) %>%
  dplyr::summarise(
    median_time = { idx <- which(estimate >= 0.5)[1]
    if (is.na(idx)) NA_real_ else time[idx] },
    ci_low = { idx <- which(conf.high >= 0.5)[1]
    if (is.na(idx)) NA_real_ else time[idx] },
    ci_high = { idx <- which(conf.low >= 0.5)[1]
    if (is.na(idx)) NA_real_ else time[idx] },
    .groups = "drop")

print(FB_medians_lineage)

#' ## By sleep condition and lineage

#' Same method as above, but stratified by sleep condition and lineage.

FB_cuminc_interaction <- cuminc(Surv(age_followup, status_competing) ~ sleep_condition + lineage, data = FB_df)

FB_medians_interaction <- tidy(FB_cuminc_interaction) %>%
  dplyr::filter(outcome == "1") %>%
  dplyr::group_by(strata) %>%
  dplyr::arrange(time) %>%
  dplyr::summarise(
    median_time = { idx <- which(estimate >= 0.5)[1]
    if (is.na(idx)) NA_real_ else time[idx] },
    ci_low = { idx <- which(conf.high >= 0.5)[1]
    if (is.na(idx)) NA_real_ else time[idx] },
    ci_high = { idx <- which(conf.low >= 0.5)[1]
    if (is.na(idx)) NA_real_ else time[idx] },
    .groups = "drop")

print(FB_medians_interaction)

#' ## By tumoral state

#' Same method as above, but stratified by tumoral state.

FB_cuminc_tum <- cuminc(Surv(age_followup, status_competing) ~ tum_state_early, data = FB_df)

FB_medians_tum <- tidy(FB_cuminc_tum) %>%
  dplyr::filter(outcome == "1") %>%
  dplyr::group_by(strata) %>%
  dplyr::arrange(time) %>%
  dplyr::summarise(
    median_time = { idx <- which(estimate >= 0.5)[1]
    if (is.na(idx)) NA_real_ else time[idx] },
    ci_low = { idx <- which(conf.high >= 0.5)[1]
    if (is.na(idx)) NA_real_ else time[idx] },
    ci_high = { idx <- which(conf.low >= 0.5)[1]
    if (is.na(idx)) NA_real_ else time[idx] },
    .groups = "drop")

print(FB_medians_tum)

#' # Cumulative incidences

#' ## By week and sleep condition

#' Convert the continuous cumulative incidence function into a weekly summary table:
#' for each stratum of sleep condition and each week, show the estimate with its 95CI.

FB_incidence_sleep <- tidy(FB_cuminc_sleep) %>%
  dplyr::filter(outcome == "1") %>%
  dplyr::mutate(week = ceiling(time)) %>%
  dplyr::group_by(strata, week) %>%
  dplyr::filter(time == max(time)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(strata, week, estimate_cuminc = estimate, conf.low, conf.high) %>%
  dplyr::arrange(strata, week)

print(FB_incidence_sleep)

#' ## By week and lineage

#' Same method as above, but stratified by lineage.

FB_incidence_lineage <- tidy(FB_cuminc_lineage) %>%
  dplyr::filter(outcome == "1") %>%
  dplyr::mutate(week = ceiling(time)) %>%
  dplyr::group_by(strata, week) %>%
  dplyr::filter(time == max(time)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(strata, week, estimate_cuminc = estimate, conf.low, conf.high) %>%
  dplyr::arrange(strata, week)

print(FB_incidence_lineage)

#' ## By week, sleep condition, and lineage

#' Same method as above, but stratified by sleep condition and lineage.

FB_incidence_interaction <- tidy(FB_cuminc_interaction) %>%
  dplyr::filter(outcome == "1") %>%
  dplyr::mutate(week = ceiling(time)) %>%
  dplyr::group_by(strata, week) %>%
  dplyr::filter(time == max(time)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(strata, week, estimate_cuminc = estimate, conf.low, conf.high) %>%
  dplyr::arrange(strata, week)

print(FB_incidence_interaction)

#' ## By week and tumoral state

#' Same method as above, but stratified by tumoral state. 

FB_incidence_tum <- tidy(FB_cuminc_tum) %>%
  dplyr::filter(outcome == "1") %>%
  dplyr::mutate(week = ceiling(time)) %>%
  dplyr::group_by(strata, week) %>%
  dplyr::filter(time == max(time)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(strata, week, estimate_cuminc = estimate, conf.low, conf.high) %>%
  dplyr::arrange(strata, week)

print(FB_incidence_tum)

#' # Cumulative proportions

#' ## By week and sleep condition

#' Compute raw cumulative proportions of first bud production per week, 
#' stratified by sleep condition.

FB_proportions_sleep <- FB_df %>%
  dplyr::mutate(week_FB_event = ceiling(age_followup)) %>%
  dplyr::group_by(sleep_condition) %>%
  dplyr::summarise(total_N = n(),
                   n1 = sum(week_FB_event <= 1 & FB_event == 1),
                   n2 = sum(week_FB_event <= 2 & FB_event == 1),
                   n3 = sum(week_FB_event <= 3 & FB_event == 1),
                   n4 = sum(week_FB_event <= 4 & FB_event == 1),
                   n5 = sum(week_FB_event <= 5 & FB_event == 1),
                   n6 = sum(week_FB_event <= 6 & FB_event == 1),
                   p1 = n1 / total_N,
                   p2 = n2 / total_N,
                   p3 = n3 / total_N,
                   p4 = n4 / total_N,
                   p5 = n5 / total_N,
                   p6 = n6 / total_N,
                   .groups = "drop") %>%
  tidyr::pivot_longer(cols = matches("^[np][1-6]$"),
                      names_to = c(".value", "week"),
                      names_pattern = "(n|p)([1-6])") %>%
  dplyr::mutate(week = as.numeric(week)) %>%
  dplyr::select(sleep_condition, week, n, total_N, raw_proportion = p)

print(FB_proportions_sleep)

#' ## By week and lineage

#' Same method as above, but stratified by lineage.

FB_proportions_lineage <- FB_df %>%
  dplyr::mutate(week_FB_event = ceiling(age_followup)) %>%
  dplyr::group_by(lineage) %>%
  dplyr::summarise(total_N = n(),
                   n1 = sum(week_FB_event <= 1 & FB_event == 1),
                   n2 = sum(week_FB_event <= 2 & FB_event == 1),
                   n3 = sum(week_FB_event <= 3 & FB_event == 1),
                   n4 = sum(week_FB_event <= 4 & FB_event == 1),
                   n5 = sum(week_FB_event <= 5 & FB_event == 1),
                   n6 = sum(week_FB_event <= 6 & FB_event == 1),
                   p1 = n1 / total_N,
                   p2 = n2 / total_N,
                   p3 = n3 / total_N,
                   p4 = n4 / total_N,
                   p5 = n5 / total_N,
                   p6 = n6 / total_N,
                   .groups = "drop") %>%
  tidyr::pivot_longer(cols = matches("^[np][1-6]$"),
                      names_to = c(".value", "week"),
                      names_pattern = "(n|p)([1-6])") %>%
  dplyr::mutate(week = as.numeric(week)) %>%
  dplyr::select(lineage, week, n, total_N, raw_proportion = p)

print(FB_proportions_lineage)

#' ## By week, sleep condition, and lineage

#' Same method as above, but stratified by sleep condition and lineage.

FB_proportions_interaction <- FB_df %>%
  dplyr::mutate(week_FB_event = ceiling(age_followup)) %>%
  dplyr::group_by(sleep_condition, lineage) %>%
  dplyr::summarise(total_N = n(),
                   n1 = sum(week_FB_event <= 1 & FB_event == 1),
                   n2 = sum(week_FB_event <= 2 & FB_event == 1),
                   n3 = sum(week_FB_event <= 3 & FB_event == 1),
                   n4 = sum(week_FB_event <= 4 & FB_event == 1),
                   n5 = sum(week_FB_event <= 5 & FB_event == 1),
                   n6 = sum(week_FB_event <= 6 & FB_event == 1),
                   p1 = n1 / total_N,
                   p2 = n2 / total_N,
                   p3 = n3 / total_N,
                   p4 = n4 / total_N,
                   p5 = n5 / total_N,
                   p6 = n6 / total_N,
                   .groups = "drop") %>%
  tidyr::pivot_longer(cols = matches("^[np][1-6]$"),
                      names_to = c(".value", "week"),
                      names_pattern = "(n|p)([1-6])") %>%
  dplyr::mutate(week = as.numeric(week)) %>%
  dplyr::select(sleep_condition, lineage, week, n, total_N, raw_proportion = p)

print(FB_proportions_interaction)

#' ## By week and tumoral state

#' Same method as above, but stratified by tumoral state.

FB_proportions_tum <- FB_df %>%
  dplyr::mutate(week_FB_event = ceiling(age_followup)) %>%
  dplyr::group_by(tum_state_early) %>%
  dplyr::summarise(total_N = n(),
                   n1 = sum(week_FB_event <= 1 & FB_event == 1),
                   n2 = sum(week_FB_event <= 2 & FB_event == 1),
                   n3 = sum(week_FB_event <= 3 & FB_event == 1),
                   n4 = sum(week_FB_event <= 4 & FB_event == 1),
                   n5 = sum(week_FB_event <= 5 & FB_event == 1),
                   n6 = sum(week_FB_event <= 6 & FB_event == 1),
                   p1 = n1 / total_N,
                   p2 = n2 / total_N,
                   p3 = n3 / total_N,
                   p4 = n4 / total_N,
                   p5 = n5 / total_N,
                   p6 = n6 / total_N,
                   .groups = "drop") %>%
  tidyr::pivot_longer(cols = matches("^[np][1-6]$"),
                      names_to = c(".value", "week"),
                      names_pattern = "(n|p)([1-6])") %>%
  dplyr::mutate(week = as.numeric(week)) %>%
  dplyr::select(tum_state_early, week, n, total_N, raw_proportion = p)

print(FB_proportions_tum)

#' # Plots of cumulative incidence of first bud production

#' ## Common themes

#' Theme of y-axis, and its presence or absence on the plots
scale_y <- scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), 
                              labels = percent, position = "right")

#' Fix to hide y-axis numbers, ticks, and titles on the right and left
hide_y <- theme(axis.text.y.right = element_blank(), 
                axis.ticks.y.right = element_blank(),
                axis.title.y.right = element_blank(),
                axis.text.y = element_blank(), 
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank())

#' Fix to show y-axis numbers and ticks on the right
show_y <- theme(axis.text.y.right = element_text(color = "black", size = 14), 
                axis.ticks.y.right = element_line(),
                axis.text.y = element_blank(), 
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank())

#' Theme of x-axis, and its presence or absence on the plots.
scale_x <- scale_x_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2))

hide_x <- theme(axis.text.x = element_blank(), 
                axis.ticks.x = element_blank(), 
                axis.title.x = element_blank())

show_x <- theme(axis.text.x = element_text(color = "black", size = 14), 
                axis.ticks.x = element_line())

#' Color scale for tumoral state
tum_scale <- list(
  scale_color_manual(values = tum_colors),
  scale_fill_manual(values = tum_colors))

#' Color scales for sleep condition
sleep_scale <- list(
  scale_color_manual(values = sleep_condition_colors, breaks = c("Control", "Sleep deprivation")),
  scale_fill_manual(values = sleep_condition_colors, breaks = c("Control", "Sleep deprivation")))

#' Theme shared by panels A, B, and C (angle = 90 for bottom-to-top reading)
theme_common <- theme(
  axis.title.x = element_text(color = "black", size = 14, margin = margin(t = 10)),
  axis.title.y.right = element_text(color = "black", size = 14, angle = 90, margin = margin(l = 15)),
  plot.margin = margin(l = 5, r = 5, t = 5, b = 5))

#' Specific label configurations for precise placement
labs_base   <- labs(x = NULL, y = NULL, title = NULL)
labs_x_only <- labs(x = "Time (weeks)", y = NULL, title = NULL)
labs_y_only <- labs(x = NULL, y = "Cumulative incidence of\nfirst bud production", title = NULL)
labs_xy     <- labs(x = "Time (weeks)", y = "Cumulative incidence of\nfirst bud production", title = NULL)

#' ## Panel A: interaction between sleep condition and lineage

#' ### HO_SPT
A1 <- cuminc(Surv(age_followup, factor(status_competing)) ~ sleep_condition, 
             data = FB_df %>% filter(lineage == "HO_SPT") %>%
               mutate(sleep_condition = factor(sleep_condition, levels = c("Control", "Sleep deprivation")))) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_base + scale_y + scale_x + sleep_scale +
  theme_bw() + theme_common + hide_y +
  theme(axis.text.x = element_text(color = "black", size = 14),
        legend.position = "none") +
  annotate("text", x = 5.8, y = 0.05, label = "HO_SPT",
           color = "black", fontface = "bold", size = 5, hjust = 1, vjust = 0)

#' ### HO_VLN
A2 <- cuminc(Surv(age_followup, factor(status_competing)) ~ sleep_condition, 
             data = FB_df %>% filter(lineage == "HO_VLN") %>%
               mutate(sleep_condition = factor(sleep_condition, levels = c("Control", "Sleep deprivation")))) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_x_only + scale_y + scale_x + sleep_scale +
  theme_bw() + theme_common + hide_y +
  theme(axis.text.x = element_text(color = "black", size = 14),
        legend.position = "none") +
  annotate("text", x = 5.8, y = 0.05, label = "HO_VLN",
           color = "black", fontface = "bold", size = 5, hjust = 1, vjust = 0)

#' ### HC_MT
A3 <- cuminc(Surv(age_followup, factor(status_competing)) ~ sleep_condition, 
             data = FB_df %>% filter(lineage == "HC_MT") %>%
               mutate(sleep_condition = factor(sleep_condition, levels = c("Control", "Sleep deprivation")))) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_y_only + scale_y + scale_x + sleep_scale +
  theme_bw() + theme_common + show_y +
  theme(axis.text.x = element_text(color = "black", size = 14),
        legend.position = "none") +
  annotate("text", x = 5.8, y = 0.05, label = "HC_MT",
           color = "black", fontface = "bold", size = 5, hjust = 1, vjust = 0)

#' Extract legend from the plot A3
shared_treatment_legend <- cowplot::get_legend(
  A3 + theme(legend.position = "bottom",
             legend.direction = "horizontal",
             legend.text = element_text(size = 13),
             legend.title = element_text(size = 14, face = "bold")))

#' Combine A1, A2, and A3 into panel "A"
A <- (A1 + A2 + A3 + plot_layout(nrow = 1)) +
  plot_annotation(title = "Interaction between sleep condition and lineage",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10))))

#' ## Panel B: by tumoral state

#' Cumulative incidence of first bud production by tumoral state.
B <- cuminc(Surv(age_followup, factor(status_competing)) ~ tum_state_early, data = FB_df) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_xy + scale_y + tum_scale + scale_x +
  theme_bw() + theme_common + show_y +
  theme(axis.text.x = element_text(color = "black", size = 14),
        axis.text.y.right = element_text(color = "black", size = 14),
        axis.ticks.y.right = element_line(),
        legend.position = "none") +
  plot_annotation(title = "Tumoral state effect",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10))))

#' Extract legend from this panel
shared_tumor_legend <- cowplot::get_legend(
  B + theme(legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 14, face = "bold")))

#' ## Panel C: by lineage

#' ### HO_MT
C1 <- cuminc(Surv(age_followup, factor(status_competing)) ~ 1, 
             data = FB_df %>% filter(lineage == "HO_MT")) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_base + scale_y + scale_x +
  theme_bw() + theme_common + hide_y + hide_x +
  annotate("text", x = 5.8, y = 0.05, label = "HO_MT",
           color = "black", fontface = "bold", size = 5, hjust = 1, vjust = 0)

#' ### HO_SPC
C2 <- cuminc(Surv(age_followup, factor(status_competing)) ~ 1, 
             data = FB_df %>% filter(lineage == "HO_SPC")) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_base + scale_y + scale_x +
  theme_bw() + theme_common + hide_y + hide_x +
  annotate("text", x = 5.8, y = 0.05, label = "HO_SPC",
           color = "black", fontface = "bold", size = 5, hjust = 1, vjust = 0)

#' ### HO_SPT
C3 <- cuminc(Surv(age_followup, factor(status_competing)) ~ 1, 
             data = FB_df %>% filter(lineage == "HO_SPT")) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_base + scale_y + scale_x +
  theme_bw() + theme_common + hide_y + hide_x + 
  theme(axis.text.y.right = element_text(color = "black", size = 14), 
        axis.ticks.y.right = element_line()) +
  annotate("text", x = 5.8, y = 0.05, label = "HO_SPT",
           color = "black", fontface = "bold", size = 5, hjust = 1, vjust = 0)

#' ### HO_VLN
C4 <- cuminc(Surv(age_followup, factor(status_competing)) ~ 1, 
             data = FB_df %>% filter(lineage == "HO_VLN")) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_base + scale_y + scale_x +
  theme_bw() + theme_common + hide_y + show_x +
  annotate("text", x = 5.8, y = 0.05, label = "HO_VLN",
           color = "black", fontface = "bold", size = 5, hjust = 1, vjust = 0)

#' ### HC_MT
C5 <- cuminc(Surv(age_followup, factor(status_competing)) ~ 1, 
             data = FB_df %>% filter(lineage == "HC_MT")) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_x_only + scale_y + scale_x +
  theme_bw() + theme_common + hide_y + show_x +
  annotate("text", x = 5.8, y = 0.05, label = "HC_MT",
           color = "black", fontface = "bold", size = 5, hjust = 1, vjust = 0)

#' ### HV_GAL
C6 <- cuminc(Surv(age_followup, factor(status_competing)) ~ 1, 
             data = FB_df %>% filter(lineage == "HV_GAL")) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_y_only + scale_y + scale_x +
  theme_bw() + theme_common + show_y + show_x +
  annotate("text", x = 5.8, y = 0.05, label = "HV_GAL",
           color = "black", fontface = "bold", size = 5, hjust = 1, vjust = 0)

#' Create panel C combining plots C1 to C6.
C <- (C1 + C2 + C3) / (C4 + C5 + C6) + 
  plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "Lineage effect",
                  theme = theme(plot.title = element_text(
                    hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10))))

#' ## Assembly

#' Add a letter tag to the first sub-plot of each panel.
A[[1]] <- A[[1]] + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold"))
B <- B + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold"))
C[[1]][[1]] <- C[[1]][[1]] + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold"))

#' Convert each panel into a graphical object.
A <- patchwork::patchworkGrob(A)
#' Wrap B in patchwork before conversion to force perfect axis alignment with A and C
B <- patchwork::patchworkGrob(B + patchwork::plot_layout(ncol = 1)) 
C <- patchwork::patchworkGrob(C)

#' Combine panels and shared legends vertically into a single column, aligning axes.
final_FB_panel <- cowplot::plot_grid(A, 
                                     shared_treatment_legend, 
                                     B, 
                                     shared_tumor_legend, 
                                     C, 
                                     ncol = 1, 
                                     align = "v",
                                     axis = "lr",
                                     rel_heights = c(1, 0.1, 1, 0.1, 1.8)) +
  theme(plot.margin = margin(l = 5, r = 5, t = 5, b = 5))

#' Show the final figure.
print(final_FB_panel)

ggsave("figure_3.png", plot = final_FB_panel, width = 21, height = 29.7, units = "cm", dpi = 300)
