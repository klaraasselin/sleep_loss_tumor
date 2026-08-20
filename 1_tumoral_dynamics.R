#' ---
#' title: "1. Tumoral dynamics"
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

tum_df <- 
  
  #' Load the raw data frame.
  df %>%
  
  #' Restrict the data frame to the three lineages that developed tumors.
  filter(lineage %in% c("HO_MT", "HO_SPC", "HO_SPT")) %>%
  
  #' Define the follow-up time for each individual:
  #' if the hydra developed a tumor -> time of tumor onset (event),
  #' if the hydra died -> time of death (censor),
  #' if the hydra did not develop a tumor and stayed alive
  #' -> 42 days (no event, censor at the end of the experiment).
  mutate(age_followup_days = case_when(
    tum_state == "Tumoral" ~ age_tum,    
    tum_state == "Healthy" & death_state == 1 ~ age_death,
    TRUE ~ 42),
    
    age_followup = age_followup_days / 7,
    tum_event = ifelse(tum_state == "Tumoral", 1, 0)) %>%
  
  #' Build a competing-risk variable for the cumulative incidence functions:
  #' if the hydra did not develop a tumor and stayed alive -> 0 (no event),
  #' if the hydra developed a tumor -> 1 (event),
  #' if the hydra died -> 2 (competing risk of death).
  mutate(status_competing = factor(case_when(tum_event == 1 ~ 1, 
                                             death_state == 1 ~ 2, 
                                             TRUE ~ 0))) %>%
  
  #' Remove missing data.
  filter(!is.na(age_followup)) %>%
  
  #' Set factor levels for later use in the plots.
  mutate(sleep_condition = factor(sleep_condition, 
                                  levels = c("Control", "Sleep deprivation"))) %>%
  
  #' Keep only the variables needed for the analysis.
  select(age_followup, tum_event, status_competing, sleep_condition, 
         lineage, parent_unique, replicate_unique)

#' Set reference levels for the categorical variables.
tum_df$lineage <- relevel(factor(tum_df$lineage), ref = "HO_MT")
tum_df$sleep_condition <- relevel(factor(tum_df$sleep_condition), ref = "Control")

#' # Survival object

#' Create a survival object with follow-up time and event indicator.
tum_surv_object <- Surv(tum_df$age_followup, tum_df$tum_event)

#' # Proportional hazards assumption check : Schoenfeld residuals

#' Fit a preliminary full Cox model (fixed effects only).
tum_cox <- coxph(Surv(age_followup, tum_event) ~ sleep_condition * lineage, data = tum_df)

#' Test the proportional hazards assumption.
print(cox.zph(tum_cox))

#' The proportional hazards assumption is violated for at least one term, 
#' an accelerated failure time mixed-effects model will be used. 
#' 
#' # AIC comparison of different distributions
#' 
#' Compare intercept-only accelerated failure time models fitted with different 
#' parametric distributions to identify the best-fitting one.

dist_weib <- SurvregME(Surv(age_followup, tum_event) ~ 1, data = tum_df, 
                       dist = "weibull")

dist_lnorm <- SurvregME(Surv(age_followup, tum_event) ~ 1, data = tum_df, 
                        dist = "lognormal")

dist_llog <- SurvregME(Surv(age_followup, tum_event) ~ 1, data = tum_df, 
                       dist = "loglogistic")

dist_exp <- SurvregME(Surv(age_followup, tum_event) ~ 1, data = tum_df, 
                      dist = "exponential")

AIC(dist_weib, dist_lnorm, dist_llog, dist_exp) %>% arrange(AIC)

#' The model with a Weibull distribution has the lowest AIC, we'll use it for the models. 
#'
#' # Model selection
#' 
#' ## Random effects selection
#' 
#' With the full fixed-effect structure held constant, compare different
#' random-effect structures to find the optimal one.

tum_R_0 <- SurvregME(Surv(age_followup, tum_event) ~ sleep_condition * lineage, 
                     data = tum_df, dist = "weibull")

tum_R_1 <- SurvregME(Surv(age_followup, tum_event) ~ sleep_condition * lineage + 
                       (1 | parent_unique), 
                     data = tum_df, dist = "weibull")

tum_R_2 <- SurvregME(Surv(age_followup, tum_event) ~ sleep_condition * lineage +
                       (1 | replicate_unique), 
                     data = tum_df, dist = "weibull")

tum_R_3 <- SurvregME(Surv(age_followup, tum_event) ~ sleep_condition * lineage + 
                       (1 | parent_unique) + (1 | replicate_unique), 
                     data = tum_df, dist = "weibull") 

AIC(tum_R_0, tum_R_1, tum_R_2, tum_R_3) %>% arrange(AIC)

#' The model without any random effect (tum_R_0) has the lowest AIC, this random
#' effect structure is therefore selected.
#' SurvregME is nonetheless kept as the modeling function (rather than switching to a simpler survreg call) 
#' since it also supports models with no random effects.
#' 
#' ## Fixed effects selection
#' 
#' With the selected random-effect structure fixed, compare different 
#' fixed-effect structures to find the optimal one.

tum_F_0 <- SurvregME(Surv(age_followup, tum_event) ~ 1, 
                     data = tum_df, dist = "weibull")

tum_F_1 <- SurvregME(Surv(age_followup, tum_event) ~ sleep_condition, 
                     data = tum_df, dist = "weibull")

tum_F_2 <- SurvregME(Surv(age_followup, tum_event) ~ lineage, 
                     data = tum_df, dist = "weibull")

tum_F_3 <- SurvregME(Surv(age_followup, tum_event) ~ sleep_condition + lineage, 
                     data = tum_df, dist = "weibull")

tum_F_4 <- SurvregME(Surv(age_followup, tum_event) ~ sleep_condition * lineage, 
                     data = tum_df, dist = "weibull")

AIC(tum_F_0, tum_F_1, tum_F_2, tum_F_3, tum_F_4) %>% arrange(AIC)

#' The model with the interaction between the sleep condition and the lineage has 
#' the lowest AIC and is therefore selected (tum_F_4).
#' 
#' ## Goodness-of-fit evaluation
#' 
#' Assess goodness-of-fit of the selected model using Cox-Snell residuals.
plot_cox_snell_survregme(model = tum_F_4, data = tum_df, 
                         time_col = "age_followup", event_col = "tum_event")

#' The curve follows a 45 degree exponential, we validate the fit of the selected model. 
#' 
#' ## Selected model summary
#' 
#' Releveling of "lineage" to obtain all pairwise comparisons from the same selected model.

tum_df$lineage <- relevel(tum_df$lineage, ref = "HO_MT")
summary(SurvregME(Surv(age_followup, tum_event) ~ sleep_condition * lineage, 
                  data = tum_df, dist = "weibull"))

tum_df$lineage <- relevel(tum_df$lineage, ref = "HO_SPC")
summary(SurvregME(Surv(age_followup, tum_event) ~ sleep_condition * lineage, 
                  data = tum_df, dist = "weibull"))

tum_df$lineage <- relevel(tum_df$lineage, ref = "HO_SPT")
summary(SurvregME(Surv(age_followup, tum_event) ~ sleep_condition * lineage, 
                  data = tum_df, dist = "weibull"))

#' # Median times

#' ## By sleep condition

#' Estimate cumulative incidence functions for the tumor event,
#' accounting for the competing risk of death, stratified by sleep condition.
tum_cuminc_sleep <- cuminc(Surv(age_followup, status_competing) ~ 
                             sleep_condition, data = tum_df)

#' For each stratum, extract the median time to event and its 95CI.
tum_medians_sleep <- tidy(tum_cuminc_sleep) %>%
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

print(tum_medians_sleep)

#' ## By lineage

#' Same method as above, but stratified by lineage.

tum_cuminc_lineage <- cuminc(Surv(age_followup, status_competing) ~ lineage, data = tum_df)

tum_medians_lineage <- tidy(tum_cuminc_lineage) %>%
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

print(tum_medians_lineage)

#' ## By sleep condition and lineage

#' Same method as above, but stratified by sleep condition and lineage.

tum_cuminc_interaction <- cuminc(Surv(age_followup, status_competing) ~ 
                                   sleep_condition + lineage, data = tum_df)

tum_medians_interaction <- tidy(tum_cuminc_interaction) %>%
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

print(tum_medians_interaction)

#' # Cumulative incidences

#' ## By week and sleep condition

#' Convert the continuous cumulative incidence function into a weekly summary table:
#' for each stratum of sleep condition and each week, show the estimate with its 95CI.
tum_incidence_sleep <- tidy(tum_cuminc_sleep) %>%
  dplyr::filter(outcome == "1") %>%
  dplyr::mutate(week = ceiling(time)) %>%
  dplyr::group_by(strata, week) %>%
  dplyr::filter(time == max(time)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(strata, week, estimate_cuminc = estimate, conf.low, conf.high) %>%
  dplyr::arrange(strata, week)

print(tum_incidence_sleep)

#' ## By week and lineage

#' Same method as above, but stratified by lineage.

tum_incidence_lineage <- tidy(tum_cuminc_lineage) %>%
  dplyr::filter(outcome == "1") %>%
  dplyr::mutate(week = ceiling(time)) %>%
  dplyr::group_by(strata, week) %>%
  dplyr::filter(time == max(time)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(strata, week, estimate_cuminc = estimate, conf.low, conf.high) %>%
  dplyr::arrange(strata, week)

print(tum_incidence_lineage)

#' ## By week, sleep condition, and lineage

#' Same method as above, but stratified by sleep condition and lineage.

tum_incidence_interaction <- tidy(tum_cuminc_interaction) %>%
  dplyr::filter(outcome == "1") %>%
  dplyr::mutate(week = ceiling(time)) %>%
  dplyr::group_by(strata, week) %>%
  dplyr::filter(time == max(time)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(strata, week, estimate_cuminc = estimate, conf.low, conf.high) %>%
  dplyr::arrange(strata, week)

print(tum_incidence_interaction)

#' # Cumulative proportions

#' ## By week and sleep condition

#' Compute raw cumulative proportions of tumor occurrence per week, 
#' stratified by sleep condition.
tum_proportions_sleep <- tum_df %>%
  dplyr::mutate(week_tum_event = ceiling(age_followup)) %>%
  dplyr::group_by(sleep_condition) %>%
  dplyr::summarise(total_N = n(),
                   n1 = sum(week_tum_event <= 1 & tum_event == 1),
                   n2 = sum(week_tum_event <= 2 & tum_event == 1),
                   n3 = sum(week_tum_event <= 3 & tum_event == 1),
                   n4 = sum(week_tum_event <= 4 & tum_event == 1),
                   n5 = sum(week_tum_event <= 5 & tum_event == 1),
                   n6 = sum(week_tum_event <= 6 & tum_event == 1),
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

print(tum_proportions_sleep)

#' ## By week and lineage

#' Same method as above, but stratified by lineage.

tum_proportions_lineage <- tum_df %>%
  dplyr::mutate(week_tum_event = ceiling(age_followup)) %>%
  dplyr::group_by(lineage) %>%
  dplyr::summarise(total_N = n(),
                   n1 = sum(week_tum_event <= 1 & tum_event == 1),
                   n2 = sum(week_tum_event <= 2 & tum_event == 1),
                   n3 = sum(week_tum_event <= 3 & tum_event == 1),
                   n4 = sum(week_tum_event <= 4 & tum_event == 1),
                   n5 = sum(week_tum_event <= 5 & tum_event == 1),
                   n6 = sum(week_tum_event <= 6 & tum_event == 1),
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

print(tum_proportions_lineage)

#' ## By week, sleep condition, and lineage

#' Same method as above, but stratified by sleep condition and lineage.

tum_proportions_interaction <- tum_df %>%
  dplyr::mutate(week_tum_event = ceiling(age_followup)) %>%
  dplyr::group_by(sleep_condition, lineage) %>%
  dplyr::summarise(total_N = n(),
                   n1 = sum(week_tum_event <= 1 & tum_event == 1),
                   n2 = sum(week_tum_event <= 2 & tum_event == 1),
                   n3 = sum(week_tum_event <= 3 & tum_event == 1),
                   n4 = sum(week_tum_event <= 4 & tum_event == 1),
                   n5 = sum(week_tum_event <= 5 & tum_event == 1),
                   n6 = sum(week_tum_event <= 6 & tum_event == 1),
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

print(tum_proportions_interaction)

#' # Plots of cumulative incidence of tumoral development

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

#' Theme of x-axis
scale_x <- scale_x_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2))

#' Color scales for sleep condition
sleep_scale <- list(
  scale_color_manual(values = sleep_condition_colors, breaks = c("Control", "Sleep deprivation")),
  scale_fill_manual(values = sleep_condition_colors, breaks = c("Control", "Sleep deprivation")))

#' Theme shared by all panels (angle = 90 for bottom-to-top reading)
theme_common <- theme(
  axis.text.x = element_text(color = "black", size = 14),
  axis.title.x = element_text(color = "black", size = 14, margin = margin(t = 10)),
  axis.title.y.right = element_text(color = "black", size = 14, angle = 90, margin = margin(l = 15)),
  plot.margin = margin(l = 5, r = 5, t = 5, b = 5))

#' Specific label configurations for precise placement
labs_base   <- labs(x = NULL, y = NULL, title = NULL)
labs_x_only <- labs(x = "Time (weeks)", y = NULL, title = NULL)
labs_y_only <- labs(x = NULL, y = "Cumulative incidence\nof tumor development", title = NULL)
labs_xy     <- labs(x = "Time (weeks)", y = "Cumulative incidence\nof tumor development", title = NULL)

#' ## Panel A: by sleep condition

#' Cumulative incidence of tumor onset by sleep condition.
A <- cuminc(Surv(age_followup, factor(status_competing)) ~ sleep_condition, 
            data = tum_df) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_xy + scale_y + sleep_scale + scale_x +
  theme_bw() + theme_common + show_y +
  theme(legend.position = "none") +
  plot_annotation(title = "Sleep condition effect",
                  theme = theme(plot.title = element_text(
                    hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10))))

#' Extract legend from plot A to use later on the full figure.
shared_legend <- cowplot::get_legend(
  A + theme(legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 14, face = "bold")))

#' ## Panel B: by sleep condition for each lineage

#' ### HO_MT

#' Cumulative incidence of tumor onset by sleep condition for the lineage HO_MT.
B1 <- cuminc(Surv(age_followup, factor(status_competing)) ~ sleep_condition, 
             data = tum_df %>% filter(lineage == "HO_MT")) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_base + scale_y + scale_x + sleep_scale +
  theme_bw() + theme_common + hide_y +
  theme(legend.position = "none") +
  annotate("text", x = 0.2, y = 0.95, label = "HO_MT",
           color = "black", fontface = "bold", size = 5, hjust = 0, vjust = 1)

#' ### HO_SPC

#' Same method as above, but for the lineage HO_SPC.
B2 <- cuminc(Surv(age_followup, factor(status_competing)) ~ sleep_condition,
             data = tum_df %>% filter(lineage == "HO_SPC")) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_x_only + scale_y + scale_x + sleep_scale + 
  theme_bw() + theme_common + hide_y +
  theme(legend.position = "none") +
  annotate("text", x = 0.2, y = 0.95, label = "HO_SPC",
           color = "black", fontface = "bold", size = 5, hjust = 0, vjust = 1)

#' ### HO_SPT

#' Same method as above, but for the lineage HO_SPT.
B3 <- cuminc(Surv(age_followup, factor(status_competing)) ~ sleep_condition,
             data = tum_df %>% filter(lineage == "HO_SPT")) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_y_only + scale_y + scale_x + sleep_scale +
  theme_bw() + theme_common + show_y +
  theme(legend.position = "none") +
  annotate("text", x = 0.2, y = 0.95, label = "HO_SPT",
           color = "black", fontface = "bold", size = 5, hjust = 0, vjust = 1)

#' Create panel B combining plots B1, B2, and B3
B <- (B1 + B2 + B3 + plot_layout(nrow = 1)) +
  plot_annotation(title = "Interaction between sleep condition and lineage",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", 
                                                          size = 16, margin = margin(b = 10))))

#' ## Panel C: by lineage

#' ### HO_MT

#' Cumulative incidence of tumor onset for the lineage HO_MT.
C1 <- cuminc(Surv(age_followup, factor(status_competing)) ~ 1, 
             data = tum_df %>% filter(lineage == "HO_MT")) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_base + scale_y + scale_x +
  theme_bw() + theme_common + hide_y +
  annotate("text", x = 0.2, y = 0.95, label = "HO_MT",
           color = "black", fontface = "bold", size = 5, hjust = 0, vjust = 1)

#' ### HO_SPC

#' Same method as above, but for the lineage HO_SPC.
C2 <- cuminc(Surv(age_followup, factor(status_competing)) ~ 1, 
             data = tum_df %>% filter(lineage == "HO_SPC")) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_x_only + scale_y + scale_x +
  theme_bw() + theme_common + hide_y +
  annotate("text", x = 0.2, y = 0.95, label = "HO_SPC",
           color = "black", fontface = "bold", size = 5, hjust = 0, vjust = 1)

#' ### HO_SPT

#' Same method as above, but for the lineage HO_SPT.
C3 <- cuminc(Surv(age_followup, factor(status_competing)) ~ 1, 
             data = tum_df %>% filter(lineage == "HO_SPT")) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  labs_y_only + scale_y + scale_x +
  theme_bw() + theme_common + show_y +
  annotate("text", x = 0.2, y = 0.95, label = "HO_SPT",
           color = "black", fontface = "bold", size = 5, hjust = 0, vjust = 1)

#' Create panel C combining plots C1, C2, and C3
C <- (C1 + C2 + C3 + plot_layout(nrow = 1)) +
  plot_annotation(title = "Lineage effect",
                  theme = theme(plot.title = element_text(
                    hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10))))

#' ## Assembly

#' Add a letter tag (A/B/C) to the first sub-plot of each panel.
A <- A + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold"))
B[[1]] <- B[[1]] + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold"))
C[[1]] <- C[[1]] + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold"))

#' Convert each panel into a graphical object.
#' (Wrapping A in patchwork to force perfect alignment with the B and C grids).
A <- patchwork::patchworkGrob(A + patchwork::plot_layout(ncol = 1))
B <- patchwork::patchworkGrob(B)
C <- patchwork::patchworkGrob(C)

#' Combine the panels and the shared legend vertically.
#' The legend is placed between B and C.
final_tum_panel <- cowplot::plot_grid(A, 
                                   B, 
                                   shared_legend, 
                                   C, 
                                   ncol = 1, 
                                   align = "v", 
                                   axis = "lr", 
                                   rel_heights = c(1, 1, 0.1, 1)) +
  theme(plot.margin = margin(l = 5, r = 5, t = 5, b = 5))

#' Show the final figure.
print(final_tum_panel)

ggsave("figure_2.png", plot = final_tum_panel, width = 21, height = 29.7, units = "cm", dpi = 300)
