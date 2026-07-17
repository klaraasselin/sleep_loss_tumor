#' ---
#' title: "4. Bud rate"
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

#' Run scripts with the needed packages, functions, and the raw data frame.

source("0a_packages.R")
source("0b_script_data.R")
source("0c_functions.R")

#+ include=TRUE, warning=FALSE, message=FALSE, results='show'

#' # Data formatting

#' Create a new data frame adapted for generalized linear mixed-effects models (GLMM).

buds_df <- 
  
  #' Load the raw data frame.
  df %>%
  
  #' Reshape the weekly bud counts from wide to long format.
  pivot_longer(cols = starts_with("buds_"),
               names_to = c(".value", "week"), 
               names_sep = "_",) %>%
  
  #' Convert week and bud count to integers.
  mutate(week = as.integer(week),
         buds = as.integer(buds)) %>% 
  
  #' Remove weeks with missing bud count.
  filter(!is.na(buds)) %>%
  
  #' Keep only the variables needed for the analysis.
  select(sleep_condition, lineage, tum_state, week, buds, id_unique, parent_unique, replicate_unique) 

#' Set reference levels for the categorical variables.
buds_df$lineage <- relevel(factor(buds_df$lineage), ref = "HO_MT")
buds_df$sleep_condition <- relevel(factor(buds_df$sleep_condition), ref = "Control")
buds_df$tum_state <- relevel(factor(buds_df$tum_state), ref = "Healthy")

#' # Model selection
#' 
#' ## Random effects selection
#' 
#' With the full fixed-effect structure held constant, compare different
#' random-effect structures to find the optimal one.

buds_pois_R_00 <- glmmTMB(buds ~ week * lineage * sleep_condition, REML = TRUE, family = poisson, data = buds_df)
buds_pois_R_01 <- glmmTMB(buds ~ week * lineage * sleep_condition + (1 | id_unique), REML = TRUE, family = poisson, data = buds_df)
buds_pois_R_02 <- glmmTMB(buds ~ week * lineage * sleep_condition + (week | id_unique), REML = TRUE, family = poisson, data = buds_df)
buds_pois_R_03 <- glmmTMB(buds ~ week * lineage * sleep_condition + (1 | parent_unique), REML = TRUE, family = poisson, data = buds_df)
buds_pois_R_04 <- glmmTMB(buds ~ week * lineage * sleep_condition + (1 | replicate_unique), REML = TRUE, family = poisson, data = buds_df)
buds_pois_R_05 <- glmmTMB(buds ~ week * lineage * sleep_condition + (1 | id_unique) + (1 | parent_unique), REML = TRUE, family = poisson, data = buds_df)
buds_pois_R_06 <- glmmTMB(buds ~ week * lineage * sleep_condition + (1 | id_unique) + (1 | replicate_unique), REML = TRUE, family = poisson, data = buds_df)
buds_pois_R_07 <- glmmTMB(buds ~ week * lineage * sleep_condition + (1 | id_unique) + (1 | parent_unique) + (1 | replicate_unique), REML = TRUE, family = poisson, data = buds_df)

buds_pois_R_names <- paste0("buds_pois_R_", sprintf("%02d", c(0:7)))
buds_pois_R_objects <- mget(buds_pois_R_names)
print(head(do.call(model.sel, buds_pois_R_objects), 5))

#' The model with the replicate id as a random effect (buds_pois_R_04) has the lowest AIC, 
#' this random effect structure is therefore selected.
#' 
#' ## Fixed effects selection
#' 
#' With the selected random-effect structure fixed, compare different 
#' fixed-effect structures to find the optimal one.

buds_pois_F_00 <- glmmTMB(buds ~ 1 + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_01 <- glmmTMB(buds ~ week + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_02 <- glmmTMB(buds ~ lineage + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_03 <- glmmTMB(buds ~ sleep_condition + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_04 <- glmmTMB(buds ~ tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_05 <- glmmTMB(buds ~ week + lineage + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_06 <- glmmTMB(buds ~ week + sleep_condition + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_07 <- glmmTMB(buds ~ week + tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_08 <- glmmTMB(buds ~ lineage + sleep_condition + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_09 <- glmmTMB(buds ~ lineage + tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_10 <- glmmTMB(buds ~ sleep_condition + tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_11 <- glmmTMB(buds ~ week * lineage + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_12 <- glmmTMB(buds ~ week * sleep_condition + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_13 <- glmmTMB(buds ~ week * tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_14 <- glmmTMB(buds ~ lineage * sleep_condition + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_15 <- glmmTMB(buds ~ lineage * tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_16 <- glmmTMB(buds ~ sleep_condition * tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_17 <- glmmTMB(buds ~ week + lineage + sleep_condition + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_18 <- glmmTMB(buds ~ week + lineage + tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_19 <- glmmTMB(buds ~ week + sleep_condition + tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_20 <- glmmTMB(buds ~ lineage + sleep_condition + tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_21 <- glmmTMB(buds ~ week * lineage + sleep_condition + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_22 <- glmmTMB(buds ~ week * sleep_condition + lineage + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_23 <- glmmTMB(buds ~ week * tum_state + lineage + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_24 <- glmmTMB(buds ~ week * lineage + tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_25 <- glmmTMB(buds ~ week * sleep_condition + tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_26 <- glmmTMB(buds ~ week * tum_state + sleep_condition + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_27 <- glmmTMB(buds ~ lineage * sleep_condition + week + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_28 <- glmmTMB(buds ~ lineage * tum_state + week + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_29 <- glmmTMB(buds ~ sleep_condition * tum_state + week + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_30 <- glmmTMB(buds ~ lineage * sleep_condition + tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_31 <- glmmTMB(buds ~ lineage * tum_state + sleep_condition + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_32 <- glmmTMB(buds ~ sleep_condition * tum_state + lineage + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_33 <- glmmTMB(buds ~ week * lineage * sleep_condition + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_34 <- glmmTMB(buds ~ week * lineage * tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_35 <- glmmTMB(buds ~ week * sleep_condition * tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)
buds_pois_F_36 <- glmmTMB(buds ~ lineage * sleep_condition * tum_state + (1 | replicate_unique), REML = FALSE, family = poisson, data = buds_df)

buds_pois_F_names <- paste0("buds_pois_F_", sprintf("%02d", c(0:36)))
buds_pois_F_objects <- mget(buds_pois_F_names)
print(head(do.call(model.sel, buds_pois_F_objects), 5))

#' The model with the interaction between the week, the lineage and the sleep condition
#' has the lowest AIC and is therefore selected (buds_pois_F_33).
#' 
#' ## Goodness-of-fit evaluation

#' Refit the selected model.
buds_pois_mod_33 <- glmmTMB(buds ~ week * lineage * sleep_condition + (1 | replicate_unique), REML = TRUE, family = poisson, data = buds_df)

#' Assess goodness-of-fit of the selected model using simulated residuals.
buds_pois_mod_33_res <- simulateResiduals(buds_pois_mod_33, plot = TRUE)

#' ## Selected model summary

# Releveling of "lineage" to obtain all pairwise comparisons from the same model.

buds_df$lineage <- relevel(buds_df$lineage, ref = "HO_MT")
summary(glmmTMB(buds ~ week * lineage * sleep_condition + (1 | replicate_unique), REML = TRUE, family = poisson, data = buds_df))

buds_df$lineage <- relevel(buds_df$lineage, ref = "HO_SPC")
summary(glmmTMB(buds ~ week * lineage * sleep_condition + (1 | replicate_unique), REML = TRUE, family = poisson, data = buds_df))

buds_df$lineage <- relevel(buds_df$lineage, ref = "HO_SPT")
summary(glmmTMB(buds ~ week * lineage * sleep_condition + (1 | replicate_unique), REML = TRUE, family = poisson, data = buds_df))

buds_df$lineage <- relevel(buds_df$lineage, ref = "HO_VLN")
summary(glmmTMB(buds ~ week * lineage * sleep_condition + (1 | replicate_unique), REML = TRUE, family = poisson, data = buds_df))

buds_df$lineage <- relevel(buds_df$lineage, ref = "HC_MT")
summary(glmmTMB(buds ~ week * lineage * sleep_condition + (1 | replicate_unique), REML = TRUE, family = poisson, data = buds_df))

buds_df$lineage <- relevel(buds_df$lineage, ref = "HV_GAL")
summary(glmmTMB(buds ~ week * lineage * sleep_condition + (1 | replicate_unique), REML = TRUE, family = poisson, data = buds_df))

#' # Estimated marginal means

#' ## Table

#' Compute the estimated marginal means of weekly bud production for every
#' combination of week, sleep condition, and lineage, based on the selected model.

buds_emmeans <- 
  as.data.frame(ggemmeans(buds_pois_mod_33,
                          terms = c("week [all]", "sleep_condition", "lineage"),
                          bias_correction = TRUE)) %>%
  arrange(facet, group, x)

#' ## Plots

#' ## Common functions

#' Function to plot, for a single lineage, the predicted weekly bud
#' production by sleep condition.
buds_plot_treatment <- function(lineage, show_y_text = FALSE, show_x_text = TRUE, show_y_title = FALSE, show_x_title = FALSE) {
  data_sub <- subset(buds_emmeans, facet == lineage)
  
  ggplot(data_sub, aes(x = x, y = predicted, color = group, fill = group)) +
    geom_line(linewidth = 0.8) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
    labs(x = if(show_x_title) "Age (weeks)" else NULL, 
         y = if(show_y_title) "Estimated marginal means\nof the number of buds produced" else NULL, 
         title = NULL) +
    theme_bw() +
    theme(axis.text.x = if(show_x_text) element_text(color = "black", size = 14) else element_blank(),
          axis.ticks.x = if(show_x_text) element_line() else element_blank(),
          axis.text.y = if(show_y_text) element_text(color = "black", size = 14) else element_blank(),
          axis.ticks.y = if(show_y_text) element_line() else element_blank(),
          axis.title.x = element_text(color = "black", size = 14, margin = margin(t = 10)),
          axis.title.y.right = element_text(color = "black", size = 14, angle = 90, margin = margin(l = 25)),
          axis.line.y.right = element_blank(),
          plot.margin = margin(l = 5, r = 5, t = 5, b = 5),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 6), breaks = c(2, 4, 6)) +
    scale_y_continuous(position = "right") +
    coord_cartesian(ylim = c(0, 3.5)) +
    scale_color_manual(values = sleep_condition_colors) +
    scale_fill_manual(values = sleep_condition_colors) +
    annotate("text", x = 1.2, y = 3.3, label = lineage, 
             color = "black", fontface = "bold", size = 5, hjust = 0, vjust = 1)}

#' Function to plot, for a single lineage, the predicted weekly bud
#' production pooled across sleep conditions
buds_plot_lineage <- function(lineage, show_y_text = FALSE, show_x_text = TRUE, show_y_title = FALSE, show_x_title = FALSE) {
  data_sub <- subset(buds_emmeans, facet == lineage)
  data_lineage <- aggregate(cbind(predicted, conf.low, conf.high) ~ x, data = data_sub, FUN = mean)
  data_lineage <- data_lineage[order(data_lineage$x), ]
  
  ggplot(data_lineage, aes(x = x, y = predicted)) +
    geom_line(linewidth = 0.8, color = "black") +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, fill = "black", linewidth = 0.3) +
    labs(x = if(show_x_title) "Age (weeks)" else NULL, 
         y = if(show_y_title) "Estimated marginal means\nof the number of buds produced" else NULL, 
         title = NULL) +
    theme_bw() +
    theme(axis.text.x = if(show_x_text) element_text(color = "black", size = 14) else element_blank(),
          axis.ticks.x = if(show_x_text) element_line() else element_blank(),
          axis.text.y = if(show_y_text) element_text(color = "black", size = 14) else element_blank(),
          axis.ticks.y = if(show_y_text) element_line() else element_blank(),
          axis.title.x = element_text(color = "black", size = 14, margin = margin(t = 10)),
          axis.title.y.right = element_text(color = "black", size = 14, angle = 90, margin = margin(l = 25)),
          axis.line.y.right = element_blank(),
          plot.margin = margin(l = 5, r = 5, t = 5, b = 5)) +
    scale_x_continuous(limits = c(1, 6), breaks = c(2, 4, 6)) +
    scale_y_continuous(position = "right") +
    coord_cartesian(ylim = c(0, 3.5)) +
    annotate("text", x = 1.2, y = 3.3, label = lineage, 
             color = "black", fontface = "bold", size = 5, hjust = 0, vjust = 1) }

#' ## Panel A: interaction between sleep condition and lineage

A1 <- buds_plot_treatment("HO_MT", show_y_text = FALSE, show_x_text = FALSE)
A2 <- buds_plot_treatment("HO_SPC", show_y_text = FALSE, show_x_text = FALSE)
A3 <- buds_plot_treatment("HO_SPT", show_y_text = TRUE, show_x_text = FALSE) 
A4 <- buds_plot_treatment("HO_VLN", show_y_text = FALSE, show_x_text = TRUE)
A5 <- buds_plot_treatment("HC_MT", show_y_text = FALSE, show_x_text = TRUE)
A6 <- buds_plot_treatment("HV_GAL", show_y_text = TRUE, show_x_text = TRUE, 
                            show_x_title = TRUE, show_y_title = TRUE) 

#' Combine A1 to A6 into panel "A"
A <- (A1 + A2 + A3) / (A4 + A5 + A6) + plot_layout(heights = c(1, 1))

#' ## Panel B: lineage effect

B1 <- buds_plot_lineage("HO_MT",  show_y_text = FALSE, show_x_text = FALSE)
B2 <- buds_plot_lineage("HO_SPC", show_y_text = FALSE, show_x_text = FALSE)
B3 <- buds_plot_lineage("HO_SPT", show_y_text = TRUE,  show_x_text = FALSE) 
B4 <- buds_plot_lineage("HO_VLN", show_y_text = FALSE, show_x_text = TRUE)
B5 <- buds_plot_lineage("HC_MT",  show_y_text = FALSE, show_x_text = TRUE)
B6 <- buds_plot_lineage("HV_GAL", show_y_text = TRUE,  show_x_text = TRUE, show_x_title = TRUE, show_y_title = TRUE)

#' Combine B1 to B6 into panel "B"
B <- (B1 + B2 + B3) / (B4 + B5 + B6) + plot_layout(heights = c(1, 1))

#' Extract legend from the plot A6 to use later on the full figure. 
legend_treatment <- cowplot::get_legend(A6 +
                                          labs(color = NULL, fill = NULL) +
                                          theme(legend.position = "right", legend.text = element_text(size = 13)))

#' ## Assembly

#' Add a letter tag (A/B) to the first sub-plot of each panel.
A[[1]][[1]] <- A[[1]][[1]] + labs(tag = "A") + 
  theme(plot.tag = element_text(size = 20, face = "bold"))
B[[1]][[1]] <- B[[1]][[1]] + labs(tag = "B") + 
  theme(plot.tag = element_text(size = 20, face = "bold"))

#' Convert each panel into a graphical object.
A <- patchwork::patchworkGrob(A)
B <- patchwork::patchworkGrob(B)

#' Build standalone title graphical objects for each panel.
interaction_title_grob <- cowplot::ggdraw() + 
  cowplot::draw_label("        Interaction between sleep condition and lineage", 
                      fontface = "bold", size = 16, x = 0.5, hjust = 0.5)

lineage_title_grob <- cowplot::ggdraw() + 
  cowplot::draw_label("     Lineage effect", 
                      fontface = "bold", size = 16, x = 0.5, hjust = 0.5)

#' Stack the two titles and their corresponding panels vertically.
left_column <- cowplot::plot_grid(interaction_title_grob,
                                  A, 
                                  lineage_title_grob, 
                                  B, 
                                  ncol = 1, 
                                  align = "v", 
                                  axis = "r", 
                                  rel_heights = c(0.05, 1, 0.05, 1))

#' Shift the sleep-condition legend to the right, to align it with the plots.
legend_treatment_shifted <- cowplot::plot_grid(NULL, legend_treatment, ncol = 2, rel_widths = c(0.6, 0.1))

#' Place the legend in a right-hand column.
right_column <- cowplot::plot_grid(legend_treatment_shifted, NULL, ncol = 1, rel_heights = c(0.9, 4))

#' Combine the panels and the legend. 
final_buds_panel <- cowplot::plot_grid(left_column, right_column, ncol = 2, rel_widths = c(1, 0.16)) +
  theme(plot.margin = margin(l = 5, r = 80, t = 5, b = 5))

#' Show the final figure.
print(final_buds_panel)
