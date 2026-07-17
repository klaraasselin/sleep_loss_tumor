#' # Functions to generate Cox-Snell residual plots for survival models

#' ## For mixed-effects accelerated failure models

plot_cox_snell_survregme <- 
  function(model, data, time_col = "age_followup", event_col = "event") {
  
  #' Step 1: Get the survival probability S(t) predicted by the model,
  #' evaluated at every unique observed time in the dataset (single call
  #' for efficiency instead of predicting one time point per individual)
    
  time_values <- sort(unique(data[[time_col]]))
  
  surv_matrix <- predict(model, type = "survivor", newdata = data, q = time_values)
  
  #' Step 2: For each individual, extract the predicted survival probability
  #' S(t_i) at their own observed time t_i (rather than at all time points)
  
  time_points <- as.numeric(rownames(surv_matrix))
  
  S_ti <- sapply(seq_len(nrow(data)), function(i) surv_matrix[
    which.min(abs(time_points - data[[time_col]][i])), i])
  
  #' Step 3: Compute Cox-Snell residuals as r_i = -log(S(t_i))
  
  cs_resid <- -log(S_ti)
  
  #' Step 4: Treat the Cox-Snell residuals themselves as "survival times"
  #' (using the original event indicator) and estimate their cumulative
  #' hazard via the Nelson-Aalen estimator. 
  
  cs_surv <- survfit(Surv(cs_resid, data[[event_col]]) ~ 1, 
                     type = "fleming-harrington")
  
  cs_df <- data.frame(time = cs_surv$time, 
                      H_hat = -log(cs_surv$surv), 
                      H_lower = -log(cs_surv$upper), 
                      H_upper = -log(cs_surv$lower))
  
  #' Add a starting point at (0,0) so the first step of the curve starts from the origin
  
  cs_df <- rbind(data.frame(time = 0, H_hat = 0, H_lower = 0, H_upper = 0), 
                 cs_df)
  
  #' Converts the confidence ribbon into a "staircase" shape
  #' by duplicating each point to create horizontal segments followed
  #' by vertical segments (matches the step function of the estimate)
  
  stepify_ribbon <- function(df, time_col = "time") {
    df <- df[order(df[[time_col]]), ]
    n  <- nrow(df)
    if (n < 2) return(df)
    out <- df[rep(seq_len(n), each = 2), ]
    out[[time_col]][seq(2, 2 * n, by = 2)] <- c(df[[time_col]][-1], 
                                                df[[time_col]][n])
    return(out) 
    }
  
  cs_ribbon <- stepify_ribbon(cs_df, time_col = "time")
  
  #' Plot: cumulative hazard of the Cox-Snell residuals against the
  #' 45-degree reference line expected under a perfect fit
  
  ggplot(cs_df, aes(x = time, y = H_hat)) +
    geom_ribbon(data = cs_ribbon, 
                aes(x = time, ymin = H_lower, ymax = H_upper, fill = "95% CI"), 
                alpha = 0.2) +
    geom_step(aes(color = "Nelson-Aalen estimate"), 
              direction = "hv", 
              linewidth = 0.9) +
    geom_abline(aes(linetype = "Unit exponential (perfect fit)"), 
                intercept = 0, 
                slope = 1, 
                color = "red", 
                linewidth = 0.8) +
    scale_color_manual(name = NULL, 
                       values = c("Nelson-Aalen estimate" = "steelblue")) +
    scale_fill_manual(name = NULL, 
                      values = c("95% CI" = "steelblue")) +
    scale_linetype_manual(name = NULL, 
                          values = c("Unit exponential (perfect fit)" = "dashed")) +
    labs(x = "Cox-Snell Residuals", 
         y = "Cumulative Hazard H(r)", 
         title = "Cox-Snell Residual Plot") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
          legend.position = "top", 
          legend.justification = "left") 
  
  }

#' ## For mixed-effects cox models

plot_cox_snell_cox <- 
  function(model, data, time_col = "age_followup", event_col = "event") {
  
  #' Step 1: Extract the linear predictor from the coxme model
  #' (combines fixed and random effects)
    
  lp <- predict(model, type = "lp")
  
  #' Step 2: Fit an auxiliary coxph model with the linear predictor as
  #' an offset (fixed, not estimated). This is a workaround to obtain
  #' martingale residuals, which coxme does not provide directly
  
  pseudo_fit <- coxph(Surv(data[[time_col]], data[[event_col]]) ~ offset(lp))
  
  #' Step 3: Derive Cox-Snell residuals from the martingale residuals
  #' (r_i = event_i - M_i). 
  cs_resid <- data[[event_col]] - residuals(pseudo_fit, type = "martingale")
  
  #' Step 4: Estimate the cumulative hazard of the Cox-Snell residuals
  #' via the Nelson-Aalen estimator.
  
  cs_surv <- survfit(Surv(cs_resid, data[[event_col]]) ~ 1, 
                     type = "fleming-harrington")
  
  cs_df <- data.frame(time = cs_surv$time, 
                      H_hat = -log(cs_surv$surv), 
                      H_lower = -log(cs_surv$upper), 
                      H_upper = -log(cs_surv$lower))
  
  #' Add a starting point at (0,0) so the first step of the curve starts from the origin
  
  cs_df <- rbind(data.frame(time = 0, H_hat = 0, H_lower = 0, H_upper = 0), 
                 cs_df)
  
  #' Converts the confidence ribbon into a "staircase" shape
  
  stepify_ribbon <- function(df, time_col = "time") {
    df <- df[order(df[[time_col]]), ]
    n  <- nrow(df)
    if (n < 2) return(df)
    out <- df[rep(seq_len(n), each = 2), ]
    out[[time_col]][seq(2, 2 * n, by = 2)] <- c(df[[time_col]][-1], 
                                                df[[time_col]][n])
    return(out) 
    }
  
  cs_ribbon <- stepify_ribbon(cs_df, time_col = "time")
  
  #' Plot: cumulative hazard of the Cox-Snell residuals against the
  #' 45-degree reference line expected under a perfect fit
  
  ggplot(cs_df, aes(x = time, y = H_hat)) +
    geom_ribbon(data = cs_ribbon,
                aes(x = time, ymin = H_lower, ymax = H_upper, fill = "95% CI"),
                alpha = 0.2) +
    geom_step(aes(color = "Nelson-Aalen estimate"), 
              direction = "hv", 
              linewidth = 0.9) +
    geom_abline(aes(linetype = "Unit exponential (perfect fit)"), 
                intercept = 0, 
                slope = 1, 
                color = "red", 
                linewidth = 0.8) +
    scale_color_manual(name = NULL, 
                       values = c("Nelson-Aalen estimate" = "steelblue")) +
    scale_fill_manual(name = NULL, 
                      values = c("95% CI" = "steelblue")) +
    scale_linetype_manual(name = NULL, 
                          values = c("Unit exponential (perfect fit)" = "dashed")) +
    labs(x = "Cox-Snell Residuals", 
         y = "Cumulative Hazard H(r)", 
         title = "Cox-Snell Residual Plot") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top", 
          legend.justification = "left")
  
  }