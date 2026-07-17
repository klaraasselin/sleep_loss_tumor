#' Data import and cleaning

df <- 
  
  #' Import the excel file
  read_excel("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/these/travail/3_annexes/ix_sommeil/260715/data.xlsx", 
             col_types = c("numeric", "numeric", "numeric", "date", "text", 
                           "text", "date", "date", "date", "numeric", "numeric",
                           "numeric", "numeric", "numeric", "numeric", "numeric",
                           "numeric", "numeric", "numeric", "numeric", "numeric")) %>%
  
  mutate(
    #' Replace textual "NA" by a real missing value NA
    across(where(is.character), ~ na_if(.x, "NA")),
    
    #' Convert to numeric type
    across(starts_with("buds_"), as.numeric),
    across(starts_with("tenta_"), as.numeric),
    
    #' Convert to date type
    date_birth = as.Date(date_birth),
    date_FB = as.Date(date_FB),
    date_tum = as.Date(date_tumor),
    date_death = as.Date(date_death),
    
    #' Create binomial variables for events: death, tumor, and first bud
    FB_state = ifelse(!is.na(date_FB), 1, 0),
    tum_state = ifelse(!is.na(date_tum), 1, 0),
    death_state = ifelse(!is.na(date_death), 1, 0),
    
    #' Calculate age in days for events
    age_FB = as.numeric(date_FB - date_birth),
    age_tum = as.numeric(date_tum - date_birth),
    age_death = as.numeric(date_death - date_birth),
    
    #' Homogenize sleep condition names
    sleep_condition = forcats::fct_recode(
      factor(treatment),
      "Control" = "Sleep",
      "Control" = "sleep",
      "Sleep deprivation" = "sleep_deprivation"),
    
    #' Homogenize lineages names 
    lineage = recode(
      lineage,
      "SPbC" = "HO_SPC", 
      "SPbT" = "HO_SPT", 
      "OMT" = "HO_MT", 
      "CMT" = "HC_MT", 
      "VGAL" = "HV_GAL", 
      "VLN" = "HO_VLN"), 
    
    #' Convert tum_state to a factor with labels
    tum_state = factor(tum_state, levels = c(0, 1), 
                            labels = c("Healthy", "Tumoral")),
    
    #' Convert categorical variables to factors
    sleep_condition = factor(as.character(sleep_condition)),
    lineage = factor(as.character(lineage)),
    replicate = factor(as.numeric(replicate)),
    parent = factor(as.numeric(parent)),
    id = factor(as.numeric(id)),
    
    #' Limit data to the experimental period: 6 weeks = 42 days
    age_tum = ifelse(age_tum > 42, NA, age_tum),
    age_death = ifelse(age_death > 42, NA, age_death),
    age_FB = ifelse(age_FB > 42, NA, age_FB),
    
    #' If an age was set to NA above (event outside the experimental period),
    #' the corresponding state is reset to "no event"
    tum_state = case_when(is.na(age_tum) ~ "Healthy",
                          TRUE ~ as.character(tum_state)),
    tum_state = factor(tum_state, levels = c("Healthy", "Tumoral")),
    death_state = ifelse(is.na(age_death), 0, death_state),
    FB_state = ifelse(is.na(age_FB), 0, FB_state),
    
    #' Create unique ids, used later as random effects
    id_unique = paste0(as.character(id), "_", as.character(replicate), "_", sleep_condition, "_", as.character(lineage)),
    parent_unique = paste0(as.character(parent), "_", as.character(lineage)),
    replicate_unique = paste0(as.character(parent), "_", as.character(replicate), "_", as.character(lineage)))

#' Color palettes for the plots

sleep_condition_colors <- c(
  "Control" = "cyan",          
  "Sleep deprivation" = "orange")

tum_colors <- c(
  "Healthy"  = "green", 
  "Tumoral"  = "red")