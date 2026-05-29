#------------------------------------------------------------------------------#
#------------------Created by Luke Barry & Roch Nianogo------------------------#
#-----------------------------Date:07/23/2025----------------------------------#
#--------------------------Purpose: DCEA Tutorial------------------------------#
#------------------------Look-up List of Parameters----------------------------#
#------------------------------------------------------------------------------#

# Variables required as part of data construction and cleaning
param_vars   <- c(
  "pr_no_cvd_history_mi", "pr_mi_cvd_death", "pr_cvd_history_mi", "pr_noncvd_death",
  "c_cvd_history", "c_no_cvd_history", "c_mi", "u_mi",
  "u_cvd_history", "u_no_cvd_history", "hc_cost", "utility"
)

# Generate ASCVD weights each cycle using values from control group 
# so that weights do not change between treatment scenarios
param_data <- average_ascvd_risk(mydata = nhanes_sample) 

# Calculate transition probabilities for each cycle
# Treated individuals
trt <- param_data %>%
  Probs(mydata = .,
        Trt    = TRUE,
        .delta_sys_bp = delta_sys_bp)

# Untreated individuals
ntrt <- param_data %>%
  Probs(mydata = .,
        Trt    = FALSE,
        .delta_sys_bp = 0)

# Append the two datasets with a treatment indicator
# and calculate costs and utilities for different model states
# Where data is unavailable in NHANES for cost/utility equations, 
# assume the same for everyone (e.g. cci = 0)
combined <- bind_rows(
  trt  %>% mutate(trt = 1),
  ntrt %>% mutate(trt = 0)) %>%
  Costs(mydata = .,
        cci = 0, 
        mi = 0,
        cardiac_dysrhythmia = 0,      
        peripheral_artery_disease = 0) %>% 
  mutate(c_no_cvd_history = hc_cost) %>%
  Costs(mydata = .,
        cci = 0,
        mi = 1,
        cardiac_dysrhythmia = 0,      
        peripheral_artery_disease = 0) %>% 
  mutate(c_cvd_history = hc_cost) %>%
  Costs(mydata = .,
        cci = 0,
        cardiac_dysrhythmia = 0,      
        peripheral_artery_disease = 0) %>% 
  mutate(c_mi = hc_cost*(1-delta_c_mi)) %>%
  Effs(mydata = .,
       cci = 0,
       mi = 0,
       cardiac_dysrhythmia = 0,      
       peripheral_artery_disease = 0) %>% 
  mutate(u_no_cvd_history = utility) %>%
  Effs(mydata = .,
       cci = 0,
       mi = 1,
       cardiac_dysrhythmia = 0,      
       peripheral_artery_disease = 0) %>% 
  mutate(u_cvd_history = utility) %>%
  Effs(mydata = .,
       cci = 0,
       cardiac_dysrhythmia = 0,      
       peripheral_artery_disease = 0) %>% 
  mutate(u_mi = utility*(1-delta_u_mi)) %>%
  group_by(across(all_of(c(group_vars, "trt")))) %>% 
  summarise(across(
    all_of(c(param_vars)), 
    ~mean(.x, na.rm = TRUE),
    .names = "{.col}"),
    .groups = "drop"
    , trt = first(trt)) 


# Expand the data by crossing with the state vector
lookup_list <- combined %>%
  mutate(row_id = row_number()) %>%  # Optional: keep a unique identifier for original rows
  crossing(state = state_names) %>%  # Cartesian product with states
  select(-row_id) %>%                # Remove if you don't need the row ID
  # update probabilities, costs and effects according state and treatment status
  mutate(
    across(
      all_of(param_vars), 
      ~ if_else(state %in% c("cvd_death", "noncvd_death"), 0, .)
    ),
    across(
      all_of(c("pr_no_cvd_history_mi", "pr_cvd_history_mi")), 
      ~ if_else(state %in% c("mi"), 0, .)
    ),
    across(
      all_of(c("pr_no_cvd_history_mi", "pr_mi_cvd_death")), 
      ~ if_else(state %in% c("history_cvd"), 0, .)
    ),
    across(
      all_of(c("pr_cvd_history_mi", "pr_mi_cvd_death")), 
      ~ if_else(state %in% c("no_cvd"), 0, .)
    ),
    cost = case_when(state == "history_cvd" ~ c_cvd_history,
                     state == "no_cvd" ~ c_no_cvd_history,
                     state == "mi" ~ c_mi, 
                     TRUE ~ hc_cost),
    qaly = case_when(state == "history_cvd" ~ u_cvd_history,
                        state == "no_cvd" ~ u_no_cvd_history,
                        state == "mi" ~ u_mi, 
                        TRUE ~ utility)
  ) %>%
  select(- c_cvd_history,
         - c_no_cvd_history,
         - c_mi,
         - u_mi,
         - u_cvd_history,
         - u_no_cvd_history,
         - hc_cost,
         - utility) %>%
  mutate(state     = factor(as.character(state), levels = state_names),
         state_idx = as.integer(state))

# remove unecessary data
rm(trt, ntrt, combined, param_data)

# save the look-up list
saveRDS(lookup_list, here("data", "output_data", "lookup_list.rds"))
