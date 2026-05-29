#------------------------------------------------------------------------------#
#------------------Created by Luke Barry & Roch Nianogo------------------------#
#-----------------------------Date:07/23/2025----------------------------------#
#--------------------------Purpose: DCEA Tutorial------------------------------#
#---------------------------------Functions------------------------------------#
#------------------------------------------------------------------------------#

#create new folers
create_folders <- function(){
  
  
  folder_names <- c("data", "rmd", "docs", "figures", "tables", "scripts")
  purrr::walk(folder_names, dir.create)
  
  sub_folder_data <- c(here("data", "raw_data"),
                       here("data", "output_data"))
  
  purrr::walk(sub_folder_data, dir.create)
  
  
  sub_folder_docs <- c(here("docs", "Lit"),
                       here("docs", "Manuscript")
  )
  
  purrr::walk(sub_folder_docs, dir.create)
  
}
#e.g. create_folders()

#### Risk Functions ----

  #### ASCVD weighting function ----
  # The ASCVD function predicts an individuals 10 year ASCVD risk and estimates their average risk within age and sex strata
  
  average_ascvd_risk <- function(mydata) {
    
    mydata <- mydata %>%
      
      # Estimate 10 year ASCVD risk
      mutate(ascvd_risk_10yr = case_when(
        # When using NHANES create "black_race" variable and apply this to black individuals as all non-black races are ascribed the white-race risk 
        # see (https://www.framinghamheartstudy.org/fhs-risk-functions/cardiovascular-disease-10-year-risk/) for details
        (race != "black" & female == 0) ~ (1 - 0.9144) * exp(12.344 * log(age) + 
                                                               0      * log(age) * log(age) + 
                                                               11.853 * log(total_chol) + 
                                                               -2.664 * log(age) * log(total_chol) +
                                                               -7.990 * log(hdl_chol) + 
                                                               1.769  * log(age) * log(hdl_chol) + 
                                                               (1.797 * log(sys_bp) + 
                                                                  0      * log(age) * log(sys_bp)) * hyperten_trt + 
                                                               (1.764 * log(sys_bp) + 
                                                                  0      * log(age) * log(sys_bp)) * (1 - hyperten_trt) + 
                                                               7.837  * smoking +
                                                               -1.795 * log(age) * smoking + 
                                                               0.658  * diabetes + 
                                                               -(61.18)),
        (race != "black" & female == 1) ~ (1 - 0.9665) * exp(-29.799   * log(age) + 
                                                               4.884   * log(age)*log(age) + 
                                                               13.54   * log(total_chol) + 
                                                               -3.114  * log(age) * log(total_chol) +
                                                               -13.578 * log(hdl_chol) + 
                                                               3.149   * log(age) * log(hdl_chol) + 
                                                               (2.019  * log(sys_bp) + 
                                                                  0       * log(age) * log(sys_bp)) * hyperten_trt + 
                                                               (1.957  * log(sys_bp) + 
                                                                  0       * log(age) * log(sys_bp)) * (1 - hyperten_trt) + 
                                                               7.574   * smoking +
                                                               -1.665  * log(age) * smoking + 
                                                               0.661   * diabetes + 
                                                               (-(-29.18))),
        (race == "black" & female == 0) ~ (1 - 0.8954) * exp(2.469    * log(age) + 
                                                               0      * log(age) * log(age) + 
                                                               0.302  * log(total_chol) + 
                                                               0      * log(age) * log(total_chol) +
                                                               -0.307 * log(hdl_chol) + 
                                                               0      * log(age) * log(hdl_chol) + 
                                                               (1.916 * log(sys_bp) + 
                                                                  0      * log(age) * log(sys_bp)) * hyperten_trt + 
                                                               (1.809 * log(sys_bp) + 
                                                                  0      * log(age) * log(sys_bp)) * (1 - hyperten_trt) + 
                                                               0.549  * smoking +
                                                               0      * log(age) * smoking + 
                                                               0.645  * diabetes + 
                                                               -(19.54)),
        (race == "black" & female == 1) ~ (1 - 0.9533) * exp(17.114    * log(age) + 
                                                               0       * log(age) * log(age) + 
                                                               0.940   * log(total_chol) + 
                                                               0       * log(age) * log(total_chol) +
                                                               -18.920 * log(hdl_chol) + 
                                                               4.475   * log(age)*log(hdl_chol) + 
                                                               (29.291 * log(sys_bp) +
                                                                  -6.432  * log(age) * log(sys_bp)) * hyperten_trt + 
                                                               (27.820 * log(sys_bp) + 
                                                                  -6.087  * log(age) * log(sys_bp)) * (1 - hyperten_trt) + 
                                                               0.691   * smoking +
                                                               0       * log(age) * smoking + 
                                                               0.874   * diabetes + 
                                                               -(86.61)),
        TRUE~NA_real_),
        
        ascvd_risk_10yr   = ifelse(ascvd_risk_10yr < 1, 
                                   # Equation can result in probabilities > 1 for very high ages; this puts a limit on it [not relevant if restricting model to age < 65]
                                   ascvd_risk_10yr, 0.999999999),               
        # Estimate individual weights to be applied to FHS MI or stroke risk within strata 
        ascvd_rate   = ((-log(1 - ascvd_risk_10yr))/(10)), 
        ascvd_risk   = 1 - exp(-(ascvd_rate) * 1)) %>%
      # Estimate average 10 year ASCVD risk within strata according to age and sex 
      group_by(across(all_of(c("age", "female")))) %>%
      mutate(avg_ascvd_risk = mean(ascvd_risk, na.rm = T)) %>%
      ungroup() %>%
      select(-ascvd_risk
             , -ascvd_rate
             , -ascvd_risk_10yr
             )
    
    return(mydata) 
    
  }
  
  #### Probability function ----
  # The Probs function that updates the transition probabilities of every cycle is shown below.
  
  # Predict annual mi risk for those with no history of CVD from FHS monthly risk re-weighted by ASCVD event risk
  Probs <- function(mydata,
                    Trt = FALSE,
                    .delta_sys_bp = 0) {
    
    # add non-CVD mortality data
    mydata <- left_join(mydata, nonCVD_mort, by = c("age", "female"), relationship = "many-to-many")
    
    # Apply treatment effect to systolic blood pressure if Trt is TRUE
    if (Trt == TRUE) {
      mydata <- mydata %>%
        mutate(sys_bp = sys_bp - (sys_bp * .delta_sys_bp))
    }
    
    mydata  <- mydata %>%
      
      # Estimate 10 year ASCVD risk
      mutate(
        ascvd_risk_10yr = case_when(
        # When using NHANES create "black_race" variable and apply this to black individuals as all non-black races are ascribed the white-race risk 
        # see (https://www.framinghamheartstudy.org/fhs-risk-functions/cardiovascular-disease-10-year-risk/) for details
        (race != "black" & female == 0) ~ (1 - 0.9144) * exp(12.344 * log(age) + 
                                                               0      * log(age) * log(age) + 
                                                               11.853 * log(total_chol) + 
                                                               -2.664 * log(age) * log(total_chol) +
                                                               -7.990 * log(hdl_chol) + 
                                                               1.769  * log(age) * log(hdl_chol) + 
                                                               (1.797 * log(sys_bp) + 
                                                                  0      * log(age) * log(sys_bp)) * hyperten_trt + 
                                                               (1.764 * log(sys_bp) + 
                                                                  0      * log(age) * log(sys_bp)) * (1 - hyperten_trt) + 
                                                               7.837  * smoking +
                                                               -1.795 * log(age) * smoking + 
                                                               0.658  * diabetes + 
                                                               -(61.18)),
        (race != "black" & female == 1) ~ (1 - 0.9665) * exp(-29.799   * log(age) + 
                                                               4.884   * log(age)*log(age) + 
                                                               13.54   * log(total_chol) + 
                                                               -3.114  * log(age) * log(total_chol) +
                                                               -13.578 * log(hdl_chol) + 
                                                               3.149   * log(age) * log(hdl_chol) + 
                                                               (2.019  * log(sys_bp) + 
                                                                  0       * log(age) * log(sys_bp)) * hyperten_trt + 
                                                               (1.957  * log(sys_bp) + 
                                                                  0       * log(age) * log(sys_bp)) * (1 - hyperten_trt) + 
                                                               7.574   * smoking +
                                                               -1.665  * log(age) * smoking + 
                                                               0.661   * diabetes + 
                                                               (-(-29.18))),
        (race == "black" & female == 0) ~ (1 - 0.8954) * exp(2.469    * log(age) + 
                                                               0      * log(age) * log(age) + 
                                                               0.302  * log(total_chol) + 
                                                               0      * log(age) * log(total_chol) +
                                                               -0.307 * log(hdl_chol) + 
                                                               0      * log(age) * log(hdl_chol) + 
                                                               (1.916 * log(sys_bp) + 
                                                                  0      * log(age) * log(sys_bp)) * hyperten_trt + 
                                                               (1.809 * log(sys_bp) + 
                                                                  0      * log(age) * log(sys_bp)) * (1 - hyperten_trt) + 
                                                               0.549  * smoking +
                                                               0      * log(age) * smoking + 
                                                               0.645  * diabetes + 
                                                               -(19.54)),
        (race == "black" & female == 1) ~ (1 - 0.9533) * exp(17.114    * log(age) + 
                                                               0       * log(age) * log(age) + 
                                                               0.940   * log(total_chol) + 
                                                               0       * log(age) * log(total_chol) +
                                                               -18.920 * log(hdl_chol) + 
                                                               4.475   * log(age)*log(hdl_chol) + 
                                                               (29.291 * log(sys_bp) +
                                                                  -6.432  * log(age) * log(sys_bp)) * hyperten_trt + 
                                                               (27.820 * log(sys_bp) + 
                                                                  -6.087  * log(age) * log(sys_bp)) * (1 - hyperten_trt) + 
                                                               0.691   * smoking +
                                                               0       * log(age) * smoking + 
                                                               0.874   * diabetes + 
                                                               -(86.61)),
        TRUE~NA_real_),
        ascvd_risk_10yr   = ifelse(ascvd_risk_10yr < 1, 
                                   # Equation can result in probabilities > 1 for very high ages; this puts a limit on it [not relevant if restricting model to age < 65]
                                   ascvd_risk_10yr, 0.9999999),               
        # Estimate individual weights to be applied to FHS MI or stroke risk within strata 
        ascvd_rate   = ((-log(1 - ascvd_risk_10yr))/(10)), 
        ascvd_risk   = 1 - exp(-(ascvd_rate) * 1),
        
        ## Estimate monthly FHS MI risk, convert to annual rate, apply 10 year ASCVD weights, convert to re-weighted annual MI risk for those without a history of CVD
        # Estimate monthly MI risk as a function of age within categories of gender
        monthly_mi_risk       = ifelse(female == 0,
                                       0.0001   * exp(0.0312 * age),      
                                       0.000008 * exp(0.0599 * age)),
        
        # Equation can result in probabilities > 1 for very high ages; this puts a limit on it
        monthly_mi_risk       = ifelse(monthly_mi_risk < 1, 
                                       monthly_mi_risk, 1),  
        
        # Convert monthly risk to annual rate (assumes constant rate over time)
        annual_mi_rate        = ((-log(1 - monthly_mi_risk))/(1 / 12)),  
        
        # Convert annual rate to annual risk (assumes constant rate over time)  
        annual_mi_risk        = 1 - exp(-(annual_mi_rate) * 1),   
        
        # Using Bayes formula to estimate the probability of an MI given an ASCVD event for age and sex categories [avg_ascvd_risk grouping (age and sex) should be the same as annual_mi_risk parameters (age and sex)]
        annual_ascvd_mi_risk  = annual_mi_risk / avg_ascvd_risk, 
        
        # Use Bayes formula to then estimate the probability of an ASCVD and MI event [this is more granular as the probability of an ASCVD event uses more information, importantly systolic BP]  
        pr_no_cvd_history_mi  = annual_ascvd_mi_risk * ascvd_risk,          

        # Estimate MI mortality risk as a function of age within categories of gender
        pr_mi_cvd_death       = ifelse(female == 0, 
                                       0.0289 * exp(0.0269 * age),        
                                       0.0004 * exp(0.0706 * age)), 
        
        # Equation can result in probabilities > 1 for very high ages; this puts a limit on it
        pr_mi_cvd_death       = ifelse(pr_mi_cvd_death < 1, 
                                       pr_mi_cvd_death, 1),                        

        # Convert annual MI risk to rate, apply "history of CVD multiplier" (see Basu, 2017)  and convert back to prob
        pr_cvd_history_mi     = 1 - exp((-(((-log(1 - pr_no_cvd_history_mi))/(1)) * cvd_multiply)) * 1), # Conversion assumes constant rate over time

        # calculate annual non-CVD mortality rate (per person) as difference between all-cause mortality rate from CVD mortality rate in CDC WONDER database
        noncvd_mortality_rate  = rate,
        
        # Convert non-CVD rate to annual risk
        pr_noncvd_death        = 1 - exp((-noncvd_mortality_rate) * 1), # Conversion assumes constant rate over time
        
        # Equation can result in probabilities > 1 for very high ages; this puts a simple limit on it
        pr_noncvd_death        = ifelse(pr_noncvd_death < 1, 
                                        pr_noncvd_death, 1)                               
        
      ) %>%  
      # Keep only variables for transitions matrices
      select(-rate, -se,
             , -noncvd_mortality_rate
             , -ascvd_risk_10yr   
             , -ascvd_rate  
             , -ascvd_risk 
             , -monthly_mi_risk    
             , -monthly_mi_risk    
             , -annual_mi_rate     
             , -annual_mi_risk     
             , -annual_ascvd_mi_risk)
    
    return(mydata) 
    
  }

#### Cost function ----
# The Costs function estimates the costs at every cycle.

Costs <- function(mydata,
                  mi = NULL,
                  stroke = NULL,
                  heart_failure = NULL,
                  cardiac_dysrhythmia = NULL,
                  angina = NULL,
                  peripheral_artery_disease = NULL,
                  diabetes = NULL,
                  age_cat = NULL,
                  female = NULL,
                  race = NULL,
                  insurance = NULL,
                  fam_income = NULL,
                  education = NULL,
                  bmi_cat = NULL,
                  cci = NULL) {
  
  mydata  <- mydata %>%
    mutate(
      odds_nonzero_hc = (1.653 *
                           ifelse(coalesce({{ mi }}, mi) == 1, 2.937, 1) *
                           ifelse(coalesce({{ stroke }}, stroke) == 1, 0.948, 1) *
                           ifelse(coalesce({{ heart_failure }}, heart_failure) == 1, 1.830, 1) *
                           ifelse(coalesce({{ cardiac_dysrhythmia }}, cardiac_dysrhythmia) == 1, 2.389, 1) *
                           ifelse(coalesce({{ angina }}, angina) == 1, 1.472, 1) *
                           ifelse(coalesce({{ peripheral_artery_disease }}, peripheral_artery_disease) == 1, 1.334, 1) *
                           ifelse(coalesce({{ diabetes }}, diabetes) == 1, 4.382, 1) *
                           ifelse(coalesce({{ age_cat }}, age_cat) == "25-44", 1.171, 1) *
                           ifelse(coalesce({{ age_cat }}, age_cat) == "45-64", 2.121, 1) *
                           ifelse(coalesce({{ age_cat }}, age_cat) == "65+", 3.728, 1) *
                           ifelse(coalesce({{ female }}, female) == 1, 2.369, 1) *
                           ifelse(coalesce({{ race }}, race) == "white", 1.627, 1) *
                           ifelse(coalesce({{ race }}, race) == "black", 1.013, 1) *
                           ifelse(coalesce({{ race }}, race) == "asian", 1.025, 1) *
                           ifelse(coalesce({{ race }}, race) == "other_race", 1.461, 1) *
                           ifelse(coalesce({{ insurance }}, insurance) == "medicare", 1.714, 1) *
                           ifelse(coalesce({{ insurance }}, insurance) %in% c("medicaid", "other_plan"), 0.959, 1) *
                           ifelse(coalesce({{ insurance }}, insurance) == "uninsured", 0.355, 1) *
                           ifelse(coalesce({{ fam_income }}, fam_income) == "near_poor", 0.905, 1) *
                           ifelse(coalesce({{ fam_income }}, fam_income) == "low", 0.919, 1) *
                           ifelse(coalesce({{ fam_income }}, fam_income) == "medium", 0.960, 1) *
                           ifelse(coalesce({{ fam_income }}, fam_income) == "high", 1.344, 1) *
                           ifelse(coalesce({{ education }}, education) == "ged_hs", 1.178, 1) *
                           ifelse(coalesce({{ education }}, education) == "associate_bachelor", 1.726, 1) *
                           ifelse(coalesce({{ education }}, education) == "master_doctorate", 2.124, 1) *
                           ifelse(coalesce({{ bmi_cat }}, bmi_cat) == "normal_weight", 0.985, 1) *
                           ifelse(coalesce({{ bmi_cat }}, bmi_cat) == "overweight", 1.015, 1) *
                           ifelse(coalesce({{ bmi_cat }}, bmi_cat) == "obese", 1.233, 1) *
                           ifelse(coalesce({{ cci }}, cci) == 1, 2.037, 1) *
                           ifelse(coalesce({{ cci }}, cci) == 2, 3.920, 1) *
                           ifelse(coalesce({{ cci }}, cci) > 2, 8.287, 1)
      ),
      
      prob_nonzero_hc = odds_nonzero_hc / (1 + odds_nonzero_hc),
      
      nonzero_hc_cost = (2600.846 *
                           ifelse(coalesce({{ mi }}, mi) == 1, 1.177, 1) *
                           ifelse(coalesce({{ stroke }}, stroke) == 1, 1.016, 1) *
                           ifelse(coalesce({{ heart_failure }}, heart_failure) == 1, 0.948, 1) *
                           ifelse(coalesce({{ cardiac_dysrhythmia }}, cardiac_dysrhythmia) == 1, 1.449, 1) *
                           ifelse(coalesce({{ angina }}, angina) == 1, 1.442, 1) *
                           ifelse(coalesce({{ peripheral_artery_disease }}, peripheral_artery_disease) == 1, 1.425, 1) *
                           ifelse(coalesce({{ diabetes }}, diabetes) == 1, 1.676, 1) *
                           ifelse(coalesce({{ age_cat }}, age_cat) == "25-44", 1.403, 1) *
                           ifelse(coalesce({{ age_cat }}, age_cat) == "45-64", 2.104, 1) *
                           ifelse(coalesce({{ age_cat }}, age_cat) == "65+", 2.240, 1) *
                           ifelse(coalesce({{ female }}, female) == 1, 1.145, 1) *
                           ifelse(coalesce({{ race }}, race) == "white", 1.180, 1) *
                           ifelse(coalesce({{ race }}, race) == "black", 1.079, 1) *
                           ifelse(coalesce({{ race }}, race) == "asian", 1.034, 1) *
                           ifelse(coalesce({{ race }}, race) == "other_race", 1.224, 1) *
                           ifelse(coalesce({{ insurance }}, insurance) == "medicare", 1.122, 1) *
                           ifelse(coalesce({{ insurance }}, insurance) %in% c("medicaid", "other_plan"), 1.308, 1) *
                           ifelse(coalesce({{ insurance }}, insurance) == "uninsured", 0.619, 1) *
                           ifelse(coalesce({{ fam_income }}, fam_income) == "near_poor", 0.958, 1) *
                           ifelse(coalesce({{ fam_income }}, fam_income) == "low", 1.020, 1) *
                           ifelse(coalesce({{ fam_income }}, fam_income) == "medium", 0.919, 1) *
                           ifelse(coalesce({{ fam_income }}, fam_income) == "high", 0.960, 1) *
                           ifelse(coalesce({{ education }}, education) == "ged_hs", 1.076, 1) *
                           ifelse(coalesce({{ education }}, education) == "associate_bachelor", 1.114, 1) *
                           ifelse(coalesce({{ education }}, education) == "master_doctorate", 1.185, 1) *
                           ifelse(coalesce({{ bmi_cat }}, bmi_cat) == "normal_weight", 0.822, 1) *
                           ifelse(coalesce({{ bmi_cat }}, bmi_cat) == "overweight", 0.794, 1) *
                           ifelse(coalesce({{ bmi_cat }}, bmi_cat) == "obese", 0.922, 1) *
                           ifelse(coalesce({{ cci }}, cci) == 1, 1.558, 1) *
                           ifelse(coalesce({{ cci }}, cci) == 2, 2.582, 1) *
                           ifelse(coalesce({{ cci }}, cci) > 2, 2.991, 1)
      ),
      
      hc_cost = prob_nonzero_hc * nonzero_hc_cost
    )
  
  return(mydata)
}


#### Utilities function ---- 
# The Effs function to update the utilities at every cycle.

Effs <- function(mydata,
                 cardiac_dysrhythmia = NULL,
                 peripheral_artery_disease = NULL,
                 cci = NULL,
                 mi = NULL,
                 stroke = NULL,
                 heart_failure = NULL,
                 angina = NULL,
                 diabetes = NULL,
                 age_cat = NULL,
                 female = NULL,
                 race = NULL,
                 insurance = NULL,
                 fam_income = NULL,
                 education = NULL,
                 bmi_cat = NULL) {
  
  mydata <- mydata %>%
    mutate(
      utility = 0.817 +
        coalesce({{ mi }}, mi)                     * (-0.012) +
        coalesce({{ stroke }}, stroke)             * (-0.026) +
        coalesce({{ heart_failure }}, heart_failure) * (-0.045) +
        coalesce({{ cardiac_dysrhythmia }}, cardiac_dysrhythmia) * (-0.024) +
        coalesce({{ angina }}, angina)             * (-0.051) +
        coalesce({{ peripheral_artery_disease }}, peripheral_artery_disease) * (-0.035) +
        coalesce({{ diabetes }}, diabetes)         * (-0.033) +
        
        (coalesce({{ age_cat }}, age_cat) == "<25")     * 0 +
        (coalesce({{ age_cat }}, age_cat) == "25-44")   * (-0.026) +
        (coalesce({{ age_cat }}, age_cat) == "45-64")   * (-0.047) +
        (coalesce({{ age_cat }}, age_cat) == "65+")     * (-0.052) +
        
        coalesce({{ female }}, female) * (-0.024) +
        
        (coalesce({{ race }}, race) == "white")         * (-0.021) +
        (coalesce({{ race }}, race) == "black")         * (-0.003) +
        (coalesce({{ race }}, race) == "asian")         * (-0.008) +
        (coalesce({{ race }}, race) == "other_race")    * (-0.030) +
        
        (coalesce({{ insurance }}, insurance) == "medicare")       * (-0.007) +
        (coalesce({{ insurance }}, insurance) %in% c("medicaid", "other_plan")) * (-0.066) +
        (coalesce({{ insurance }}, insurance) == "uninsured")      * (-0.016) +
        
        (coalesce({{ fam_income }}, fam_income) == "near_poor")   * (0.009) +
        (coalesce({{ fam_income }}, fam_income) == "low")         * (0.023) +
        (coalesce({{ fam_income }}, fam_income) == "medium")      * (0.040) +
        (coalesce({{ fam_income }}, fam_income) == "high")        * (0.063) +
        
        (coalesce({{ education }}, education) == "ged_hs")             * (0.013) +
        (coalesce({{ education }}, education) == "associate_bachelor") * (0.016) +
        (coalesce({{ education }}, education) == "master_doctorate")   * (0.022) +
        
        (coalesce({{ bmi_cat }}, bmi_cat) == "normal_weight")      * (0.007) +
        (coalesce({{ bmi_cat }}, bmi_cat) == "overweight")         * (0.006) +
        (coalesce({{ bmi_cat }}, bmi_cat) == "obese")              * (-0.015) +
        
        (coalesce({{ cci }}, cci) == 1)        * (-0.039) +
        (coalesce({{ cci }}, cci) == 2)        * (-0.052) +
        (coalesce({{ cci }}, cci) > 2)         * (-0.076)
    )
  
  return(mydata)
}


#test <- Effs(mydata = param_data,
#             cci = 0,
#             cardiac_dysrhythmia = 0,      
#             peripheral_artery_disease = 0,
#             heart_failure = 1,
#             female = 1,
#             age_cat = "25-44",
#             race = "black",
#             insurance = "private",
#             fam_income = "high",
#             education = "associate_bachelor",
#             bmi_cat = "obese"
#             )
