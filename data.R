#------------------------------------------------------------------------------#
#---------------------------Created by Luke Barry------------------------------#
#-----------------------------Date:07/11/2025----------------------------------#
#--------------------------Purpose: MSM CVD DCEA ------------------------------#
#-------------------------------Dataset Prep-----------------------------------#
#------------------------------------------------------------------------------#

#### Load the data ----
  
  nhanes1 <- read_rds(here("data", "output_data",
                          "appended_nhanes_data.rds")) %>%
    as_tibble() %>%
    select(age_years, 
           stroke_dx, 
           heart_attack_dx, 
           chd_dx,
           cong_heart_fail_dx,
           angina_dx,
           race_5cat, 
           ratio_poverty_inc, 
           edu_lev_5cat, 
           hdl_chol_mgdl,
           total_chol_mgdl,
           health_ins_cover,
           health_ins_private,
           medicare,
           medicaid,
           smoke_3cat,
           bmi_cat,
           hyperten_dx,
           meds_hyperten,
           bp_systolic,
           diabetes_dx,
           male, 
           wave_year) %>% 
    print(n = 10, width = Inf)
  
  # Import non-CVD death rates [taken from CDC WONDER]
  nonCVD_mort <- read.xlsx(here("data", "raw_data", 
                                "AllCause-MI_Deaths_Rates_2005-16.xlsx"), 
                           sheet = "ImportR_nCVD") %>%
    rename(age = Age) %>%
    mutate(female = ifelse(Gender == "Female", 1, ifelse(Gender == "Male", 0, NA))) %>%
    select(-Gender,
           -Gender.Code)
  
  saveRDS(nonCVD_mort, here("data", "output_data", "nonCVD_mort.rds"))
  
  #### Format the data ----

  nhanes2 <- nhanes1 %>% 
    # create categorical variables to correspond to MEPS cost and utility calculations
    rename(age = age_years,
           hdl_chol = hdl_chol_mgdl,
           total_chol = total_chol_mgdl,
           sys_bp = bp_systolic) %>%
    mutate(female = if_else(male == 1, 0, 1), 
           mi            = as.numeric(heart_attack_dx),
           stroke        = as.numeric(stroke_dx),
           chd           = as.numeric(chd_dx),
           heart_failure = as.numeric(cong_heart_fail_dx),
           angina        = as.numeric(angina_dx),
           diabetes      = as.numeric(diabetes_dx),
           fpl           = as.numeric(ratio_poverty_inc),
           hyperten_trt  = case_when(hyperten_dx == 0 ~ 0,
                                     hyperten_dx == 1 & meds_hyperten == 0 ~ 0,
                                     hyperten_dx == 1 & meds_hyperten == 1 ~ 1,
                                     TRUE ~ NA),
           insurance = case_when(health_ins_private == "1"  ~ "private",          
                                 medicare           == "1"  ~ "medicare",
                                 medicaid           == "1"  ~ "medicaid",
                                 health_ins_cover   == "1"  ~ "other_plan",
                                 health_ins_cover   == "0"  ~ "uninsured",
                                 TRUE                       ~ NA),
           fam_income = case_when(ratio_poverty_inc <= 1     ~ "poor",                          
                                  ratio_poverty_inc <= 1.25  ~ "near_poor",
                                  ratio_poverty_inc <= 2     ~ "low",
                                  ratio_poverty_inc <= 4     ~ "medium",
                                  ratio_poverty_inc  > 4     ~ "high",
                                  TRUE         ~ NA),
           fam_income = factor(fam_income,
                               levels = c("poor", "near_poor", "low", "medium", "high"),
                               ordered = TRUE),
           education = case_when((edu_lev_5cat == "1" | edu_lev_5cat == "2")  ~ "no_degree",              
                                 edu_lev_5cat == "3"  ~ "ged_hs",
                                 edu_lev_5cat == "4"  ~ "associate_bachelor",
                                 edu_lev_5cat == "5"  ~ "master_doctorate",
                                 TRUE                 ~ NA), 
           education = factor(education,
                              levels = c("no_degree", "ged_hs", "associate_bachelor", "master_doctorate"),
                              ordered = TRUE),
           race = case_when(race_5cat == "1"  ~ "hispanic",                       
                            race_5cat == "2"  ~ "white",
                            race_5cat == "3"  ~ "black",
                            race_5cat == "4"  ~ "asian",
                            race_5cat == "5"  ~ "other_race",
                            TRUE              ~ NA),   
           bmi_cat = case_when(bmi_cat == "1" ~ "underweight", 
                               bmi_cat == "2" ~ "normal_weight",
                               bmi_cat == "3" ~ "overweight",  
                               bmi_cat == "4" ~ "obese",
                               TRUE              ~ NA),
           smoking = case_when(smoke_3cat         == "1"  ~ 1,                
                               smoke_3cat         == "2"  ~ 0,
                               smoke_3cat         == "0"  ~ 0,
                               TRUE                       ~ NA),
           # create vector of starting states per individual according to their history of a diagnosis of any CVD conditions in bootstrapped NHANES data
           cvd = ifelse((mi == 1 | stroke == 1 | chd == 1 | heart_failure == 1),
                        "history_cvd", "no_cvd"),  
           # categorize age into specified categories each cycle that age is updated - 
           # Same categories used for cost and utility prediction as for weighting MI/stroke risk by ASCVD risk
           age_cat = cut(age, breaks = c(0, 24, 44, 64, Inf)
                         , labels = c('<25', '25-44', '45-64', '65+')),
           across(all_of(c("bmi_cat", "cvd", "race", "education"
                           , "fam_income", "insurance")), to_factor)
           ) %>%
    # drop unused variables after formatting
    select(-smoke_3cat
           , -health_ins_private
           , -medicare          
           , -medicaid          
           , -health_ins_cover  
           , -health_ins_cover
           , -race_5cat
           , -ratio_poverty_inc
           , -edu_lev_5cat
           , -hyperten_dx
           , -meds_hyperten
           , -heart_attack_dx
           , -stroke_dx
           , -chd_dx
           , -cong_heart_fail_dx
           , -diabetes_dx
           , -angina_dx
           , -male
           , -ratio_poverty_inc
           ) %>% 
    modify_if(is.labelled, to_factor) %>%
    set_variable_labels(age           = "Age (Years)",
                        age_cat       = "Age Category (Years)",
                        hdl_chol      = "Direct HDL-Cholesterol (mg/dL)",
                        total_chol    = "Total Cholesterol (mg/dL)",
                        bmi_cat       = "BMI Category",
                        sys_bp        = "Mean Systolic Blood Pressure (mm Hg)",
                        female        = "Female (Y/N)",
                        mi            = "Myocardial Infarction (Y/N)",
                        stroke        = "Stroke (Y/N)",
                        chd           = "Coronary Heart Disease (Y/N)",
                        heart_failure = "Heart Failure (Y/N)",
                        angina        = "Angina (Y/N)",
                        diabetes      = "Diabetes (Y/N)",
                        hyperten_trt  = "Hypertension Medication (Y/N)",
                        insurance     = "Health Insurance Status",
                        fam_income    = "Annual Family Income",
                        education     = "Highest Educational Achievement",
                        race          = "Race/Ethnicity",
                        smoking       = "Current Smoker (Y/N)",
                        cvd           = "History of CVD (Y/N)",
                        wave_year     = "NHANES wave and survey year",
                        fpl           = "Income to poverty ratio") 


  #glimpse(nhanes2)

  # Drop missing variables and select survey year for main data
  nhanes_sample <- nhanes2 %>%
    # drop missing variables [for simplicity - not recommended]
    na.omit() 

  saveRDS(nhanes_sample, here("data", "output_data", "nhanes_sample.rds"))

  