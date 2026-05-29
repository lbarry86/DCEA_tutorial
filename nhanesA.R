#------------------------------------------------------------------------------#
#---------------------Purpose: MSM DCEA Medicaid-CVD---------------------------#
#----------------------------------1-LOAD--------------------------------------#
#------------------------------------------------------------------------------#

#### Use nhanesA package to create NHANES dataset ----

# Define a function to use different variable names across different nhanes waves
get_variable <- function(year_letter) {
  # To account for different race variables used in different waves, i.e. separating "Asian" from "Other Race"
  RIDRETH <- if (year_letter %in% c("G", "H", "I", "J")) {
    "RIDRETH3"
  } else {
    "RIDRETH1"
  }
  
  # To account for different health insurance cover (any) variable name used in different waves
  H_INS_ANY <- if (year_letter %in% c("A", "B", "C")) {
    "HID010"
  } else {
    "HIQ011"
  }
  
  # To account for different health insurance cover (private) variables name used in different waves
  H_INS_PRIV <- if (year_letter %in% c("A", "B", "C")) {
    "HID030A"
  } else {
    "HIQ031A"
  }
  
  # To account for different health insurance cover (medicare) variables name used in different waves
  H_INS_MCARE <- if (year_letter %in% c("A", "B", "C")) {
    "HID030B"
  } else {
    "HIQ031B"
  }
  
  # To account for different health insurance cover (medicaid) variables name used in different waves
  H_INS_MCAID <- if (year_letter %in% c("A", "B", "C")) {
    "HID030C"
  } else {
    "HIQ031D"
  }
  
  # To account for different cholesterol variables name used in different waves
  HDL_CHOL <- if (year_letter %in% c("C")) {
    "LBXHDD"
  } else if (year_letter %in% c("B")) {
    "LBDHDL"
  } else {
    "LBDHDD"
  }
  
  # To account for different pregnancy variables name used in different waves
  PREG <- if (year_letter %in% c("A", "B")) {
    "RHQ141" # Do you think you are pregnant?
  } else {
    "RHD143" # Are you are pregnant?
  }
  
  # Store the INMB for this perturbation
  name_list <- list(RIDRETH     = RIDRETH,
                    H_INS_ANY   = H_INS_ANY,
                    H_INS_PRIV  = H_INS_PRIV,
                    H_INS_MCARE = H_INS_MCARE,
                    H_INS_MCAID = H_INS_MCAID,
                    HDL_CHOL    = HDL_CHOL,
                    PREG        = PREG)
  
  return(name_list)
}

# Define a function to get the year corresponding to a letter
get_year_from_letter <- function(letter) {
  letter_values <- c('B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J')
  year_values <- seq(2001, 2017, by = 2)
  return(year_values[letter == letter_values])
}

# Define a function to get the letter corresponding to a year
get_letter_from_year <- function(year) {
  letter_values <- c('B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J')
  year_values <- seq(2001, 2017, by = 2)
  return(letter_values[year == year_values])
}

# Function for creating NHANES pooled data
create_nhanes_data <- function(year_letter) {
  # List of NHANES datasets to merge
  merge_data <- list()
  
  # List of NHANES datasets to merge
  prefixes <- c('HIQ', 'INQ', 'MCQ', 'BPQ', 'SMQ', 'DIQ', 'RHQ', 'BPX', 'BMX', 'HDL', 'TCHOL', 'GLU', 'GHB', 'L13', 'L10AM', 'L10')
  # Update list according to wave being merged
  for (prefix in prefixes) {
    dataset_name <- paste0(prefix, '_', year_letter)
    merge_data[[dataset_name]] <- nhanes(dataset_name, includelabels = TRUE, translate=FALSE)
  }
  
  # Fetch the DEMO dataset (always available)
  nhanes_data <- nhanes(paste0('DEMO_', year_letter), includelabels = TRUE, translate=FALSE)
  
  # Merge the datasets (only the ones that exist)
  for (dataframe in merge_data) {
    if (!is.null(dataframe)) {
      nhanes_data <- left_join(nhanes_data, dataframe, by = "SEQN")
    }
  }  
  
  # Update variable names according to wave changes
  name_list <- get_variable(year_letter)
  
  # Select variables to include
  nhanes_data <- nhanes_data %>%
    select('SEQN', 'RIDSTATR',                  
           'RIDAGEYR', 'RIAGENDR', !!name_list$RIDRETH, 
           'DMDEDUC2', 'WTINT2YR', 'WTMEC2YR',                                       # from Table DEMO
           !!name_list$H_INS_ANY, !!name_list$H_INS_PRIV, 
           !!name_list$H_INS_MCARE, !!name_list$H_INS_MCAID,                         # from Table HIQ
           'INDFMPIR',                                                               # from Table INQ
           'DIQ010',                                                                 # from Table DIQ
           'MCQ160B', 'MCQ160C', 'MCQ160D', 'MCQ160E', 'MCQ160F',                    # from Table MCQ
           'BPQ020', 'BPQ040A', 'BPQ050A',                                           # from Table BPQ
           'SMQ020', 'SMQ040',                                                       # from Table SMQ
           !!name_list$PREG,                                                         # from Table RHQ
           'BPXSY1', 'BPXSY2', 'BPXSY3', 'BPXSY4',                                   # from Table BPX
           'BMXBMI',                                                                 # from Table BMX
           !!name_list$HDL_CHOL,                                                     # from Table HDL or I13
           'LBXTC',                                                                  # from Table TCHOL
           'LBXGLU', 'WTSAF2YR',                                                     # from Table GLU
           'LBXGH',                                                                  # from Table GHB
           'SDMVPSU', 'SDMVSTRA'
    ) 
  
  # Age in years; truncates at 80 years 
  nhanes_data <- nhanes_data %>%
    mutate(age_years = ifelse(RIDAGEYR == ".", NA, RIDAGEYR))
  
  # Gender        
  nhanes_data <- nhanes_data %>%
    mutate(male = case_when(RIAGENDR == 1 ~ 1, RIAGENDR == 2 ~ 0, TRUE ~ NA))
  
  # Race; Need to combine 'Mexican American' with 'Other Hispanic'
  nhanes_data <- nhanes_data %>%
    mutate(across(!!name_list$RIDRETH, 
                  ~ case_when(
                    (. == 1 | . == 2) ~ 1,  # Mexican American OR Other Hispanic
                    . == 3 ~ 2,             # Non-Hispanic White
                    . == 4 ~ 3,             # Non-Hispanic Black
                    . == 5 ~ 10,            # Other Race - Including Multi-Racial (2009 and earlier; includes Asian)
                    . == 6 ~ 4,             # Non-Hispanic Asian 
                    . == 7 ~ 5,             # Other Race - Including Multi-Racial (2011 and later; excludes Asian) 
                    TRUE ~ NA
                  ),
                  .names = "race_5cat"))
  
  # Highest education
  nhanes_data <- nhanes_data %>%
    mutate(edu_lev_5cat = case_when(DMDEDUC2 < 6 ~ DMDEDUC2,
                                    (DMDEDUC2 == 7 | DMDEDUC2 == 9) ~ NA,
                                    TRUE ~ NA))
  
  # ratio of annual family income to HHS poverty level
  nhanes_data <- nhanes_data %>%
    mutate(ratio_poverty_inc = ifelse(INDFMPIR <= 5,
                                      INDFMPIR, 
                                      NA))
  
  # Any health insurance
  nhanes_data <- nhanes_data %>%
    mutate(across(!!name_list$H_INS_ANY, 
                  ~ case_when(
                    . == 1 ~ 1,  # Covered
                    . == 2 ~ 0,  # Not covered
                    TRUE ~ NA
                  ),
                  .names = "health_ins_cover"))
  
  # private insurance
  nhanes_data <- nhanes_data %>%
    mutate(across(!!name_list$H_INS_PRIV, 
                  ~ case_when(
                    (. == 14 | . == 1) ~ 1,  # has private insurance
                    TRUE ~ NA
                  ),
                  .names = "health_ins_private"))
  
  # Medicare 
  nhanes_data <- nhanes_data %>%
    mutate(across(!!name_list$H_INS_MCARE, 
                  ~ case_when(
                    (. == 15 | . == 1) ~ 1,  # has medicare
                    TRUE ~ NA
                  ),
                  .names = "medicare"))
  
  # Medicaid 
  nhanes_data <- nhanes_data %>%
    mutate(across(!!name_list$H_INS_MCAID, 
                  ~ case_when(
                    (. == 17 | . == 1) ~ 1,  # has medicaid
                    TRUE ~ NA
                  ),
                  .names = "medicaid"))
  
  # Pregnant 
  nhanes_data <- nhanes_data %>%
    mutate(across(!!name_list$PREG, 
                  ~ case_when(
                    . == 1 ~ 1,  # pregnant
                    TRUE ~ NA
                  ),
                  .names = "pregnant"))
  
  # diagnosed diabetes                 
  nhanes_data <- nhanes_data %>%
    mutate(diabetes_dx = case_when(DIQ010 == 1 ~ 1,
                                   DIQ010 == 2 ~ 0,
                                   TRUE ~ NA))
  
  # congestive heart failure                 
  nhanes_data <- nhanes_data %>%
    mutate(cong_heart_fail_dx = case_when(MCQ160B == 1 ~ 1,
                                          MCQ160B == 2 ~ 0,
                                          TRUE ~ NA))
  
  # coronary heart disease
  nhanes_data <- nhanes_data %>%
    mutate(chd_dx = case_when(MCQ160C == 1 ~ 1,
                              MCQ160C == 2 ~ 0,
                              TRUE ~ NA))
  
  # angina pectoris
  nhanes_data <- nhanes_data %>%
    mutate(angina_dx = case_when(MCQ160D == 1 ~ 1,
                                 MCQ160D == 2 ~ 0,
                                 TRUE ~ NA))
  # myocardial infarction
  nhanes_data <- nhanes_data %>%
    mutate(heart_attack_dx = case_when(MCQ160E == 1 ~ 1,
                                       MCQ160E == 2 ~ 0,
                                       TRUE ~ NA))
  
  # stroke
  nhanes_data <- nhanes_data %>%
    mutate(stroke_dx = case_when(MCQ160F == 1 ~ 1,
                                 MCQ160F == 2 ~ 0,
                                 TRUE ~ NA))
  
  # high blood pressure
  nhanes_data <- nhanes_data %>%
    mutate(hyperten_dx = case_when(BPQ020 == 1 ~ 1,
                                     BPQ020 == 2 ~ 0,
                                     TRUE ~ NA))
  
  # hypertension treatment (currently)
  # Assumes:
  # those who self-report no doctor diagnosed high blood pressure (BPQ020 == 2) are not currently in receipt of hypertension medication
  # those who self-report no prescriptions for hypertension ever (BPQ040A == 2) are not currently in receipt of hypertension medication
  # those who self-report no prescriptions for hypertension currently (BPQ050A == 2) are not currently in receipt of hypertension medication
  nhanes_data <- nhanes_data %>%
    mutate(meds_hyperten = case_when(BPQ050A == 1 ~ 1,
                                     (BPQ020 == 2 | BPQ040A == 2 | BPQ050A == 2) ~ 0,
                                     TRUE ~ NA))
  
  # Combination of:
  # Ever smoked (100 cigarettes - used as threshold for 'never')
  # Current smoker
  nhanes_data <- nhanes_data %>%
    mutate(smoke_3cat = case_when(SMQ020 == 2 ~ 0,
                                  (SMQ040 == 1 | SMQ040 == 2) ~ 1,
                                  SMQ040 == 3 ~ 2,
                                  TRUE ~ NA))
  
  # Body Mass Index (kg/m**2)
  # https://www.cdc.gov/obesity/basics/adult-defining.html
  nhanes_data <- nhanes_data %>%
    mutate(bmi_cat = case_when(BMXBMI < 18.5 ~ 1,
                               (BMXBMI >= 18.5 & BMXBMI < 25) ~ 2,
                               (BMXBMI >= 25 & BMXBMI < 30) ~ 3,
                               BMXBMI >= 30 ~ 4,
                               TRUE ~ NA))
  
  # Direct HDL-Cholesterol (mg/dL)
  nhanes_data <- nhanes_data %>%
    mutate(across(!!name_list$HDL_CHOL, 
                  ~ ifelse(. != ".", ., NA),  # Replace non-missing values with themselves, and missing values with NA
                  .names = "hdl_chol_mgdl"))
  
  # Total Cholesterol (mg/dL)
  nhanes_data <- nhanes_data %>%
    mutate(total_chol_mgdl = ifelse(LBXTC == ".", 
                                    NA, 
                                    LBXTC))
  
  # Use this for full diabetes approach... https://www.cdc.gov/nchs/data/databriefs/db319.pdf
  # Fasting Glucose (mg/dL)
  nhanes_data <- nhanes_data %>%
    mutate(fast_glucose_mgdl = ifelse(LBXGLU == ".", 
                                      NA, 
                                      LBXGLU))
  
  # Glycohemoglobin (%); hba1c
  nhanes_data <- nhanes_data %>%
    mutate(hba1c = ifelse(LBXGH == ".",
                          NA, 
                          LBXGH))
  
  '
 # to compare with the nhanes 2015-2016 box folder dataset use this configuration of bp_systolic
    mutate(bp_systolic = ifelse(BPXSY1 == ".", 
                                NA, 
                                BPXSY1)) %>%
'
  
  # systolic blood pressure (mm Hg)
  # https://wwwn.cdc.gov/nchs/nhanes/1999-2000/BPX.htm
  nhanes_data <- nhanes_data %>%
    mutate(bp_systolic = case_when(any(is.na(BPXSY1), 
                                       is.na(BPXSY2), 
                                       is.na(BPXSY3)) ~ rowMeans(select(., starts_with("BPXSY")),
                                                                 na.rm = TRUE),
                                   TRUE ~ rowMeans(select(., starts_with("BPXSY"), -BPXSY4),
                                                   na.rm = TRUE)),
           bp_systolic = replace(bp_systolic, is.nan(bp_systolic), NA)) 
  
  return(nhanes_data)
}

### Create new nhanes dataset for each defined wave ###
for (letter in c('D', 'E', 'F', 'G', 'H', 'I')) {
  # Get the corresponding year
  year <- get_year_from_letter(letter)
  
  # Call the create_nhanes_data function and assign the result to a dataset with the year
  assign(paste0("nhanes_", year), create_nhanes_data(letter))
}

# Initialize an empty dataset to store the appended data
appended_nhanes_data <- NULL

# Loop through the years from 2001 to 2017 in steps of 2
for (year in seq(2001, 2017, by = 2)) {
  # Create the dataset name based on the year
  dataset_name <- paste0("nhanes_", year)
  
  # Check if the dataset exists
  if (exists(dataset_name)) {
    # Get the corresponding letter for the year
    letter <- get_letter_from_year(year)
    
    # Create a character variable for the wave year range
    wave_year <- paste0(letter, "_", substring(year, 3, 4))
    
    # Add the wave_year_range variable to the dataset
    dataset <- get(dataset_name) 
    dataset$wave_year <- wave_year
    dataset <- dataset %>%
      # select variables to keep
      select(SEQN,
             age_years, 
             stroke_dx, 
             heart_attack_dx, 
             chd_dx,
             cong_heart_fail_dx,
             diabetes_dx,
             angina_dx,
             hyperten_dx,
             pregnant,
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
             meds_hyperten,
             bp_systolic,
             hba1c,
             male,
             fast_glucose_mgdl,
             wave_year,
             WTINT2YR,
             WTSAF2YR,
             WTMEC2YR,
             SDMVPSU, 
             SDMVSTRA,
             RIDSTATR,
             BMXBMI) 
    
    # Append the dataset to the previously appended data
    if (is.null(appended_nhanes_data)) {
      appended_nhanes_data <- dataset
    } else {
      appended_nhanes_data <- rbind(appended_nhanes_data, dataset)
    }
  }
  
  # delete dataset once appended 
  rm(dataset_name)
}


### Save R file
# output data before applying exclusion criteria
write_rds(appended_nhanes_data, here("data", "output_data", "appended_nhanes_data.rds"))

