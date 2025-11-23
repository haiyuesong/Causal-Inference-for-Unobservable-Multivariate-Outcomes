library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(gtsummary)

## ------------- 1. Load data -------------
df_demo  <- read.csv("AMYLOID_FMRI_Cohort_Study_My_Table_20Nov2025.csv",
                     stringsAsFactors = FALSE)

df_amy   <- read.csv("AMYLOID_FMRI_Cohort_Study_UCBERKELEY_AMY_6MM_20Nov2025.csv",
                     stringsAsFactors = FALSE)

df_fmri  <- read.csv("AMYLOID_FMRI_Cohort_Study_Functional_MRI_Images_20Nov2025.csv",
                     stringsAsFactors = FALSE)

## ------------- 2. Clean demographics (df_demo) -------------
df_demo_clean <- df_demo %>%
  arrange(subject_id, visit) %>%
  group_by(subject_id) %>%
  # keep gender, APOE, education and birth year
  tidyr::fill(PTGENDER, GENOTYPE, PTEDUCAT, PTDOBYY, .direction = "downup") %>%
  slice(1) %>% # only keep one record for each subject
  ungroup() %>%
  mutate(
    sex = ifelse(PTGENDER %in% c("Male","M","1"), "Male","Female"),
    # APOE4 indicator: 1 if any "4" allele present
    apoe4 = ifelse(str_detect(toupper(GENOTYPE), "4"), 1, 0)
  ) %>%
  select(subject_id, sex, education = PTEDUCAT, birth_year = PTDOBYY, apoe4)

## ------------- 3. Clean Amyloid PET table (df_amy) -------------
df_amy_clean <- df_amy %>%
  rename(subject_id = PTID) %>%
  mutate(
    # adjust date format
    PET_date = as.Date(SCANDATE), 
    amyloid_status = AMYLOID_STATUS
  ) %>%
  select(subject_id, PET_date, amyloid_status)

## ------------- 4. Clean fMRI info (df_fmri) -------------
df_fmri_clean <- df_fmri %>%
  mutate(
    fmri_date = as.Date(fmri_date),  # adjust format if needed
    rsfmri = ifelse(str_detect(fmri_description,"rsfMRI"),1,0)
  ) %>%
  filter(rsfmri == 1,
         fmri_n_images >= 100, # for a valid rsfMRI scan
         fmri_date >= as.Date("2017-01-01")) %>% # for ADNI 3 and 4
  rename(fmri_visit = fmri_visit) %>%
  # only keep subject_id, fMRI date, image_id, and visit_label
  select(subject_id, fmri_visit, fmri_date, image_id)

## ------------- 5. Match PET and fMRI by scan dates -------------
df_amy_fmri <- df_amy_clean %>%
  inner_join(df_fmri_clean, by = "subject_id") %>% # All subject-wise PET–fMRI combinations
  mutate(
    diff_days = as.numeric(fmri_date - PET_date)
  ) %>%
  filter(diff_days == 0) %>%
  dplyr::select(-diff_days)

## ------------- 6. Merge with demographics and compute age at fMRI/PET scan -------------
df_merged <- df_amy_fmri %>%
  left_join(df_demo_clean, by = "subject_id") %>%
  mutate(
    age = year(fmri_date) - birth_year
  )

## ------------- 7. Final tidy cohort dataset -------------
df_final <- df_merged %>%
  select(
    subject_id,
    fmri_visit,
    fmri_date,
    age,
    sex,
    education,
    apoe4,
    amyloid_status,
    image_id
  ) %>%
  group_by(subject_id) %>%
  slice_max(fmri_date, with_ties = FALSE) %>%
  ungroup() %>%
  select(subject_id, age, sex, education, apoe4, amyloid_status, image_id,fmri_date)

View(df_final)
write.csv(df_final,"cohort.csv")


## ------------- 8. Summary table by amyloid status -------------
df_final_tab <- df_final %>%
  mutate(
    amyloid_status = factor(amyloid_status,
                            levels = c(0, 1),
                            labels = c("Amyloid negative", "Amyloid positive")),
    sex = factor(sex),
    apoe4 = factor(apoe4,
                   levels = c(0, 1),
                   labels = c("APOE4 non-carrier", "APOE4 carrier"))
  )
tbl_amyloid <- df_final_tab %>%
  tbl_summary(
    by = amyloid_status,
    include = c(age, sex, education, apoe4),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no"
  ) %>%
  add_p() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Baseline characteristics by amyloid status**")

tbl_amyloid
