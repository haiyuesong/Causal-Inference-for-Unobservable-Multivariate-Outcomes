library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(gtsummary)

## ------------------------------------------------------------
## 1. Load data
## ------------------------------------------------------------
df_demo  <- read.csv("AMYLOID_FMRI_Cohort_Study_My_Table_20Nov2025.csv",
                     stringsAsFactors = FALSE)

df_amy   <- read.csv("AMYLOID_FMRI_Cohort_Study_UCBERKELEY_AMY_6MM_20Nov2025.csv",
                     stringsAsFactors = FALSE)

df_fmri  <- read.csv("AMYLOID_FMRI_Cohort_Study_Functional_MRI_Images_20Nov2025.csv",
                     stringsAsFactors = FALSE)

df_mri <- read.csv("ADNI3_Amyloid-Positive_fMRI_Key_MRI_04Dec2025.csv",
                   stringsAsFactors = FALSE)

## ------------------------------------------------------------
## 2. Clean DEMOGRAPHICS (df_demo)
##    - Fill PTEDUCAT, PTDOBYY, etc. within subject
##    - One row per subject
##    - Code sex and APOE4 indicator
## ------------------------------------------------------------
df_demo_clean <- df_demo %>%
  arrange(subject_id, visit) %>%
  group_by(subject_id) %>%
  tidyr::fill(PTGENDER, GENOTYPE, PTEDUCAT, PTDOBYY, .direction = "downup") %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    sex = ifelse(PTGENDER %in% c("Male","M","1"), "Male","Female"),
    # APOE4 indicator: 1 if any "4" allele present (e.g. "3/4", "4/4")
    apoe4 = ifelse(str_detect(toupper(GENOTYPE), "4"), 1, 0)
  ) %>%
  select(subject_id, sex, education = PTEDUCAT, birth_year = PTDOBYY, apoe4)

## ------------------------------------------------------------
## 3. Clean AMYLOID PET table (df_amy)
##    - We keep subject_id, PET scan date, amyloid_status
## ------------------------------------------------------------
df_amy_clean <- df_amy %>%
  rename(subject_id = PTID) %>%
  mutate(
    # adjust format if needed; ADNI is usually "YYYY-MM-DD"
    PET_date = as.Date(SCANDATE), 
    amyloid_status = AMYLOID_STATUS
  ) %>%
  select(subject_id, PET_date, amyloid_status)

## ------------------------------------------------------------
## 4. Clean fMRI info (df_fmri)
##    - We keep subject_id, fMRI date, image_id, and visit label
## ------------------------------------------------------------
df_fmri_clean <- df_fmri %>%
  mutate(
    fmri_date = as.Date(fmri_date),  # adjust format if needed
    rsfmri = ifelse(str_detect(fmri_description,"rsfMRI"),1,0)
  ) %>%
  filter(rsfmri == 1,
         fmri_n_images >= 100,
         fmri_date >= as.Date("2017-01-01")) %>% # for ADNI 3 and 4
  rename(fmri_visit = fmri_visit) %>%
  select(subject_id, fmri_visit, fmri_date, image_id)

## ------------------------------------------------------------
## 5. Clean MRI info (df_mri)
##    - We keep subject_id, MRI date, t1_image_id
##    - We keep Accelerated Sagittal MPRAGE MRI (T1w)
## ------------------------------------------------------------

df_mri_clean <- df_mri %>%
  filter(series_description == "Accelerated Sagittal MPRAGE") %>%
  mutate(
    mri_date = as.Date(image_date),
    t1_image_id = image_id
  ) %>%
  select(subject_id, mri_date, t1_image_id)

## ------------------------------------------------------------
## 6. Match PET and fMRI by scan dates
##    For each fMRI scan:
##    - Find PET scans from same subject
##    - PET_date == fmri_date
## ------------------------------------------------------------

# All subject-wise PET–fMRI combinations
df_pairs <- df_amy_clean %>%
  inner_join(df_fmri_clean, by = "subject_id") %>%
  mutate(
    diff_days = as.numeric(fmri_date - PET_date)
  ) %>%
  filter(diff_days == 0) %>%
  dplyr::select(-diff_days)


# For each fMRI scan (subject_id + image_id + fmri_date), keep same day PET
df_amy_fmri <- df_pairs %>%
  group_by(subject_id, image_id, fmri_date) %>%
  ungroup()

## ------------------------------------------------------------
## 6. Match MRI and fMRI by scan dates
##    For each fMRI scan:
##    - Find MRI scans from same subject
##    - mri_date == fmri_date
## ------------------------------------------------------------
df_mri_pairs <- df_mri_clean %>%
  inner_join(df_amy_fmri, by = "subject_id") %>%
  mutate(diff_days = as.numeric(fmri_date - mri_date)) %>%
  filter(diff_days == 0) %>%
  dplyr::select(-diff_days)

## ------------------------------------------------------------
## 7. Merge with demographics and compute age at fMRI
## ------------------------------------------------------------
df_merged <- df_mri_pairs %>%
  left_join(df_demo_clean, by = "subject_id") %>%
  mutate(
    age = year(fmri_date) - birth_year
  )

## ------------------------------------------------------------
## 7. Final tidy dataset
## ------------------------------------------------------------
df_final <- df_merged %>%
  select(
    subject_id,
    t1_image_id,
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
  rename(image_date = fmri_date) %>%
  ungroup() %>%
  # drop fmri_date if you don't want it in the final analysis table:
  select(subject_id, age, sex, education, apoe4, amyloid_status, image_id,t1_image_id,image_date)
View(df_final)
write.csv(df_final,"cohort.csv")

## ------------------------------------------------------------
## 8. Prepare variables for gtsummary
## ------------------------------------------------------------
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

## ------------------------------------------------------------
## 9. gtsummary table by amyloid status
## ------------------------------------------------------------
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
  modify_header(label ~ "**Characteristic**") %>%
  modify_caption("**Baseline characteristics by amyloid status**")

tbl_amyloid

## ------------------------------------------------------------
## 10. get the image IDs for image downloading
## ------------------------------------------------------------
paste(unique(df_final$image_id),collapse = ",")
paste(unique(df_final$t1_image_id),collapse = ",")
