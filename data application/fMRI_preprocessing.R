## ============================================================
## fMRI + MRI preprocessing for ADNI
## - For each fMRI run:
##   1) DICOM -> 4D NIfTI
##   2) Remove first N_discard timepoints
##   3) Slice timing correction (slicetimer)
##   4) Motion correction (mcflirt)
##   5) Trim to T_limit = 197 (for stability)
##   6) EPI ref: N4 + BET
##   7) T1 DICOM -> NIfTI + N4 + BET
##   8) EPI->T1 and T1->MNI (ANTS SyN)
##   9) Apply EPI->T1->MNI to full 4D fMRI (antsApplyTransforms -e 3)
##  10) Smoothing (sigma~2mm => FWHM~4.7mm)
##  11) AAL2 ROI extraction
##  12) Nuisance regression (Friston 24 motion)
##  13) Linear trend removal
##  14) Band-pass filter (0.01–0.1 Hz)
##  15) Weak stationary check
##  16) Save cleaned ROI ts as subjectid_imageid.csv
## ============================================================

## ==============================
## User settings
## ==============================

base_dir      <- "/home/data_application/ADNI"
template_file <- "/home/data_application/template/MNI152_T1_2mm.nii"
aal2_file     <- "/home/data_application/AAL2/AAL2.nii.gz"
aal2_txt_file <- "/home/data_application/AAL2/aal2_ROI_names.txt"

cohort_file   <- "/home/data_application/cohort.csv"

processed_out_dir <- "/home/data_application/processed_fmri_denoised"
dir.create(processed_out_dir, recursive = TRUE, showWarnings = FALSE)

N_discard <- 10   # remove first 10 volumes (DPARSF-style)
T_limit   <- 197  # max timepoints after discard
low_cut   <- 0.01 # Hz
high_cut  <- 0.10 # Hz

## ==============================
## Libraries
## ==============================
suppressPackageStartupMessages({
  library(oro.nifti)
  library(neurobase)
  library(signal) 
  library(tseries)
})

## ==============================
## Environment: threads
## ==============================
threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
if (!is.na(threads) && threads > 0) {
  Sys.setenv(ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS = as.character(threads))
}

## ==============================
## Helpers
## ==============================
run_cmd <- function(cmd, args = character()) {
  cat("Running:", cmd, paste(args, collapse = " "), "\n")
  res <- system2(cmd, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(res, "status")
  if (!is.null(status) && status != 0) {
    cat(res, sep = "\n")
    stop("Command failed: ", cmd)
  }
  invisible(res)
}

## Find candidate fMRI runs (I* dirs with DICOM)
find_fmri_runs <- function(subj_dir) {
  subj_dir <- normalizePath(subj_dir, mustWork = TRUE)
  all_I_dirs <- system2(
    "find",
    args = c(subj_dir, "-maxdepth", "6", "-type", "d", "-name", "I*"),
    stdout = TRUE
  )
  if (length(all_I_dirs) == 0) {
    stop("No I* DICOM directories found under ", subj_dir)
  }
  has_dicom <- function(d) {
    d <- normalizePath(d, mustWork = TRUE)
    n_files <- length(list.files(d, pattern = "\\.(dcm|ima|ima2)$",
                                 ignore.case = TRUE, recursive = TRUE))
    n_files > 0
  }
  good_dirs <- Filter(has_dicom, all_I_dirs)
  if (length(good_dirs) == 0) {
    stop("Found I* dirs but none contain DICOM files under ", subj_dir)
  }
  fmri_like <- grep("rfMRI|rsfMRI|rest|eyes", good_dirs,
                    value = TRUE, ignore.case = TRUE)
  if (length(fmri_like) > 0) {
    return(unique(normalizePath(fmri_like)))
  } else {
    warning("No obvious rsfMRI pattern; using all DICOM-containing I* dirs.")
    return(unique(normalizePath(good_dirs)))
  }
}

## Find DICOM dir for a given image_id (for T1)
find_image_dir <- function(subj_dir, image_id) {
  subj_dir <- normalizePath(subj_dir, mustWork = TRUE)
  img_name <- if (grepl("^I", image_id)) image_id else paste0("I", image_id)
  cat("  Searching for DICOM dir of image ID", img_name, "\n")
  out <- system2(
    "find",
    args = c(subj_dir, "-maxdepth", "6", "-type", "d", "-name", img_name),
    stdout = TRUE
  )
  if (length(out) == 0) {
    stop("No directory found for image ID ", img_name, " under ", subj_dir)
  }
  if (length(out) > 1) {
    warning("Multiple directories found for ", img_name, "; using first:\n",
            paste(out, collapse = "\n"))
  }
  normalizePath(out[1], mustWork = TRUE)
}

## Motion nuisance regression (Friston 24)
regress_motion <- function(Y, motion_par) {
  Y <- as.matrix(Y)
  motion <- as.matrix(motion_par)
  motion_deriv    <- rbind(motion[1, ], diff(motion))
  motion_sq       <- motion^2
  motion_deriv_sq <- motion_deriv^2
  C_motion <- cbind(motion, motion_deriv, motion_sq, motion_deriv_sq)
  C_motion <- scale(C_motion, center = TRUE, scale = TRUE)
  
  Tlen <- nrow(Y); p <- ncol(Y)
  Y_resid <- matrix(NA_real_, Tlen, p)
  for (j in seq_len(p)) {
    yj <- Y[, j]
    if (all(is.na(yj))) {
      Y_resid[, j] <- NA_real_
    } else {
      fit <- lm(yj ~ C_motion, na.action = na.exclude)
      Y_resid[, j] <- resid(fit)
    }
  }
  colnames(Y_resid) <- colnames(Y)
  Y_resid
}

## Linear detrend (per ROI time series)
detrend_ts <- function(Y) {
  Y <- as.matrix(Y)
  Tlen <- nrow(Y)
  tvec <- seq_len(Tlen)
  Y_detr <- apply(Y, 2, function(y) {
    if (all(is.na(y))) return(y)
    resid(lm(y ~ tvec))
  })
  colnames(Y_detr) <- colnames(Y)
  Y_detr
}

## Bandpass filter using TR
bandpass_filter <- function(Y, TR, low = 0.01, high = 0.10) {
  Y <- as.matrix(Y)
  fs <- 1 / TR
  nyq <- fs / 2
  Wn  <- c(low, high) / nyq
  if (any(Wn <= 0) || any(Wn >= 1)) {
    warning("Invalid bandpass frequencies for TR=", TR,
            "; skipping temporal filtering.")
    return(Y)
  }
  bf <- butter(2, Wn, type = "pass")
  Y_filt <- apply(Y, 2, function(x) {
    if (all(is.na(x))) return(x)
    signal::filtfilt(bf, x)
  })
  colnames(Y_filt) <- colnames(Y)
  Y_filt
}
## Weak stationarity check: ALL ROIs must pass ADF (p < alpha)
check_stationarity_all_rois <- function(Y, alpha = 0.05) {
  Y <- as.matrix(Y)
  p <- ncol(Y)
  if (p == 0) {
    cat("  No ROIs to test; treating as non-stationary.\n")
    return(FALSE)
  }
  
  pvals <- numeric(p)
  for (j in seq_len(p)) {
    y <- Y[, j]
    if (all(is.na(y))) {
      pvals[j] <- 1   # definitely non-stationary
      next
    }
    out <- tryCatch(
      tseries::adf.test(y, k = 1),   # k=1 lag, reasonable after bandpass
      error = function(e) NULL
    )
    if (is.null(out)) {
      pvals[j] <- 1   # treat errors as non-stationary
    } else {
      pvals[j] <- out$p.value
    }
  }
  
  n_stationary <- sum(pvals < alpha, na.rm = TRUE)
  cat("  ADF stationarity check: ", n_stationary, "/", p,
      " ROIs stationary (p < ", alpha, ")\n", sep = "")
  
  all(pvals < alpha, na.rm = TRUE)
}

## ==============================
## Load cohort once
## ==============================
cohort <- read.csv(cohort_file, stringsAsFactors = FALSE)
if (!all(c("subject_id", "image_id", "t1_image_id") %in% names(cohort))) {
  stop("cohort.csv must have columns: subject_id, image_id, t1_image_id")
}

## ==============================
## Preprocess ONE subject
## ==============================
preprocess_subject <- function(subject_id) {
  cat("\n============================================\n")
  cat("Processing subject:", subject_id, "\n\n")
  
  subj_dir <- file.path(base_dir, subject_id)
  subj_dir <- normalizePath(subj_dir, mustWork = TRUE)
  
  cohort_subj <- cohort[cohort$subject_id == subject_id, , drop = FALSE]
  
  fmri_dirs <- find_fmri_runs(subj_dir)
  cat("Found", length(fmri_dirs), "candidate fMRI runs\n")
  
  for (d in fmri_dirs) {
    cat("\n--------------------------------------------\n")
    cat("Run:", d, "\n")
    
    image_id <- basename(d)              # e.g. "I1422666"
    image_id_clean <- sub("^I", "", image_id)
    cat("  fMRI image_id_clean:", image_id_clean, "\n")
    
    # Find corresponding T1 in cohort
    row_idx <- which(cohort_subj$image_id == image_id_clean)
    if (length(row_idx) == 0) {
      warning("No matching row in cohort for subject=", subject_id,
              " fMRI image_id=", image_id_clean, "; skipping this run.")
      next
    }
    if (length(row_idx) > 1) {
      warning("Multiple cohort rows for subject ", subject_id,
              " fMRI image_id=", image_id_clean, "; using first.")
      row_idx <- row_idx[1]
    }
    mri_image_id <- as.character(cohort_subj$t1_image_id[row_idx])
    if (is.na(mri_image_id) || mri_image_id == "") {
      warning("t1_image_id missing for subject ", subject_id,
              " fMRI image_id=", image_id_clean, "; skipping this run.")
      next
    }
    cat("  Matched T1 image_id:", mri_image_id, "\n")
    
    out_root <- file.path(subj_dir, paste0("proc_", image_id_clean))
    dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
    out_root <- normalizePath(out_root, mustWork = TRUE)
    setwd(out_root)
    
    fmri_dicom_dir <- d
    mri_dicom_dir  <- find_image_dir(subj_dir, mri_image_id)
    
    cat("  fMRI dicom:", fmri_dicom_dir, "\n")
    cat("  T1   dicom:", mri_dicom_dir,  "\n")
    
    ## 1. DICOM → NIfTI
    cat("Converting fMRI DICOM → NIfTI with dcm2niix...\n")
    res <- system2("dcm2niix",
                   args = c("-z","y","-f","fmri","-o",out_root,fmri_dicom_dir),
                   stdout = TRUE, stderr = TRUE)
    if (!file.exists("fmri.nii.gz")) {
      warning("No fmri.nii.gz for run ", fmri_dicom_dir, "; skipping.")
      next
    }
    fmri_img <- readNIfTI("fmri.nii.gz", reorient = FALSE)
    dims <- dim(fmri_img)
    if (length(dims) != 4L) {
      warning("fmri.nii.gz is not 4D; skipping.")
      next
    }
    ntime  <- dims[4]
    nslice <- dims[3]
    TR     <- fmri_img@"pixdim"[5]
    cat("  fMRI dims:", paste(dims, collapse=" x "),
        " (T=", ntime, ", slices=", nslice, ", TR=", TR, "s)\n")
    if (ntime <= N_discard + 20 || nslice < 20 || nslice > 200) {
      warning("  Too few timepoints or suspicious slices; skipping run ", d)
      next
    }
    
    ## 2. Remove first N_discard volumes
    cat("Removing first", N_discard, "timepoints (signal equilibration/adaptation)...\n")
    run_cmd("fslroi", c("fmri.nii.gz", "fmri_trim.nii.gz",
                        as.character(N_discard),
                        as.character(ntime - N_discard)))
    # Update ntime after discard
    fmri_trim_img <- readNIfTI("fmri_trim.nii.gz", reorient = FALSE)
    ntime_trim <- dim(fmri_trim_img)[4]
    cat("  Remaining timepoints:", ntime_trim, "\n")
    if (ntime_trim < 50) {
      warning("Too few timepoints after discard; skipping.")
      next
    }
    
    ## 3. Slice timing correction
    cat("Slice timing correction (slicetimer)...\n")
    # NOTE: adjust slicetimer options (e.g. --odd, --down) to match ADNI slice order
    run_cmd("slicetimer", c(
      "-i", "fmri_trim.nii.gz",
      "-o", "fmri_stc.nii.gz"
    ))
    
    ## 4. Motion correction (mcflirt)
    cat("Motion correction with mcflirt...\n")
    run_cmd("mcflirt", c(
      "-in",  "fmri_stc.nii.gz",
      "-out", "fmri_mc",
      "-refvol", "0",
      "-plots"
    ))
    
    ## 5. Trim to T_limit
    base <- "fmri_mc"
    if (ntime_trim > T_limit) {
      cat("Truncating motion-corrected series to first", T_limit, "timepoints...\n")
      run_cmd("fslroi", c("fmri_mc", "fmri_mc_trim", "0", as.character(T_limit)))
      base <- "fmri_mc_trim"
      ntime_trim <- T_limit
    }
    
    ## 6. Reference + BET + N4 (EPI)
    cat("Extracting and skull-stripping reference volume...\n")
    run_cmd("fslroi", c(base, "fmri_ref", "0", "1"))
    run_cmd("N4BiasFieldCorrection", c(
      "-d","3","-i","fmri_ref.nii.gz","-o","fmri_ref_N4.nii.gz"))
    run_cmd("bet", c("fmri_ref_N4.nii.gz", "fmri_ref_brain.nii.gz"))
    
    ## 7. T1 DICOM → NIfTI, N4, BET
    cat("Converting T1 DICOM → NIfTI with dcm2niix + N4 + BET...\n")
    res_mri <- system2("dcm2niix",
                       args = c("-z","y","-f","T1","-o",out_root, mri_dicom_dir),
                       stdout = TRUE, stderr = TRUE)
    if (!file.exists("T1.nii.gz")) {
      warning("No T1.nii.gz for run ", d, "; skipping.")
      next
    }
    T1_img <- readNIfTI("T1.nii.gz", reorient = FALSE)
    cat("  T1 dims:", paste(dim(T1_img), collapse=" x "), "\n")
    run_cmd("N4BiasFieldCorrection", c("-d","3","-i","T1.nii.gz","-o","T1_N4.nii.gz"))
    run_cmd("bet", c("T1_N4.nii.gz", "T1_brain.nii.gz"))
    
    ## 8. EPI→T1 registration
    cat("Registering EPI ref to T1 (EPI→T1)...\n")
    run_cmd("antsRegistrationSyN.sh", c(
      "-d","3",
      "-f","T1_brain.nii.gz",
      "-m","fmri_ref_brain.nii.gz",
      "-o","EPI2T1_"
    ))
    
    ## 9. T1→MNI registration
    cat("Registering T1 to MNI (T1→MNI)...\n")
    run_cmd("antsRegistrationSyN.sh", c(
      "-d","3",
      "-f", template_file,
      "-m","T1_brain.nii.gz",
      "-o","T12MNI_"
    ))
    
    ## 10. Apply transforms to full 4D fMRI
    cat("Applying EPI→T1→MNI transforms to 4D fMRI...\n")
    run_cmd("antsApplyTransforms", c(
      "-d","3","-e","3","--float",
      "-i", paste0(base, ".nii.gz"),
      "-o","fmri_mni.nii.gz",
      "-t","T12MNI_1Warp.nii.gz",
      "-t","T12MNI_0GenericAffine.mat",
      "-t","EPI2T1_1Warp.nii.gz",
      "-t","EPI2T1_0GenericAffine.mat",
      "-r", template_file
    ))
    
    ## 11. Smoothing
    cat("Smoothing in MNI space...\n")
    desired_fwhm <- 4.7
    sigma <- desired_fwhm / 2.3548
    cat("  Using sigma =", sigma, "for FWHM ≈", desired_fwhm, "mm\n")
    run_cmd("fslmaths", c(
      "fmri_mni.nii.gz",
      "-s", sprintf("%.3f", sigma),
      "fmri_mni_smooth.nii.gz"
    ))
    
    ## 12. ROI extraction
    cat("Extracting AAL2 ROI time series...\n")
    img4d <- readNIfTI("fmri_mni_smooth.nii.gz", reorient = FALSE)
    aal   <- readNIfTI(aal2_file,           reorient = FALSE)
    aal_info <- read.table(aal2_txt_file, stringsAsFactors = FALSE)
    
    if (!all(dim(aal)[1:3] == dim(img4d)[1:3])) {
      warning("AAL2 and fMRI dims mismatch; skipping run ", d)
      next
    }
    
    Tlen <- dim(img4d)[4]
    nroi <- max(aal_info$V1)
    ts_mat <- matrix(NA_real_, Tlen, nroi)
    
    for (k in seq_len(nroi)) {
      roi_label <- aal_info$V3[k]
      mask <- (aal == roi_label)
      if (!any(mask)) {
        warning("ROI ", k, " has no voxels; skipping")
        next
      }
      vox <- matrix(img4d[mask], nrow = sum(mask), ncol = Tlen)
      ts_mat[, k] <- colMeans(vox, na.rm = TRUE)
    }
    ts_df <- as.data.frame(ts_mat)
    colnames(ts_df) <- aal_info$V2
    
    ## 13. Motion nuisance regression
    cat("Regressing motion (Friston 24)...\n")
    motion_file <- file.path(out_root, "fmri_mc.par")
    if (!file.exists(motion_file)) {
      warning("No motion par file for ", d, "; skipping nuisance regression.")
      ts_clean1 <- ts_df
    } else {
      motion <- as.matrix(read.table(motion_file))
      if (nrow(motion) >= Tlen) {
        motion <- motion[1:Tlen, , drop = FALSE]
      } else {
        warning("Motion par rows < Tlen; recycling last row.")
        motion <- rbind(motion,
                        matrix(rep(motion[nrow(motion),], Tlen - nrow(motion)),
                               ncol = ncol(motion), byrow = TRUE))
      }
      ts_clean1 <- regress_motion(ts_df, motion)
    }
    
    
    ## 14. Linear trend removal
    cat("Removing linear trend per ROI time series...\n")
    ts_detr <- detrend_ts(ts_clean1)
    
    ## 15. Bandpass filter
    cat("Applying bandpass filter (", low_cut, "-", high_cut, "Hz)...\n")
    ts_filt <- bandpass_filter(ts_detr, TR, low = low_cut, high = high_cut)
    
    ts_final <- as.data.frame(ts_filt)
    colnames(ts_final) <- colnames(ts_df)
    
    ## 16. Stationarity check for ALL ROIs
    cat("Checking weak stationarity (ADF) for all ROIs...\n")
    is_stationary <- check_stationarity_all_rois(ts_final, alpha = 0.05)
    
    if (!is_stationary) {
      warning("Time series for subject ", subject_id,
              " run ", image_id_clean,
              " failed stationarity check; dropping this run.")
      next  # skip saving this run
    }
    
    ## 17. Save cleaned, stationary ROI time series
    out_name   <- paste0(subject_id, "_", image_id_clean, ".csv")
    ts_outfile <- file.path(processed_out_dir, out_name)
    write.csv(ts_final, ts_outfile, row.names = TRUE)
    
    cat("Saved cleaned, stationary ROI time series to:", ts_outfile, "\n")
    
  }
  
  cat("\nCompleted subject:", subject_id, "\n")
}

## ==============================
## SLURM ARRAY INPUT
## ==============================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: preprocess_fmri_mri_denoised_DPARSF.R <array_index>")
}

idx <- as.integer(args[1])

all_subjects <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
all_subjects <- all_subjects[nchar(all_subjects) > 0]
all_subjects <- sort(all_subjects)

if (idx > length(all_subjects)) {
  stop("Index ", idx, " exceeds number of subject folders (", length(all_subjects), ")")
}

subject_id <- all_subjects[idx]
cat("SLURM ARRAY INDEX:", idx, "→ subject:", subject_id, "\n")

preprocess_subject(subject_id)
