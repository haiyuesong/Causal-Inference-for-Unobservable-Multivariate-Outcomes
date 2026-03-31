## ============================================================
## compute_ec_subject.R
##   - Compute EC vector for one subject's fMRI time series
## ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript compute_ec_subject.R <file_index>", call. = FALSE)
}
idx <- as.integer(args[1])
if (is.na(idx) || idx <= 0) {
  stop("file_index must be a positive integer", call. = FALSE)
}

## ------------------------------------------------------------
## 0. Paths and basic settings
## ------------------------------------------------------------

base_dir    <- "/home/data_application"
trimmed_dir <- file.path(base_dir, "processed_fmri")     # input CSVs
ec_dir      <- file.path(base_dir, "ec_denoised")                    # output RDS
aal2_txt_file <- file.path(base_dir, "AAL2/aal2_ROI_names.txt")

dir.create(ec_dir, showWarnings = FALSE, recursive = TRUE)

setwd(base_dir)

## ------------------------------------------------------------
## 1. Libraries
## ------------------------------------------------------------
suppressPackageStartupMessages({
  library(parallel)
  library(vars)
  library(data.table)
  library(dplyr)
  library(psych)
  library(lmtest)   # for grangertest
})

## ------------------------------------------------------------
## 2. Fixed parameters
## ------------------------------------------------------------

aal2_txt <- read.table(aal2_txt_file, stringsAsFactors = FALSE)
nroi     <- 120
t_corr   <- 50               # correlation window length
std_err  <- 1 / sqrt(t_corr - 3)

# cores for mclapply
mc_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
if (is.na(mc_cores) || mc_cores <= 0) mc_cores <- 4
options(mc.cores = mc_cores)

## ------------------------------------------------------------
## 3. Helper functions
## ------------------------------------------------------------

# Fisher z-transform (psych::fisherz is fine; wrapper for clarity)
fisherz_ <- function(r) psych::fisherz(r)

# high correlation test
high_corr <- function(roi1, roi2, std_fisherz, threshold = 0.5, std_err) {
  # std_fisherz: standardized fisher-z corr matrix
  # H0: corr < threshold vs corr >= threshold (one-sided)
  std_fisherz_threshold <- fisherz_(threshold) / std_err
  z_diff <- std_fisherz[roi1, roi2] - std_fisherz_threshold
  p_value <- 1 - pnorm(z_diff)
  p_value
}

# robust scaling: center, and scale if sd>0; otherwise just center
scale_column <- function(x) {
  if (all(is.na(x))) return(x)
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(x - m)
  } else {
    return(as.numeric(scale(x)))
  }
}

# Safe 2-variable Granger F-stat
safe_granger_two_var <- function(y, x, order_p) {
  dat <- data.frame(y = y, x = x)
  out <- tryCatch(
    {
      gt <- lmtest::grangertest(y ~ x, order = order_p, data = dat)
      as.numeric(gt$`F`[2])
    },
    error = function(e) {
      NA_real_
    }
  )
  out
}

# Main conditional Granger function (robust)
intersection_conditional_granger <- function(roi_1, roi_2, ts, order_p,
                                             high_corr_mat,
                                             max_cond = 40) {
  # ts: T x nroi matrix (or ts object) after scaling and dropping first t_corr rows
  # high_corr_mat: nroi x nroi logical
  # max_cond: maximum number of conditioning ROIs
  
  if (roi_1 == roi_2) return(NA_real_)
  
  cond_roi <- intersect(which(high_corr_mat[, roi_1]),
                        which(high_corr_mat[, roi_2]))
  
  # cap the size of conditioning set
  if (length(cond_roi) > max_cond) {
    cond_roi <- cond_roi[1:max_cond]  # simple choice: first few; could rank by corr if desired
  }
  
  if (length(cond_roi) > 0) {
    # Conditional Granger: Y = roi_2, X = roi_1, plus cond ROIs
    vars_in <- c(roi_2, roi_1, cond_roi)
    dat <- ts[, vars_in, drop = FALSE]
    colnames(dat) <- paste0("V", seq_along(vars_in))  # V1=y, V2=x, rest=cond
    
    res <- tryCatch(
      {
        var_fit <- vars::VAR(dat, p = order_p, type = "const")
        cg      <- vars::causality(var_fit, cause = "V2")
        as.numeric(cg$Granger$statistic)
      },
      error = function(e) {
        # Non-full rank, singular, etc.
        NA_real_
      }
    )
    return(res)
    
  } else {
    # Unconditional (2-var) Granger
    y <- ts[, roi_2]
    x <- ts[, roi_1]
    return(safe_granger_two_var(y, x, order_p))
  }
}

## ------------------------------------------------------------
## 4. Select CSV file by index and parse subject_id & image_id
## ------------------------------------------------------------

csv_files <- list.files(trimmed_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) {
  stop("No CSV files found in trimmed_dir: ", trimmed_dir)
}
csv_files <- sort(csv_files)

if (idx > length(csv_files)) {
  stop("Index ", idx, " exceeds number of files (", length(csv_files), ")")
}

csv_file <- csv_files[idx]

# Strip extension and parse IDs from names like "012_S_6073_11344002.csv"
base_name <- tools::file_path_sans_ext(basename(csv_file))
parts <- strsplit(base_name, "_")[[1]]

if (length(parts) >= 4) {
  # ADNI-style: "012_S_6073_11344002" → subject_id = "012_S_6073", image_id = "11344002"
  subject_id <- paste(parts[1:3], collapse = "_")
  image_id   <- parts[4]
} else if (length(parts) == 2) {
  # Fallback: "SubjectID_ImageID"
  subject_id <- parts[0 + 1]
  image_id   <- parts[1 + 1]
} else {
  stop("Filename '", base_name, "' does not match expected pattern <subject>_<image>.csv")
}

cat("------------------------------------------------------------\n")
cat("EC computation for index:", idx, "\n")
cat("  file      :", csv_file, "\n")
cat("  subject ID:", subject_id, "\n")
cat("  image ID  :", image_id, "\n")
cat("------------------------------------------------------------\n")


## ------------------------------------------------------------
## 5. Read time series and prepare ROI data
## ------------------------------------------------------------

data <- fread(csv_file)

# Assume first column is time index; next 120 columns are ROIs
if (ncol(data) < (nroi + 1)) {
  stop("File ", csv_file, " has fewer than ", nroi + 1, " columns; cannot extract 120 ROIs.")
}

roi_data <- as.matrix(data[, 2:(nroi + 1), with = FALSE])
T_total  <- nrow(roi_data)

cat("Total time points in file:", T_total, "\n")

if (T_total <= t_corr + 10) {
  stop("Not enough time points in ", csv_file, " for t_corr=", t_corr, " & VAR fitting.")
}

## ------------------------------------------------------------
## 6. Correlation and high-corr mask on first t_corr points
## ------------------------------------------------------------

cat("Computing correlation matrix on first", t_corr, "time points...\n")

cor_mat <- cor(roi_data[1:t_corr, ], use = "pairwise.complete.obs")
fisherz_mat <- fisherz_(cor_mat)
standard_fisherz <- fisherz_mat / std_err

high_corr_mat <- matrix(FALSE, nrow = nroi, ncol = nroi)
for (roi1 in 1:(nroi - 1)) {
  p_vec <- unlist(mclapply((roi1 + 1):nroi, function(roi2) {
    high_corr(roi1, roi2, standard_fisherz, threshold = 0.8, std_err = std_err)
  }))
  sec_index <- which(p_vec <= 0.05)
  if (length(sec_index) > 0) {
    high_corr_mat[roi1, sec_index + roi1] <- TRUE
  }
}

## ------------------------------------------------------------
## 7. Prepare scaled data for VAR / Granger
## ------------------------------------------------------------

cat("Preparing scaled data for VAR/Granger...\n")

scaled_block <- roi_data[(t_corr + 1):T_total, , drop = FALSE]
scaled_block <- apply(scaled_block, 2, scale_column)
scaled_ts    <- ts(scaled_block)

T_eff <- nrow(scaled_ts)
cat("Effective time points for VAR:", T_eff, "\n")

## ------------------------------------------------------------
## 8. Select VAR lag order (capped)
## ------------------------------------------------------------

cat("Selecting VAR lag order...\n")
lag_selection <- VARselect(scaled_ts, lag.max = 10, type = "const")
order_p_raw   <- lag_selection$selection[1]
order_p       <- min(order_p_raw, 5)   # cap at 5 for stability

cat("Selected VAR lag order p =", order_p_raw, " (capped to ", order_p, ")\n")

## ------------------------------------------------------------
## 9. Compute EC (conditional Granger) matrix
## ------------------------------------------------------------

cat("Computing EC (conditional Granger) matrix...\n")

cg_mat <- matrix(NA_real_, nrow = nroi, ncol = nroi)

# Compute column-by-column (target j)
for (j in 1:nroi) {
  cg_mat[, j] <- unlist(mclapply(1:nroi, function(i) {
    intersection_conditional_granger(i, j, scaled_ts, order_p,
                                     high_corr_mat, max_cond = 20)
  }))
}

## ------------------------------------------------------------
## 10. Vectorize EC matrix and save (with subject_id + image_id)
## ------------------------------------------------------------

cg_vec <- as.vector(t(cg_mat))               # row-major flatten
diag_idx <- seq(1, length(cg_vec), by = nroi + 1)  # diagonal entries
cg_vec <- cg_vec[-diag_idx]                  # drop self-edges

# Optional: attach IDs as attributes inside the object as well
attr(cg_vec, "subject_id") <- subject_id
attr(cg_vec, "image_id")   <- image_id

# File name: e.g. "012_S_6073_11344002_ec.rds"
out_file <- file.path(ec_dir, paste0(subject_id, "_", image_id, "_ec.rds"))
saveRDS(cg_vec, file = out_file)

cat("Saved EC vector to:", out_file, "\n")
cat("  subject_id:", subject_id, "\n")
cat("  image_id  :", image_id, "\n")
cat("Done.\n")
