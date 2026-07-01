library(tidyverse)
library(data.table)
library(vars)
library(parallel)
library(kableExtra)
library(latex2exp)
library(gtsummary)
library(ComplexHeatmap)
library(circlize)
library(viridis)

# EC results directory: one RDS per subject (fMRI image)
ec_dir     <- file.path(base_dir, "fmri_mri_ec_denoised")

# Cohort info
cohort_file <- file.path(base_dir, "cohort.csv")

# AAL2 info (ROI names)
aal2_txt_file <- file.path(base_dir, "AAL2/aal2_ROI_names.txt")



source("functions.R")

aal2_txt <- read.table(aal2_txt_file, stringsAsFactors = FALSE)
nroi <- 120

## ------------------------------------------------------------
## 1. Load EC vectors and parse subject_id + image_id
## ------------------------------------------------------------

rds_files <- list.files(ec_dir, pattern = "_ec\\.rds$", full.names = TRUE)

ec_index <- tibble(
  file = rds_files,
  base = tools::file_path_sans_ext(basename(file))
) %>%
  # split into components: site, subj, visit, image_id, tag (= "ec")
  separate(base, into = c("part1","part2","part3","image_id","tag"),
           sep = "_", fill = "right", remove = FALSE) %>%
  mutate(
    subject_id = paste(part1, part2, part3, sep = "_")
  ) %>%
  select(file, subject_id, image_id)


# Read EC vectors in the same order as ec_index$file
X_list <- lapply(ec_index$file, function(f) {
  v <- readRDS(f)
  if (!is.numeric(v)) {
    stop("Non-numeric EC vector in file: ", f, " (class: ", paste(class(v), collapse=","), ")")
  }
  v
})

# Ensure all have the same length
p_vec <- sapply(X_list, length)
if (length(unique(p_vec)) != 1) {
  stop("Not all EC vectors have the same length! Lengths: ", paste(unique(p_vec), collapse = ", "))
}
p <- unique(p_vec)

# Combine into matrix: rows = subject-image pairs, cols = edges
X <- matrix(unlist(X_list), nrow = length(X_list), byrow = TRUE)
apply(X,2,sd) %>% summary()
hist(apply(X,2,sd))

# Use a composite ID "subject_id_image_id" as rownames
id_vec <- paste(ec_index$subject_id, ec_index$image_id, sep = "_")
rownames(X) <- id_vec

cat("EC matrix dimension: ", nrow(X), " subject-image pairs × ", ncol(X), " edges\n")

## ------------------------------------------------------------
## 2. Load cohort data and align with EC (by subject_id + image_id)
## ------------------------------------------------------------

cohort <- readr::read_csv(cohort_file, show_col_types = FALSE)
covars <- cohort %>%
  mutate(
    subj_img_id   = paste(subject_id, image_id, sep = "_"),
    gender        = sex,
    apoe          = apoe4,
    amyloid       = amyloid_status,
    research_grp  = entry_research_group
  ) %>%
  select(
    subj_img_id, subject_id, image_id,
    gender, age, education, apoe, amyloid, research_grp
  )

# Align by composite ID
common_ids <- intersect(rownames(X), covars$subj_img_id)
X      <- X[common_ids, , drop = FALSE]
covars <- covars %>% filter(subj_img_id %in% common_ids)

# Reorder covars to match X
covars <- covars[match(rownames(X), covars$subj_img_id), ]


## ------------------------------------------------------------
## 3. Diagnose NA in EC and drop very bad subjects
## ------------------------------------------------------------

na_total <- sum(is.na(X))
cat("Total NA entries in X:", na_total, "\n")

# Drop subjects with at least 1000 NAs in EC results
na_per_subject <- rowSums(is.na(X))
bad_threshold <- 1000
bad_ids <- names(which(na_per_subject > bad_threshold))
if (length(bad_ids) > 0) {
  X      <- X[!rownames(X) %in% bad_ids, , drop = FALSE]
  covars <- covars[!covars$subj_img_id %in% bad_ids, , drop = FALSE]
}


## ------------------------------------------------------------
## EC magnitude outlier check (row-wise)
## ------------------------------------------------------------

# Compute a global EC magnitude per subject-image pair
ec_scale <- apply(X, 1, function(v) {
  v_num <- v[!is.na(v)]
  if (!length(v_num)) return(NA_real_)
  sqrt(mean(v_num^2))
})
print(summary(ec_scale))

# threshold: median + k * MAD
k_out <- 3
med_ec <- median(ec_scale, na.rm = TRUE)
mad_ec <- mad(ec_scale, center = med_ec, constant = 1.4826)  # approx sd if normal
thr_ec <- med_ec + k_out * mad_ec

cat("Median EC RMS:", med_ec, "  MAD:", mad_ec, "  threshold (med +", k_out, "*MAD) =", thr_ec, "\n")

# Identify outliers
is_outlier <- ec_scale > thr_ec
outlier_idx <- which(is_outlier)

if (length(outlier_idx) > 0) {
  outlier_table <- tibble(
    subj_img_id = rownames(X)[outlier_idx],
    subject_id  = covars$subject_id[outlier_idx],
    image_id    = covars$image_id[outlier_idx],
    ec_rms      = ec_scale[outlier_idx]
  ) %>%
    arrange(desc(ec_rms))  
  # drop these from analysis
  X      <- X[!is_outlier, , drop = FALSE]
  covars <- covars[!is_outlier, , drop = FALSE]
  ec_scale <- ec_scale[!is_outlier]
  cat("After dropping EC outliers, remaining n =", nrow(X), "\n")
}

## ------------------------------------------------------------
## 4. Prepare covariates & summary table
## ------------------------------------------------------------

covars <- covars %>%
  mutate(
    gender = factor(
      gender,
      levels = c(1, 2, "Male", "Female"),
      labels = c("Male", "Female", "Male", "Female")
    ),
    apoe = factor(
      apoe,
      levels = c(0, 1),
      labels = c("W/O APOE4", "W/ APOE4")
    ),
    amyloid = factor(
      amyloid,
      levels = c(0, 1),
      labels = c("Amyloid Negative", "Amyloid Positive")
    ),
    research_grp = factor(
      research_grp,
      levels = c("CN", "SMC", "EMCI", "MCI")
    )
  )

labels1 <- list(
  amyloid      = "Amyloid PET Status",
  gender       = "Gender",
  age          = "Age (years)",
  education    = "Education (years)",
  apoe         = "Apolipoprotein E4 Carrier Status",
  research_grp = "Entry Research Group"
)

tbl1 <- covars %>%
  select(amyloid, gender, age, education, apoe, research_grp) %>%
  gtsummary::tbl_summary(
    by = amyloid,
    type = list(age ~ "continuous2", education ~ "continuous2"),
    statistic = list(
      gtsummary::all_continuous()  ~ c("{mean} ({sd})","{median} ({min}, {max})"),
      gtsummary::all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    digits = gtsummary::all_continuous() ~ 2,
    label  = labels1
  ) %>%
  gtsummary::add_overall() %>%
  gtsummary::add_p(gtsummary::all_categorical() ~ "chisq.test") %>%
  modify_fmt_fun(
    p.value ~ function(x) style_pvalue(x, digits = 3)
  ) %>%
  gtsummary::bold_labels()

tbl1

## ------------------------------------------------------------
## 5. IPW ATE estimation on EC outcomes
## ------------------------------------------------------------

n <- nrow(X)
p <- ncol(X)

D <- as.numeric(covars$amyloid) - 1   # 0/1 numeric
W <- model.matrix(~ gender + age + education + apoe, data = covars)
q <- ncol(W)

# Propensity scores
psm <- glm(D ~ W, family = binomial())
propensity_scores <- as.numeric(plogis(predict(psm)))

# IPW weights (ATE form)
weight <- ifelse(D == 1,  1 / propensity_scores,
                 -1 / (1 - propensity_scores))

# ATE
tau.estimate <- colMeans(X * weight, na.rm = TRUE)

# Fisher Information and Theta
Fisher.Information <- Reduce(`+`, lapply(1:n, function(i) {
  pi_i <- propensity_scores[i]
  pi_i * (1 - pi_i) * (W[i, ] %*% t(W[i, ]))
})) / n
Theta <- solve(Fisher.Information)

# H matrix
weight.H <- D * (1 - propensity_scores) / propensity_scores +
  (1 - D) * propensity_scores / (1 - propensity_scores)

H <- matrix(0, nrow = p, ncol = q)
for (j in 1:p) {
  H[j, ] <- colMeans(weight.H * X[, j] * W, na.rm = TRUE)
}

# Influence functions
eta <- matrix(0, nrow = n, ncol = p)
for (j in 1:p) {
  for (i in 1:n) {
    xij <- X[i, j]
    if (is.na(xij)) {
      eta[i, j] <- 0  # treat missing EC as no contribution
    } else {
      eta[i, j] <- xij * weight[i] -
        H[j, ] %*% Theta %*% W[i, ] * (D[i] - propensity_scores[i])
    }
  }
}

theta.var <- colMeans((eta - matrix(rep(tau.estimate, each = n), nrow = n))^2)
# all no less than 0.1

## ------------------------------------------------------------
## 6. Marginal inference: CI, p-values, masking
## ------------------------------------------------------------

CI_upper <- tau.estimate + 1.96 * sqrt(theta.var / n)
CI_lower <- tau.estimate - 1.96 * sqrt(theta.var / n)

z_stat  <- sqrt(n) * tau.estimate / sqrt(theta.var)
p_value <- 2 * (1 - pnorm(abs(z_stat)))

mat_pval <- vec_to_square_matrix(p_value, nroi)
mat_tau  <- vec_to_square_matrix(tau.estimate, nroi)

tau_masked <- tau.estimate
tau_masked[p_value >= 0.05] <- NA
mat_tau_masked <- vec_to_square_matrix(tau_masked, nroi)

## ------------------------------------------------------------
## 6. Heatmaps (all ATE, marginal significant ATE)
## ------------------------------------------------------------

acic_application_res_hm1 <- Heatmap(
  mat_tau,
  name  = "Estimated effect sizes",
  col   = viridis(100),
  na_col = "white",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  column_names_rot = 90,
  heatmap_legend_param = list(
    direction = "horizontal",
    labels_position = "right",
    title_gp  = gpar(fontsize = 22),
    labels_gp = gpar(fontsize = 12)
  )
)

acic_application_res_hm2 <- Heatmap(
  mat_tau_masked,
  column_title = "Significant ATE (p<0.05)",
  name  = "ATE",
  col   = viridis(100),
  na_col = "white",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 3.5),
  column_names_rot = 90,
  heatmap_legend_param = list(
    direction = "horizontal",
    labels_position = "right",
    title_gp  = gpar(fontsize = 6),
    labels_gp = gpar(fontsize = 6)
  )
)
pdf(
  file = "application_res_hm1.pdf",
  width  = 19,   
  height = 19    
)

draw(acic_application_res_hm1,
     column_title = "Estimated average causal effects on each pair of ROIs",
     column_title_gp = gpar(fontsize = 30, fontface = "bold"),
     heatmap_legend_side = "bottom")
dev.off()

## ------------------------------------------------------------
## 7. Step-down + augmentation multiple testing
## ------------------------------------------------------------

Tstat <- sqrt(n) * abs(tau.estimate) / sqrt(theta.var)

set.seed(2025)
B   <- 1000
aph <- 0.05
c   <- 0.1

z <- matrix(0, nrow = p, ncol = B)
for (b in 1:B) {
  g <- rnorm(n)
  temp <- colSums((eta - matrix(rep(tau.estimate, each = n), nrow = n)) * g)
  z[, b] <- temp / sqrt(n) / sqrt(theta.var)
}

Tstat0      <- Tstat
result      <- rep(0, p)
z.initial   <- z
Tstat.initial <- Tstat0

repeat {
  Tstat.max <- max(abs(Tstat.initial))
  if (Tstat.max == 0) break
  index.temp <- which(abs(Tstat0) == Tstat.max)
  z.max <- apply(abs(z.initial), 2, max)
  z.max.quan <- quantile(z.max, 1 - aph)
  if (Tstat.max < z.max.quan) break
  Tstat.initial[index.temp] <- 0
  z.initial[index.temp, ]   <- 0
  result[index.temp]        <- 1
}

# Augmentation
size    <- sum(result)
num_add <- floor(c * size / (1 - c))
if (num_add >= 1) {
  test_replace <- sort(Tstat.initial, decreasing = TRUE)
  for (i in 1:num_add) {
    index.temp <- which(abs(Tstat0) == test_replace[i])
    result[index.temp] <- 1
  }
}

tau_masked_mt <- tau.estimate
tau_masked_mt[result == 0] <- NA
mat_tau_masked_mt <- vec_to_square_matrix(tau_masked_mt, nroi)

which(!is.na(mat_tau_masked_mt), arr.ind = T)
aal2_txt[35,]
aal2_txt[55,]
# Remove rows/cols that are all NA
mat_tau_masked_mt_cleaned <- mat_tau_masked_mt[
  !apply(is.na(mat_tau_masked_mt), 1, all),
  !apply(is.na(mat_tau_masked_mt), 2, all)
]

acic_application_res_hm3 <- Heatmap(
  mat_tau_masked_mt_cleaned,
  column_title = "Significant ATE (after multiple testing)",
  name = "ATE",
  col  = viridis(100),
  na_col = "white",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8),
  column_names_rot = 90,
  heatmap_legend_param = list(
    direction = "horizontal",
    labels_position = "bottom",
    title_gp  = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8)
  )
)

draw(acic_application_res_hm3,
     column_title = "Result Map (Multiple Testing)",
     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
     heatmap_legend_side = "bottom")

# confidence interval
mat_tau[35,55]
q  <- unname(quantile(z.max, 1 - aph)) # as we only have one component remains significant
CI_upper_simultaneous <- tau.estimate + q * sqrt(theta.var / n)
CI_lower_simultaneous <- tau.estimate - q * sqrt(theta.var / n)
CI_upper_35_55 <- vec_to_square_matrix(CI_upper_simultaneous, nroi)[35,55]
CI_lower_35_55 <- vec_to_square_matrix(CI_lower_simultaneous, nroi)[35,55]

## ------------------------------------------------------------
## 8. Save results
## ------------------------------------------------------------

tau_vec      <- tau.estimate
sig_vec      <- as.integer(result)
ATE_mat      <- vec_to_square_matrix(tau_vec, nroi)
SIG_mat      <- matrix(sig_vec, nroi, nroi, byrow = FALSE)

write.csv(tau_vec, "tau_estimate_vec.csv", row.names = FALSE)
write.csv(sig_vec, "sig_mask_vec.csv",    row.names = FALSE)
write.csv(ATE_mat, "ATE_mat.csv",         row.names = FALSE)
write.csv(SIG_mat, "SIG_mat.csv",         row.names = FALSE)


## ============================================================
## OVERLAP DIAGNOSTICS, WEIGHT DIAGNOSTICS, SENSITIVITY ANALYSIS
## ============================================================
idx_mat  <- vec_to_square_matrix(seq_along(tau.estimate), nroi)
lin_3555 <- idx_mat[35, 55]

## ------------------------------------------------------------
## 1. Overlap diagnostics
## ------------------------------------------------------------
ps_df <- tibble(ps = propensity_scores,
                grp = factor(D, levels = c(0, 1),
                             labels = c("Amyloid Negative", "Amyloid Positive")))

ps_summary <- ps_df %>%
  group_by(grp) %>%
  summarise(min = min(ps), q25 = quantile(ps, .25),
            median = median(ps), q75 = quantile(ps, .75),
            max = max(ps), .groups = "drop")
cat("\nPropensity score summary by group:\n"); print(ps_summary)

overlap_lo <- max(tapply(propensity_scores, D, min))
overlap_hi <- min(tapply(propensity_scores, D, max))
n_outside  <- sum(propensity_scores < overlap_lo | propensity_scores > overlap_hi)
cat("\nOverall PS range: [", round(min(propensity_scores), 3), ",",
    round(max(propensity_scores), 3), "]\n")
cat("Common support:  [", round(overlap_lo, 3), ",", round(overlap_hi, 3), "]\n")
cat("Subjects outside common support:", n_outside, "of", n, "\n")

# figure
p_overlap <- ggplot(ps_df, aes(x = ps, fill = grp)) +
  geom_histogram(data = subset(ps_df, grp == "Amyloid Positive"),
                 aes(y = after_stat(count)), bins = 20, alpha = .6) +
  geom_histogram(data = subset(ps_df, grp == "Amyloid Negative"),
                 aes(y = -after_stat(count)), bins = 20, alpha = .6) +
  geom_hline(yintercept = 0, linewidth = .3) +
  labs(x = "Estimated propensity score", y = "Count", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")
ggsave("ps_overlap.pdf", p_overlap, width = 6, height = 4)

## ------------------------------------------------------------
## 2. Weight diagnostics
## ------------------------------------------------------------
cat("\nMax IPW weight (abs):", round(max(abs(weight)), 2), "\n")
cat("Weight quantiles (|w|):\n")
print(round(quantile(abs(weight), c(.5, .9, .99, 1)), 2))

## ------------------------------------------------------------
## 3. Trimming sensitivity for the surviving edge (35 -> 55)
##    Re-estimate tau AND its z-statistic under symmetric trimming,
##    so we can report whether it stays significant, not just its size.
## ------------------------------------------------------------
trim_grid <- c(0, 0.01, 0.025, 0.05, 0.10)

trim_res <- lapply(trim_grid, function(t) {
  keep <- propensity_scores > t & propensity_scores < (1 - t)
  nt   <- sum(keep)
  Dt   <- D[keep]
  Wt   <- W[keep, , drop = FALSE]
  Xt   <- X[keep, , drop = FALSE]
  
  # refit propensity model on the trimmed sample for internal consistency
  psm_t <- glm(Dt ~ Wt, family = binomial())
  pst   <- as.numeric(plogis(predict(psm_t)))
  wt    <- ifelse(Dt == 1, 1 / pst, -1 / (1 - pst))
  
  tau_t <- colMeans(Xt * wt, na.rm = TRUE)
  
  # variance for edge (35,55) via its influence function on the trimmed sample
  q_t     <- ncol(Wt)
  FI_t    <- Reduce(`+`, lapply(seq_len(nt), function(i)
    pst[i] * (1 - pst[i]) * (Wt[i, ] %*% t(Wt[i, ]))))/nt
  Theta_t <- solve(FI_t)
  wH_t    <- Dt * (1 - pst)/pst + (1 - Dt) * pst/(1 - pst)
  Hj      <- colMeans(wH_t * Xt[, lin_3555] * Wt, na.rm = TRUE)
  xj      <- Xt[, lin_3555]
  eta_j   <- ifelse(is.na(xj), 0,
                    xj * wt - as.numeric(Hj %*% Theta_t %*% t(Wt)) *
                      (Dt - pst))
  var_j   <- mean((eta_j - tau_t[lin_3555])^2)
  z_j     <- sqrt(nt) * tau_t[lin_3555] / sqrt(var_j)
  
  tibble(trim = t, n_kept = nt,
         tau_3555 = tau_t[lin_3555],
         z_3555   = z_j,
         p_3555   = 2 * (1 - pnorm(abs(z_j))))
})
cat("\nTrimming sensitivity for edge (35,55):\n")
print(bind_rows(trim_res))

## ------------------------------------------------------------
## 4. E-value for the discovery
##    Sign fix: feed the absolute standardized effect; use the
##    simultaneous CI limit nearest the null.
## ------------------------------------------------------------
est_3555   <- tau.estimate[lin_3555]
ci_lo_3555 <- CI_lower_35_55        # simultaneous, from your script
ci_hi_3555 <- CI_upper_35_55
sd_3555    <- sd(X[, lin_3555], na.rm = TRUE)

library(EValue)
d_est <- abs(est_3555) / sd_3555
se_d  <- abs(ci_hi_3555 - ci_lo_3555) / (2 * 1.96) / sd_3555
ev    <- evalues.MD(est = d_est, se = se_d)

cat("\nEdge (35,55): est =", round(est_3555, 3),
    " simult. CI = [", round(ci_lo_3555, 3), ",", round(ci_hi_3555, 3), "]\n")
cat("Edge SD =", round(sd_3555, 3),
    " standardized d =", round(d_est, 3), "\n")
cat("E-values (point estimate and CI limit nearest null):\n")
print(ev)
