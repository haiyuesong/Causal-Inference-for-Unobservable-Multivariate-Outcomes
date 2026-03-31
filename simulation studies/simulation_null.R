suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppArmadillo)
  library(MASS)
  library(parallel)
})

## =========================================================
##  C++ kernels
## =========================================================
Rcpp::sourceCpp("/home/ec_kernels.cpp")

## =========================================================
##  Helpers: baseline, masks, generator (same as alt)
## =========================================================

make_stable_transition_matrix <- function(B, buffer = 0.05) {
  ev  <- eigen(B, only.values = TRUE)$values
  rho <- max(Mod(ev))
  if (rho >= 1) B <- B / (rho + buffer)
  B
}

build_baseline_B0 <- function(p,
                              block_prob_within = 0.35,
                              diag_range = c(0.72, 0.85),
                              weight_range_within = c(0.03, 0.08),
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  block_p <- p/3; stopifnot(block_p * 3 == p)
  diag_base <- runif(p, diag_range[1], diag_range[2])
  
  w1 <- matrix(runif(block_p^2, weight_range_within[1], weight_range_within[2]), block_p, block_p)
  w2 <- matrix(runif(block_p^2, weight_range_within[1], weight_range_within[2]), block_p, block_p)
  w3 <- matrix(runif(block_p^2, weight_range_within[1], weight_range_within[2]), block_p, block_p)
  
  m1 <- matrix(rbinom(block_p^2, 1, block_prob_within), block_p, block_p)
  m2 <- matrix(rbinom(block_p^2, 1, block_prob_within), block_p, block_p)
  m3 <- matrix(rbinom(block_p^2, 1, block_prob_within), block_p, block_p)
  
  Z <- matrix(0, block_p, block_p)
  W_off <- rbind(
    cbind(w1*m1, Z,       Z),
    cbind(Z,      w2*m2,  Z),
    cbind(Z,      Z,      w3*m3)
  )
  diag(W_off) <- 0
  
  B0 <- diag(diag_base, p, p) + W_off
  make_stable_transition_matrix(B0)
}

.block_ids <- function(p) rep(1:3, each = p/3)

.within_deg <- function(B0) {
  p   <- nrow(B0)
  blk <- .block_ids(p)
  W   <- (abs(B0) > 0) * 1
  deg_in  <- deg_out <- integer(p)
  for (b in 1:3) {
    idx <- which(blk == b)
    Wb  <- W[idx, idx, drop = FALSE]
    deg_in[idx]  <- colSums(Wb)
    deg_out[idx] <- rowSums(Wb)
  }
  list(deg_in = deg_in, deg_out = deg_out, blk = blk)
}

build_treat_mask_lowdeg <- function(B0, k_per_pair = 3, seed = 7L) {
  set.seed(seed)
  p   <- nrow(B0); block_p <- p/3
  info <- .within_deg(B0); blk <- info$blk
  Z   <- matrix(0, block_p, block_p)
  
  pick_no_share <- function(src_b, tgt_b, k) {
    src <- which(blk == src_b)
    tgt <- which(blk == tgt_b)
    used_s <- used_t <- integer(0)
    M <- matrix(0, block_p, block_p)
    cand <- expand.grid(s = src, t = tgt)
    cand$score <- info$deg_out[cand$s] + info$deg_in[cand$t]
    cand <- cand[order(cand$score), ]
    for (r in seq_len(nrow(cand))) {
      if (sum(M) >= k) break
      s <- cand$s[r]; t <- cand$t[r]
      if (s %in% used_s || t %in% used_t) next
      M[(t - (tgt_b-1)*block_p), (s - (src_b-1)*block_p)] <- 1
      used_s <- c(used_s, s); used_t <- c(used_t, t)
    }
    M
  }
  
  B12 <- pick_no_share(1, 2, k_per_pair)
  B23 <- pick_no_share(2, 3, k_per_pair)
  B13 <- pick_no_share(1, 3, k_per_pair)
  
  # j->i orientation
  rbind(
    cbind(Z, B12, Z),
    cbind(Z, Z,  B23),
    cbind(B13, Z, Z)
  )
}

shrink_ring <- function(B0, Streat, gamma_end = 0.20, gamma_1hop = 0.40) {
  p   <- nrow(B0); blk <- .block_ids(p)
  B   <- B0
  E   <- union(which(rowSums(Streat)!=0), which(colSums(Streat)!=0))
  if (length(E)) {
    for (v in E) {
      idx <- which(blk == blk[v])
      B[v, idx] <- gamma_end  * B[v, idx]
      B[idx, v] <- gamma_end  * B[idx, v]
    }
    N1 <- integer(0)
    for (v in E) {
      idx <- which(blk == blk[v])
      nb  <- union(which(abs(B0[v, idx])>0), which(abs(B0[idx, v])>0))
      if (length(nb)) N1 <- union(N1, idx[nb])
    }
    for (u in setdiff(N1, E)) {
      idx <- which(blk == blk[u])
      B[u, idx] <- gamma_1hop * B[u, idx]
      B[idx, u] <- gamma_1hop * B[idx, u]
    }
  }
  make_stable_transition_matrix(B)
}

build_conf_mask_disjoint <- function(p, prob = 0.30, Streat, seed = 11L) {
  set.seed(seed)
  block_p <- p/3; Z <- matrix(0, block_p, block_p)
  blk <- matrix(rbinom(block_p*block_p, 1, prob), block_p, block_p)
  S <- rbind(
    cbind(Z, blk, Z),
    cbind(blk, Z, blk),
    cbind(Z, blk, Z)
  )
  S[Streat == 1] <- 0
  S
}

VAR_generation_1 <- function(A_1, time_length = 2000, drop_length = 500,
                             lag = 1, num_var = 3, sigma = 1) {
  X <- matrix(runif(num_var*lag, -0.5, 0.5), nrow=num_var, ncol=lag)
  for (t in (lag+1):time_length) {
    X <- cbind(X, A_1 %*% X[, t-1] + matrix(rnorm(num_var, sd=sigma), nrow=num_var))
  }
  X[, -(1:drop_length), drop = FALSE]
}

vec_row_offdiag <- function(M) {
  p <- nrow(M)
  v <- as.vector(t(M))
  v[-seq(1, length(v), by = p + 1)]
}

## =========================================================
##  Step-down + augmentation
## =========================================================

stepdown_select <- function(Tstat, eta_center, var_j,
                            B = 1000, alpha = 0.05, c_aug = 0.10) {
  n <- nrow(eta_center)
  G <- matrix(rnorm(n*B), n, B)
  Z <- t(eta_center) %*% G / sqrt(n)
  Z <- Z / sqrt(var_j)
  
  T0 <- Tstat
  res  <- integer(length(T0))
  Tcur <- T0
  Zcur <- Z
  
  repeat {
    Tm <- max(abs(Tcur))
    if (Tm == 0) break
    idx <- which(abs(T0) == Tm)
    zmax <- apply(abs(Zcur), 2, max)
    zq   <- as.numeric(quantile(zmax, 1-alpha, names=FALSE))
    if (Tm < zq) break
    Tcur[idx] <- 0
    Zcur[idx,] <- 0
    res[idx] <- 1L
  }
  
  k <- sum(res)
  add <- floor(c_aug * k / (1 - c_aug))
  if (add >= 1) {
    ord  <- order(abs(Tcur), decreasing = TRUE)
    take <- head(ord[abs(Tcur[ord]) > 0], add)
    res[take] <- 1L
  }
  res
}

ipw_stat <- function(results_mat, D, ps) {
  eps <- 0.10
  ps  <- pmin(pmax(ps, eps), 1 - eps)
  w   <- ifelse(D == 1, 1/ps, -1/(1-ps))
  G   <- sweep(results_mat, 1, w, `*`)
  tau <- colMeans(G)
  Gc  <- sweep(G, 2, tau, "-")
  var <- colMeans(Gc^2) + 1e-12
  list(tau = tau, Gc = Gc, var = var)
}

## =========================================================
##  IPW τ-hat and bias summary (under null)
## =========================================================

ipw_tau <- function(results_mat, D, ps) {
  eps <- 0.10
  ps  <- pmin(pmax(ps, eps), 1 - eps)
  w   <- ifelse(D == 1, 1/ps, -1/(1-ps))
  colMeans(sweep(results_mat, 1, w, `*`))
}

bias_summ <- function(tau_vec) {
  c(
    mean_bias = mean(tau_vec),              # signed mean
    mean_abs  = mean(abs(tau_vec)),
    l2_mean   = sqrt(mean(tau_vec^2)),
    max_abs   = max(abs(tau_vec))
  )
}

## =========================================================
##  One replicate under global null
## =========================================================

simulate_once_null_stepdown_bias <- function(p, n,
                                             Time_length = 500,
                                             cap_k = 12,
                                             lambda_w = 0.10,
                                             seed = 2025) {
  set.seed(seed)
  
  q <- 5
  W <- scale(MASS::mvrnorm(n, mu = rep(0,q), Sigma = diag(q)))
  beta_true <- c(1.2, -1.0, 0.8, 0.6, -0.5)
  logits <- -0.3 + as.vector(W %*% beta_true)
  pr     <- plogis(logits)
  D      <- rbinom(n, 1, pr)
  
  ps_true  <- plogis(predict(glm(D ~ W, family = binomial())))
  ps_wrong <- plogis(predict(glm(D ~ 1, family = binomial())))
  
  # baseline & masks (no treatment term under null)
  if (p == 51) {
    B0 <- build_baseline_B0(p, block_prob_within=0.30, weight_range_within=c(0.025,0.06), seed=seed)
    cap_here   <- max(cap_k, 18)
    k_per_pair <- 4
  } else {
    B0 <- build_baseline_B0(p, block_prob_within=0.35, weight_range_within=c(0.03,0.08), seed=seed)
    cap_here   <- cap_k
    k_per_pair <- 2
  }
  
  Streat <- build_treat_mask_lowdeg(B0, k_per_pair = k_per_pair, seed = seed+7)
  Sconf  <- build_conf_mask_disjoint(p, prob = 0.30, Streat = Streat, seed = seed+11)
  B0     <- shrink_ring(B0, Streat, gamma_end = 0.20, gamma_1hop = 0.40)
  
  a_i <- plogis(as.vector(W %*% beta_true))
  
  get_vecs <- function(i) {
    Bi <- make_stable_transition_matrix(B0 + lambda_w * a_i[i] * Sconf)
    sig <- t(VAR_generation_1(
      A_1 = Bi,
      time_length = Time_length*2 + 500,
      drop_length = 500,
      lag = 1, num_var = p, sigma = 1
    ))
    X1 <- sig[1:Time_length, ]
    X2 <- sig[(Time_length+1):(2*Time_length), ]
    
    sds <- apply(X1, 2, sd)
    sds[!is.finite(sds) | sds < 1e-8] <- 1
    X1 <- sweep(X1, 2, sds, "/")
    X2 <- sweep(X2, 2, sds, "/")
    
    cond_prop  <- condset_topk_cpp(X1, r_thresh=0.1, alpha=0.05, cap_k=cap_here)
    cond_leak  <- condset_topk_cpp(X2, r_thresh=0.1, alpha=0.05, cap_k=cap_here)
    cond_empty <- replicate(p, integer(0), simplify = FALSE)
    
    F_prop  <- batched_F_cpp(X2, cond_prop)
    F_leak  <- batched_F_cpp(X2, cond_leak)
    F_empty <- batched_F_cpp(X2, cond_empty)
    
    list(
      v_prop  = vec_row_offdiag(F_prop),
      v_leak  = vec_row_offdiag(F_leak),
      v_empty = vec_row_offdiag(F_empty)
    )
  }
  
  L <- lapply(1:n, get_vecs)
  res_prop  <- do.call(rbind, lapply(L, `[[`, "v_prop"))
  res_leak  <- do.call(rbind, lapply(L, `[[`, "v_leak"))
  res_empty <- do.call(rbind, lapply(L, `[[`, "v_empty"))
  
  ## FWER via step-down (same as alt)
  stat_prop  <- ipw_stat(res_prop,  D, ps_true)
  stat_wps   <- ipw_stat(res_prop,  D, ps_wrong)
  stat_leak  <- ipw_stat(res_leak,  D, ps_true)
  stat_both  <- ipw_stat(res_empty, D, ps_wrong)
  
  sel_prop <- stepdown_select(
    sqrt(n)*abs(stat_prop$tau)/sqrt(stat_prop$var),
    stat_prop$Gc, stat_prop$var
  )
  sel_wps  <- stepdown_select(
    sqrt(n)*abs(stat_wps$tau)/sqrt(stat_wps$var),
    stat_wps$Gc, stat_wps$var
  )
  sel_leak <- stepdown_select(
    sqrt(n)*abs(stat_leak$tau)/sqrt(stat_leak$var),
    stat_leak$Gc, stat_leak$var
  )
  sel_both <- stepdown_select(
    sqrt(n)*abs(stat_both$tau)/sqrt(stat_both$var),
    stat_both$Gc, stat_both$var
  )
  
  f_prop <- as.integer(any(sel_prop == 1L))
  f_wps  <- as.integer(any(sel_wps  == 1L))
  f_leak <- as.integer(any(sel_leak == 1L))
  f_both <- as.integer(any(sel_both == 1L))
