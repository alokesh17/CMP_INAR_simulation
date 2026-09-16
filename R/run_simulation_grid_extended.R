############################################################################
# Extended simulation grid for the ZICMP-INAR(1) model: adds one MORE
# strongly overdispersed regime and one MORE strongly underdispersed regime
# to the three already in the draft (Table~\ref{tab:sim_placeholder}), per
# the request to probe the dispersion spectrum further out from nu=1 than
# the original nu in {0.5, 1.0, 1.5}.
#
# This is the SAME driver as run_simulation_grid.R (same simulator,
# comparators, Stan file, and output format) with the regime list widened
# to
#     NUS <- c(0.2, 0.5, 1.0, 1.5, 2.5)
# so the two new cells slot into the existing pipeline and the existing
# emit_latex()/plotting code (see make_sim_figures.R) without change --
# only the regime *labels* need to know about the two new rows, which is
# handled generically below (labels are generated from NUS, not
# hardcoded).
#
# Rationale for the two added values:
#   nu = 0.2  -- markedly stronger overdispersion than nu=0.5. The CMP mean
#                approximation and hybrid log-Z switch are both exercised
#                harder here (lambda^{-1/nu} shrinks faster), so this also
#                doubles as a numerical stress test of Section 2.7's hybrid
#                approximation away from its comfortable middle range.
#   nu = 2.5  -- markedly stronger underdispersion than nu=1.5, deep enough
#                into underdispersion that the CMP distribution is
#                approaching its near-deterministic (large-nu) boundary;
#                worth checking that Stan still mixes cleanly and that
#                lambda calibration (solve_lambda_for_mean) still converges
#                there.
#
# Nothing else about the design changes: alpha_true = rho_true = 0.30,
# n in {100, 200, 400, 600}, R replications per cell, same priors, same
# convergence screen (Rhat < 1.05).
#
# HOW TO USE: identical to run_simulation_grid.R -- put this file next to
# ZIINAR1-CMP-fast-reparam.stan (or switch STAN_FILE below) and run with
# RUN_MODE <- "quick" first, then "pilot", then "full". This script writes
# its own checkpoint/output files (suffixed _ext) so it will not clobber
# results from the original 3-regime run if both are run in the same
# directory.
############################################################################

library(rstan)
library(COMPoissonReg)   # for rcmp(), same as the original simulation code
rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))

set.seed(2026)

# ---------------------------------------------------------------------
# 0. Config -- adjust these for your machine / time budget
# ---------------------------------------------------------------------
RUN_MODE <- "pilot"    # <- "quick" | "pilot" | "full"

STAN_FILE   <- "ZIINAR1-CMP-fast-reparam.stan"
M           <- 300
HYBRID_TOL  <- 1e-6
FF          <- 0
ALPHA_TRUE  <- 0.30
RHO_TRUE    <- 0.30
TARGET_MU   <- 3.0

# The two NEW regimes are nu = 0.2 (stronger overdispersion) and nu = 2.5
# (stronger underdispersion). Set EXTRA_ONLY <- TRUE to run just the two
# new cells (cheaper if you already have the original 3-regime results and
# only want to add the extremes), or FALSE to redo the full 5-regime grid
# in one pass.
EXTRA_ONLY <- FALSE

if (RUN_MODE == "quick") {
  NUS <- c(2.5); NS <- c(100); R <- 3
  ITER <- 800; WARMUP <- 400; CHAINS <- 2
} else if (RUN_MODE == "pilot") {
  NUS <- c(0.2); NS <- c(600); R <- 15
  ITER <- 2000; WARMUP <- 1000; CHAINS <- 2
} else {
  NUS <- if (EXTRA_ONLY) c(0.2, 2.5) else c(0.2, 0.5, 1.0, 1.5, 2.5)
  NS  <- c(100, 200, 400, 600)
  R   <- 30
  ITER <- 2000; WARMUP <- 1000; CHAINS <- 2
}
QUICK_TEST <- (RUN_MODE != "full")

# ---------------------------------------------------------------------
# 1. Simulator (unchanged from run_simulation_grid.R)
# ---------------------------------------------------------------------
thin_operator <- function(x, alpha) sum(rbinom(x, size = 1, prob = alpha))

simul_zinarCMP <- function(n, alpha, lambda, nu, zi_prob) {
  zinp_inar <- integer(n)
  zinp_inar[1] <- rcmp(1, lambda, nu)
  for (t in 2:n) {
    thinned <- thin_operator(zinp_inar[t - 1], alpha)
    if (runif(1) < zi_prob) {
      innovation <- 0
    } else {
      innovation <- rcmp(1, lambda, nu)
    }
    zinp_inar[t] <- thinned + innovation
  }
  zinp_inar
}

# ---------------------------------------------------------------------
# 2. Exact CMP moments (unchanged)
# ---------------------------------------------------------------------
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

cmp_moments <- function(lambda, nu, kmax = 600) {
  # kmax raised from 400 to 600 relative to the original script: at
  # nu = 0.2 the CMP tail is heavier, so more terms are needed for the
  # truncated-series moment sum to converge to the same tolerance.
  k <- 0:kmax
  logw <- k * log(lambda) - nu * lgamma(k + 1)
  logZ <- logsumexp(logw)
  p <- exp(logw - logZ)
  mu <- sum(k * p)
  c(mean = mu, var = sum(k^2 * p) - mu^2)
}

solve_lambda_for_mean <- function(nu, target_mu, kmax = 600) {
  f <- function(lam) cmp_moments(lam, nu, kmax)["mean"] - target_mu
  uniroot(f, c(1e-6, 1e6), tol = 1e-8)$root
}

# ---------------------------------------------------------------------
# 3. Comparators (unchanged)
# ---------------------------------------------------------------------
cls_alpha <- function(y) {
  y1 <- y[-length(y)]; y2 <- y[-1]
  m1 <- mean(y1); m2 <- mean(y2)
  sum((y1 - m1) * (y2 - m2)) / sum((y1 - m1)^2)
}

mom_oracle_lambda_nu <- function(y, alpha_true, rho_true, kmax = 600) {
  mu_hat <- mean(y); s2_hat <- var(y)
  muU_target <- mu_hat * (1 - alpha_true) / (1 - rho_true)
  rhs <- (1 - alpha_true^2) * s2_hat - alpha_true * (1 - alpha_true) * mu_hat
  sigU2_target <- (rhs - rho_true * (1 - rho_true) * muU_target^2) / (1 - rho_true)
  if (!is.finite(sigU2_target) || sigU2_target <= 0 || muU_target <= 0) {
    return(c(lambda = NA_real_, nu = NA_real_))
  }
  resid <- function(nu) {
    lam <- tryCatch(solve_lambda_for_mean(nu, muU_target, kmax), error = function(e) NA)
    if (is.na(lam)) return(NA_real_)
    cmp_moments(lam, nu, kmax)["var"] - sigU2_target
  }
  nu_hat <- tryCatch({
    # search range widened to (0.02, 12) to cover the two new regimes
    r_lo <- resid(0.02); r_hi <- resid(12)
    if (!is.finite(r_lo) || !is.finite(r_hi) || sign(r_lo) == sign(r_hi)) return(NA_real_)
    uniroot(resid, c(0.02, 12), tol = 1e-6)$root
  }, error = function(e) NA_real_)
  if (is.na(nu_hat)) return(c(lambda = NA_real_, nu = NA_real_))
  lam_hat <- tryCatch(solve_lambda_for_mean(nu_hat, muU_target, kmax), error = function(e) NA_real_)
  c(lambda = lam_hat, nu = nu_hat)
}

# ---------------------------------------------------------------------
# 4. Calibrate lambda per regime
# ---------------------------------------------------------------------
lambdas <- setNames(sapply(NUS, solve_lambda_for_mean, target_mu = TARGET_MU), as.character(NUS))
cat("Calibrated lambdas (target CMP mean =", TARGET_MU, "):\n"); print(lambdas)
for (nu in NUS) {
  mv <- cmp_moments(lambdas[[as.character(nu)]], nu)
  cat(sprintf("  nu=%.2f  lambda=%.4f  mean=%.3f  var=%.3f  d_U=%.3f\n",
              nu, lambdas[[as.character(nu)]], mv["mean"], mv["var"], mv["var"] / mv["mean"]))
}

# ---------------------------------------------------------------------
# 5. Compile once, reuse for every dataset
# ---------------------------------------------------------------------
cat("\nCompiling", STAN_FILE, "...\n")
mod <- stan_model(STAN_FILE)

# ---------------------------------------------------------------------
# 6. Main grid
# ---------------------------------------------------------------------
results <- list()
t_start <- Sys.time()

for (nu_true in NUS) {
  lam_true <- lambdas[[as.character(nu_true)]]
  for (n in NS) {
    key <- paste0("nu=", nu_true, "_n=", n)
    cat("\n===", key, "===\n")
    rows <- vector("list", R)
    for (r in seq_len(R)) {
      y <- simul_zinarCMP(n, ALPHA_TRUE, lam_true, nu_true, RHO_TRUE)

      fit <- tryCatch(
        sampling(mod,
                 data = list(T = n, y = y, M = M, hybrid_tol = HYBRID_TOL, ff = FF),
                 chains = CHAINS, iter = ITER, warmup = WARMUP,
                 seed = 1000 + r, refresh = 0,
                 control = list(adapt_delta = 0.95)),  # raised from 0.9: the two
                 # extra regimes sit closer to the nu boundary of the CMP family
                 # (near-Bernoulli at nu=0.2's tail, near-deterministic as nu
                 # grows), so a slightly more conservative step size helps
                 # avoid divergences there.
        error = function(e) { message("  sampling failed: ", conditionMessage(e)); NULL }
      )
      if (is.null(fit)) next

      s <- summary(fit, pars = c("alpha", "lambda", "nu", "rho"))$summary
      rhat_ok <- all(s[, "Rhat"] < 1.05, na.rm = TRUE)

      post <- rstan::extract(fit, pars = c("lambda", "nu"))
      lamnu_cor <- suppressWarnings(cor(post$lambda, post$nu))

      a_cls <- cls_alpha(y)
      mom <- mom_oracle_lambda_nu(y, ALPHA_TRUE, RHO_TRUE)

      rows[[r]] <- data.frame(
        rep = r,
        alpha_mean = s["alpha", "mean"],   alpha_sd = s["alpha", "sd"],
        lambda_mean = s["lambda", "mean"], lambda_sd = s["lambda", "sd"],
        nu_mean = s["nu", "mean"],         nu_sd = s["nu", "sd"],
        rho_mean = s["rho", "mean"],       rho_sd = s["rho", "sd"],
        rhat_ok = rhat_ok,
        lambda_nu_cor = lamnu_cor,
        alpha_cls = a_cls,
        lambda_mom = unname(mom["lambda"]), nu_mom = unname(mom["nu"])
      )
      cat(sprintf("  rep %d/%d  Rhat_ok=%s  alpha=%.3f  lambda=%.3f  nu=%.3f  rho=%.3f  cor(lam,nu)=%.3f  (cls_alpha=%.3f, mom_nu=%.3f)\n",
                   r, R, rhat_ok, s["alpha", "mean"], s["lambda", "mean"],
                   s["nu", "mean"], s["rho", "mean"], lamnu_cor, a_cls, unname(mom["nu"])))
    }
    results[[key]] <- do.call(rbind, rows[!sapply(rows, is.null)])
    saveRDS(results, "sim_results_partial_ext.rds")
    cat(sprintf("  elapsed so far: %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }
}

saveRDS(results, "sim_results_final_ext.rds")
cat("\nTOTAL TIME (min):", as.numeric(difftime(Sys.time(), t_start, units = "mins")), "\n")

# ---------------------------------------------------------------------
# 7. Summarize (unchanged logic)
# ---------------------------------------------------------------------
summarize_cell <- function(df) {
  data.frame(
    n_reps = nrow(df),
    alpha_MCmean = mean(df$alpha_mean),  alpha_MCsd = sd(df$alpha_mean),
    lambda_MCmean = mean(df$lambda_mean), lambda_MCsd = sd(df$lambda_mean),
    nu_MCmean = mean(df$nu_mean),        nu_MCsd = sd(df$nu_mean),
    rho_MCmean = mean(df$rho_mean),      rho_MCsd = sd(df$rho_mean),
    alpha_cls_mean = mean(df$alpha_cls, na.rm = TRUE),
    nu_mom_mean = mean(df$nu_mom, na.rm = TRUE),
    rhat_ok_frac = mean(df$rhat_ok),
    lambda_nu_cor_mean = mean(df$lambda_nu_cor, na.rm = TRUE)
  )
}

summary_tab <- do.call(rbind, lapply(names(results), function(k) {
  cbind(cell = k, summarize_cell(results[[k]]))
}))
write.csv(summary_tab, "sim_summary_table_ext.csv", row.names = FALSE)
cat("\n--- Summary table (also written to sim_summary_table_ext.csv) ---\n")
print(summary_tab)

# ---------------------------------------------------------------------
# 8. Emit LaTeX rows for an expanded Table 1 (or a standalone appendix
#    table if you'd rather keep the original 3-regime Table 1 untouched
#    and add the two new regimes as a small supplementary table -- either
#    way, this generates rows for EVERY regime in NUS, generically, not
#    just the original three.
# ---------------------------------------------------------------------
fmt <- function(m, s) sprintf("%.3f (%.3f)", m, s)

nu_regime_label <- function(nu) {
  if (nu < 0.5) return(sprintf("Strong overdispersion regime ($\\nu=%.1f$)", nu))
  if (nu < 1.0) return(sprintf("Overdispersion regime ($\\nu=%.1f$)", nu))
  if (nu == 1.0) return("Equidispersion regime ($\\nu=1.0$)")
  if (nu <= 1.5) return(sprintf("Underdispersion regime ($\\nu=%.1f$)", nu))
  return(sprintf("Strong underdispersion regime ($\\nu=%.1f$)", nu))
}

emit_latex <- function(summary_tab, lambdas, ns, nus) {
  lines <- c()
  for (nu in sort(nus)) {
    lines <- c(lines, paste0("\\multicolumn{6}{l}{\\textbf{", nu_regime_label(nu), "}} \\\\"))
    for (par in c("alpha", "lambda", "nu", "rho")) {
      true_val <- switch(par,
        alpha = 0.30, rho = 0.30, nu = nu,
        lambda = lambdas[[as.character(nu)]])
      cells <- sapply(ns, function(n) {
        key <- paste0("nu=", nu, "_n=", n)
        row <- summary_tab[summary_tab$cell == key, ]
        if (nrow(row) == 0) return("--")
        fmt(row[[paste0(par, "_MCmean")]], row[[paste0(par, "_MCsd")]])
      })
      symb <- switch(par, alpha = "\\alpha", lambda = "\\lambda", nu = "\\nu", rho = "\\rho")
      lines <- c(lines, sprintf("$%s$ & %.2f & %s & %s & %s & %s \\\\",
                                 symb, true_val, cells[1], cells[2], cells[3], cells[4]))
    }
    lines <- c(lines, "\\midrule")
  }
  writeLines(lines, "sim_table_body_ext.tex")
  cat("\nWrote sim_table_body_ext.tex -- rows for all regimes in NUS",
      "(paste the nu=0.2 and nu=2.5 blocks into a new appendix table,",
      "or splice all five regimes into an expanded Table 1, as you prefer).\n")
}

if (!QUICK_TEST) emit_latex(summary_tab, lambdas, NS, NUS)

cat("\nDone.\n")
