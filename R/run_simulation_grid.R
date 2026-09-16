############################################################################
# Simulation study driver for the ZICMP-INAR(1) model
# Default target is ZIINAR1-CMP-fast-reparam.stan: same likelihood as your
# ZIINAR1-CMP-fast.stan, but reparameterized to (mu_cmp, nu) with
# lambda := mu_cmp^nu, because the natural (lambda, nu) parameterization
# showed cor(lambda, nu) = 0.925 in posterior draws at n=100 -- a known CMP
# ridge, not a bug. Set STAN_FILE back to the original if you'd rather work
# around it a different way. Produces the table for Section 4.2
# (Table \ref{tab:sim_placeholder}) plus two comparators that never touch
# Stan (CLS for alpha; a moment-based estimator for lambda,nu given the true
# alpha,rho) and the posterior cor(lambda,nu) diagnostic for every fit.
#
# HOW TO USE
#   1. Put this file, ZIINAR1-CMP-fast.stan, and
#      ZIINAR1-CMP-fast-reparam.stan in the same directory.
#   2. RUN_MODE <- "quick": confirms the model compiles and one cell runs.
#   3. RUN_MODE <- "pilot": nu=1.0, n=600, 15 reps -- check that
#      lambda_nu_cor_mean has dropped well below the 0.925 seen before, and
#      that nu_MCmean sits near 1.0, before paying for the full grid.
#   4. RUN_MODE <- "full": the real 3x4 grid. Adjust R/ITER/WARMUP for your
#      compute budget, then run: Rscript run_simulation_grid.R
#   5. Results checkpoint after every (nu, n) cell to sim_results_partial.rds
#      so a long run can be killed and resumed without losing progress.
#   6. Final outputs: sim_summary_table.csv (numbers) and
#      sim_table_body.tex (ready-to-paste LaTeX rows for Table 1 in the draft).
############################################################################

library(rstan)
library(COMPoissonReg)   # for rcmp(), same as your original simulation code
rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))

set.seed(2026)

# ---------------------------------------------------------------------
# 0. Config -- adjust these for your machine / time budget
# ---------------------------------------------------------------------
# "quick" : one cell, 3 reps -- just confirms the model compiles and runs
# "pilot" : nu=1.0 only, n=600, 15 reps -- checks whether the lambda-nu
#           reparameterization actually fixed the ridge before you pay for
#           the full grid (cor(lambda,nu) was 0.925 under the natural
#           parameterization at n=100 -- re-check it here)
# "full"  : the real 3x4 grid used for Table 1
RUN_MODE <- "pilot"    # <- "quick" | "pilot" | "full"

# Use the reparameterized model (mu_cmp, nu) with lambda := mu_cmp^nu.
# The original natural-parameterization file (lambda, nu) showed
# cor(lambda, nu) = 0.925 in posterior draws -- a known CMP identifiability
# ridge, not a bug -- so this is now the default.
STAN_FILE   <- "ZIINAR1-CMP-fast-reparam.stan"
M           <- 300
HYBRID_TOL  <- 1e-6
FF          <- 0        # no forecast needed for the parameter-recovery table
ALPHA_TRUE  <- 0.30
RHO_TRUE    <- 0.30
TARGET_MU   <- 3.0      # common CMP innovation mean across dispersion regimes

if (RUN_MODE == "quick") {
  NUS <- c(1.0); NS <- c(100); R <- 3
  ITER <- 800; WARMUP <- 400; CHAINS <- 2
} else if (RUN_MODE == "pilot") {
  NUS <- c(1.0); NS <- c(600); R <- 15
  ITER <- 2000; WARMUP <- 1000; CHAINS <- 2
} else {
  NUS <- c(0.5, 1.0, 1.5)
  NS  <- c(100, 200, 400, 600)
  R   <- 30              # replications per (nu, n) cell -- raise if time allows
  ITER <- 2000; WARMUP <- 1000; CHAINS <- 2
}
QUICK_TEST <- (RUN_MODE != "full")   # controls whether the LaTeX table gets emitted

# ---------------------------------------------------------------------
# 1. Simulator (your own functions, unchanged)
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
# 2. Exact CMP moments (independent of Stan; used for lambda calibration
#    and for the moment-based comparator below)
# ---------------------------------------------------------------------
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

cmp_moments <- function(lambda, nu, kmax = 400) {
  k <- 0:kmax
  logw <- k * log(lambda) - nu * lgamma(k + 1)
  logZ <- logsumexp(logw)
  p <- exp(logw - logZ)
  mu <- sum(k * p)
  c(mean = mu, var = sum(k^2 * p) - mu^2)
}

solve_lambda_for_mean <- function(nu, target_mu, kmax = 400) {
  f <- function(lam) cmp_moments(lam, nu, kmax)["mean"] - target_mu
  uniroot(f, c(1e-6, 1e6), tol = 1e-8)$root
}

# ---------------------------------------------------------------------
# 3. Comparators, independent of the Stan fit
# ---------------------------------------------------------------------

# Conditional Least Squares estimator of alpha: closed-form, holds for ANY
# innovation law (Remark on ACF invariance in the theory note) -- a check
# on the Bayesian alpha estimate that does not touch the CMP likelihood.
cls_alpha <- function(y) {
  y1 <- y[-length(y)]; y2 <- y[-1]
  m1 <- mean(y1); m2 <- mean(y2)
  sum((y1 - m1) * (y2 - m2)) / sum((y1 - m1)^2)
}

# Oracle moment-matching estimator of (lambda, nu): plug the TRUE alpha and
# rho into the stationary mean/variance identities (Proposition 1 of the
# draft) and invert for (lambda, nu) from the empirical mean/variance alone.
# This isolates whether nu is recoverable from moments, with no MCMC involved.
mom_oracle_lambda_nu <- function(y, alpha_true, rho_true, kmax = 400) {
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
    r_lo <- resid(0.05); r_hi <- resid(8)
    if (!is.finite(r_lo) || !is.finite(r_hi) || sign(r_lo) == sign(r_hi)) return(NA_real_)
    uniroot(resid, c(0.05, 8), tol = 1e-6)$root
  }, error = function(e) NA_real_)
  if (is.na(nu_hat)) return(c(lambda = NA_real_, nu = NA_real_))
  lam_hat <- tryCatch(solve_lambda_for_mean(nu_hat, muU_target, kmax), error = function(e) NA_real_)
  c(lambda = lam_hat, nu = nu_hat)
}

# ---------------------------------------------------------------------
# 4. Calibrate lambda per regime (common innovation mean across nu)
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
                 control = list(adapt_delta = 0.9)),
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
    saveRDS(results, "sim_results_partial.rds")
    cat(sprintf("  elapsed so far: %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }
}

saveRDS(results, "sim_results_final.rds")
cat("\nTOTAL TIME (min):", as.numeric(difftime(Sys.time(), t_start, units = "mins")), "\n")

# ---------------------------------------------------------------------
# 7. Summarize: Monte Carlo mean and SD of the posterior means, across reps
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
write.csv(summary_tab, "sim_summary_table.csv", row.names = FALSE)
cat("\n--- Summary table (also written to sim_summary_table.csv) ---\n")
print(summary_tab)

# ---------------------------------------------------------------------
# 8. Emit LaTeX rows matching Table \ref{tab:sim_placeholder}
# ---------------------------------------------------------------------
fmt <- function(m, s) sprintf("%.3f (%.3f)", m, s)

emit_latex <- function(summary_tab, lambdas, ns) {
  regimes <- list(
    list(nu = 0.5, label = "Overdispersion regime ($\\nu=0.5$)"),
    list(nu = 1.0, label = "Equidispersion regime ($\\nu=1.0$)"),
    list(nu = 1.5, label = "Underdispersion regime ($\\nu=1.5$)")
  )
  lines <- c()
  for (reg in regimes) {
    lines <- c(lines, paste0("\\multicolumn{6}{l}{\\textbf{", reg$label, "}} \\\\"))
    for (par in c("alpha", "lambda", "nu", "rho")) {
      true_val <- switch(par,
        alpha = 0.30, rho = 0.30, nu = reg$nu,
        lambda = lambdas[[as.character(reg$nu)]])
      cells <- sapply(ns, function(n) {
        key <- paste0("nu=", reg$nu, "_n=", n)
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
  writeLines(lines, "sim_table_body.tex")
  cat("\nWrote sim_table_body.tex -- paste these rows into Table 1 of the draft.\n")
}

if (!QUICK_TEST) emit_latex(summary_tab, lambdas, NS)

cat("\nDone.\n")
