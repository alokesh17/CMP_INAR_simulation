############################################################################
# Model comparison, EXTENDED: ZICMP-INAR(1) vs. Poisson / Negative Binomial /
# Generalized Poisson innovations (all via ZIHINAR1), across FIVE dispersion
# regimes -- the original nu in {0.5, 1.0, 1.5} plus the two new, stronger
# regimes nu in {0.2, 2.5} used in run_simulation_grid_extended.R /
# ZICMP_INAR_simulation_plot_predictive.R.
#
# This is compare_innovations.R with exactly one substantive change: NUS
# widened from c(0.5, 1.0, 1.5) to c(0.2, 0.5, 1.0, 1.5, 2.5). Everything
# else -- the DGP, the R port of the CMP transition log-likelihood used for
# DIC's point estimate, the EAIC/EBIC/DIC/WAIC1/WAIC2 formulas matched to
# ZIHINAR1's own get_mod_sel(), the win-rate computation -- is unchanged, so
# results for the original three regimes should reproduce what's already in
# the draft's Table 2 / Table 3 (mean-zero Monte Carlo noise aside).
#
# Rationale for the two additions, same as run_simulation_grid_extended.R:
#   nu = 0.2 -- markedly stronger overdispersion than nu=0.5, testing
#               whether ZICMP's edge over ZINB/ZIGP (or lack of it) holds up
#               further from equidispersion.
#   nu = 2.5 -- markedly stronger underdispersion than nu=1.5, the regime
#               where ZINB structurally cannot compete at all (its variance
#               is bounded below by its mean) -- worth checking whether
#               ZICMP's margin over ZIP/ZIGP widens, narrows, or is already
#               saturated at nu=1.5.
#
# INSTALL (once): install.packages("ZIHINAR1")  (or remotes::install_github
# if not yet on CRAN for your R version)
#
# HOW TO USE: identical to compare_innovations.R -- put this file next to
# ZIINAR1-CMP-fast-reparam.stan and run. RUN_MODE controls cost:
#   "quick" -- 1 regime (nu=2.5), 2 reps : confirms everything compiles/runs
#   "pilot" -- all 5 regimes, 10 reps, n=600 : a first honest read
#   "full"  -- all 5 regimes, 30 reps, n=600 : numbers to report
# Outputs are suffixed _ext so a run here won't overwrite results from the
# original 3-regime compare_innovations.R if both are run in the same
# directory: compare_summary_table_ext.csv, compare_win_rate_ext.csv,
# compare_results_partial_ext.rds / compare_results_final_ext.rds.
############################################################################

library(rstan)
library(ZIHINAR1)
library(COMPoissonReg)      # for rcmp(), same DGP as run_simulation_grid.R
rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))
set.seed(2026)

# ---------------------------------------------------------------------
# 0. Config
# ---------------------------------------------------------------------
RUN_MODE <- "pilot"    # <- "quick" | "pilot" | "full"

CMP_STAN_FILE <- "ZIINAR1-CMP-fast-reparam.stan"
M          <- 300
HYBRID_TOL <- 1e-6
ALPHA_TRUE <- 0.30
RHO_TRUE   <- 0.30
TARGET_MU  <- 3.0
N_FIXED    <- 600      # sample size for the model-comparison study

# Set EXTRA_ONLY <- TRUE to run just the two new regimes (cheaper if you
# already have compare_summary_table.csv/compare_win_rate.csv from the
# original 3-regime script and just want to append the extremes).
EXTRA_ONLY <- FALSE

if (RUN_MODE == "quick") {
  NUS <- c(2.5); R <- 2
  ITER <- 800; WARMUP <- 400; CHAINS <- 2
} else if (RUN_MODE == "pilot") {
  NUS <- if (EXTRA_ONLY) c(0.2, 2.5) else c(0.2, 0.5, 1.0, 1.5, 2.5)
  R <- 10
  ITER <- 2000; WARMUP <- 1000; CHAINS <- 2
} else {
  NUS <- if (EXTRA_ONLY) c(0.2, 2.5) else c(0.2, 0.5, 1.0, 1.5, 2.5)
  R <- 30
  ITER <- 2000; WARMUP <- 1000; CHAINS <- 2
}

# ---------------------------------------------------------------------
# 1. DGP (identical to run_simulation_grid.R / compare_innovations.R)
# ---------------------------------------------------------------------
thin_operator <- function(x, alpha) sum(rbinom(x, size = 1, prob = alpha))

simul_zinarCMP <- function(n, alpha, lambda, nu, zi_prob) {
  zinp_inar <- integer(n)
  zinp_inar[1] <- rcmp(1, lambda, nu)
  for (t in 2:n) {
    thinned <- thin_operator(zinp_inar[t - 1], alpha)
    innovation <- if (runif(1) < zi_prob) 0 else rcmp(1, lambda, nu)
    zinp_inar[t] <- thinned + innovation
  }
  zinp_inar
}

logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

# kmax raised from 400 to 600: nu=0.2's heavier tail needs more terms to
# converge to the same tolerance than the original nu in {0.5,1,1.5} did.
cmp_moments <- function(lambda, nu, kmax = 600) {
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
# 2. R port of the CMP transition log-likelihood (unchanged from
#    compare_innovations.R -- exact match to ZIINAR1-CMP-fast-reparam.stan's
#    transformed-parameters block), used to evaluate logphat at a point
#    estimate for DIC.
# ---------------------------------------------------------------------
cmp_log_Z_R <- function(lambda, nu, M, hybrid_tol, lgam) {
  log_lambda <- log(lambda)
  test <- exp(-log_lambda / nu)
  if (test < hybrid_tol) {
    nu * exp(log_lambda / nu) - ((nu - 1) / (2 * nu)) * log_lambda -
      ((nu - 1) / 2) * log(2 * pi) - 0.5 * log(nu)
  } else {
    r <- 0:M
    logsumexp(r * log_lambda - nu * lgam[r + 1])
  }
}

cmp_transition_loglik <- function(y, alpha, lambda, nu, rho, M = 300, hybrid_tol = 1e-6) {
  Tt <- length(y)
  lgam <- lgamma(0:(M + 1) + 1)
  log_Z <- cmp_log_Z_R(lambda, nu, M, hybrid_tol, lgam)
  log_lam <- log(lambda); log_rho <- log(rho); log1mrho <- log1p(-rho)
  out <- numeric(Tt - 1)
  for (t in 2:Tt) {
    yt <- y[t]; yt1 <- y[t - 1]; p <- min(yt1, yt)
    lbin0 <- dbinom(0, yt1, alpha, log = TRUE)
    lcmp0 <- yt * log_lam - nu * lgam[yt + 1] - log_Z
    lterm0 <- if (yt == 0) lbin0 + logsumexp(c(log_rho, log1mrho + lcmp0)) else
                            lbin0 + log1mrho + lcmp0
    if (p == 0) {
      out[t - 1] <- lterm0
    } else {
      lterms <- numeric(p + 1); lterms[1] <- lterm0
      for (j in 1:p) {
        lbinj <- dbinom(j, yt1, alpha, log = TRUE)
        diff  <- yt - j
        lcmpj <- diff * log_lam - nu * lgam[diff + 1] - log_Z
        lterms[j + 1] <- if (yt == j) lbinj + logsumexp(c(log_rho, log1mrho + lcmpj)) else
                                       lbinj + log1mrho + lcmpj
      }
      out[t - 1] <- logsumexp(lterms)
    }
  }
  out
}

get_mod_sel_cmp <- function(y, stan_fit, M = 300, hybrid_tol = 1e-6) {
  aic <- rstan::extract(stan_fit, pars = "aic")[[1]]; eaic <- mean(aic)
  bic <- rstan::extract(stan_fit, pars = "bic")[[1]]; ebic <- mean(bic)

  ph <- summary(stan_fit, pars = c("alpha", "rho", "lambda", "nu"))$summary
  logphat <- sum(cmp_transition_loglik(y, ph["alpha", "mean"], ph["lambda", "mean"],
                                        ph["nu", "mean"], ph["rho", "mean"],
                                        M, hybrid_tol))

  ll <- rstan::extract(stan_fit, pars = "ll")[[1]]
  pdic <- 2 * (logphat - mean(ll))
  dic  <- -2 * logphat + 2 * pdic

  loglik_mat <- rstan::extract(stan_fit, pars = "log_lik")[[1]][, 2:length(y), drop = FALSE]
  lik_mat    <- exp(loglik_mat)
  lppd   <- sum(log(colMeans(lik_mat)))
  pwaic1 <- 2 * sum(log(colMeans(lik_mat)) - colMeans(loglik_mat))
  pwaic2 <- sum(matrixStats::colVars(loglik_mat))
  waic1  <- -2 * (lppd - pwaic1)
  waic2  <- -2 * (lppd - pwaic2)

  data.frame(EAIC = eaic, EBIC = ebic, DIC = dic, WAIC1 = waic1, WAIC2 = waic2)
}

# ---------------------------------------------------------------------
# 3. Calibrate lambda per regime
# ---------------------------------------------------------------------
lambdas <- setNames(sapply(NUS, solve_lambda_for_mean, target_mu = TARGET_MU), as.character(NUS))
cat("Calibrated lambdas (target CMP mean =", TARGET_MU, "):\n"); print(lambdas)

# ---------------------------------------------------------------------
# 4. Compile the CMP Stan model once
# ---------------------------------------------------------------------
cat("\nCompiling", CMP_STAN_FILE, "...\n")
cmp_mod <- stan_model(CMP_STAN_FILE)

# ---------------------------------------------------------------------
# 5. Main loop (unchanged logic, just iterates over the wider NUS)
# ---------------------------------------------------------------------
results <- list()
t0 <- Sys.time()

for (nu_true in NUS) {
  lam_true <- lambdas[[as.character(nu_true)]]
  key <- paste0("nu=", nu_true)
  cat("\n===", key, "(n =", N_FIXED, ") ===\n")
  rows <- vector("list", R)

  for (r in seq_len(R)) {
    y <- simul_zinarCMP(N_FIXED, ALPHA_TRUE, lam_true, nu_true, RHO_TRUE)

    fit_cmp <- tryCatch(
      sampling(cmp_mod, data = list(T = N_FIXED, y = y, M = M,
                                     hybrid_tol = HYBRID_TOL, ff = 0),
               chains = CHAINS, iter = ITER, warmup = WARMUP,
               seed = 2000 + r, refresh = 0,
               control = list(adapt_delta = 0.95)),  # raised from 0.9, same
               # reasoning as run_simulation_grid_extended.R: the two new
               # regimes sit closer to the nu boundary of the CMP family
      error = function(e) { message("  ZICMP fit failed: ", conditionMessage(e)); NULL })

    fit_poi <- tryCatch(
      ZIHINAR1::get_stanfit(mod_type = "zi", distri = "poi", y = y,
                             n_pred = 0, chains = CHAINS, iter = ITER,
                             warmup = WARMUP, seed = 2000 + r),
      error = function(e) { message("  ZIP fit failed: ", conditionMessage(e)); NULL })

    fit_nb <- tryCatch(
      ZIHINAR1::get_stanfit(mod_type = "zi", distri = "nb", y = y,
                             n_pred = 0, chains = CHAINS, iter = ITER,
                             warmup = WARMUP, seed = 2000 + r),
      error = function(e) { message("  ZINB fit failed: ", conditionMessage(e)); NULL })

    fit_gp <- tryCatch(
      ZIHINAR1::get_stanfit(mod_type = "zi", distri = "gp", y = y,
                             n_pred = 0, chains = CHAINS, iter = ITER,
                             warmup = WARMUP, seed = 2000 + r),
      error = function(e) { message("  ZIGP fit failed: ", conditionMessage(e)); NULL })

    if (is.null(fit_cmp) || is.null(fit_poi) || is.null(fit_nb) || is.null(fit_gp)) next

    crit_cmp <- get_mod_sel_cmp(y, fit_cmp, M, HYBRID_TOL)
    crit_poi <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "poi", stan_fit = fit_poi)
    crit_nb  <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "nb",  stan_fit = fit_nb)
    crit_gp  <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "gp",  stan_fit = fit_gp)

    rows[[r]] <- rbind(
      cbind(model = "ZICMP", rep = r, crit_cmp),
      cbind(model = "ZIP",   rep = r, crit_poi),
      cbind(model = "ZINB",  rep = r, crit_nb),
      cbind(model = "ZIGP",  rep = r, crit_gp)
    )
    cat(sprintf("  rep %d/%d  EAIC[cmp,poi,nb,gp] = %.1f, %.1f, %.1f, %.1f  |  WAIC2[cmp,poi,nb,gp] = %.1f, %.1f, %.1f, %.1f\n",
                r, R, crit_cmp$EAIC, crit_poi$EAIC, crit_nb$EAIC, crit_gp$EAIC,
                crit_cmp$WAIC2, crit_poi$WAIC2, crit_nb$WAIC2, crit_gp$WAIC2))
  }

  results[[key]] <- do.call(rbind, rows[!sapply(rows, is.null)])
  saveRDS(results, "compare_results_partial_ext.rds")
  cat(sprintf("  elapsed so far: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

saveRDS(results, "compare_results_final_ext.rds")

# ---------------------------------------------------------------------
# 6. Summarize: mean criteria per model per regime, and win rate
# ---------------------------------------------------------------------
summary_tab <- do.call(rbind, lapply(names(results), function(k) {
  df <- results[[k]]
  agg <- aggregate(cbind(EAIC, EBIC, DIC, WAIC1, WAIC2) ~ model, df, mean)
  cbind(regime = k, agg)
}))
write.csv(summary_tab, "compare_summary_table_ext.csv", row.names = FALSE)

win_rate <- do.call(rbind, lapply(names(results), function(k) {
  df <- results[[k]]
  reps <- unique(df$rep)
  winners <- sapply(reps, function(rr) {
    sub <- df[df$rep == rr, ]
    sub$model[which.min(sub$WAIC2)]
  })
  data.frame(regime = k, t(prop.table(table(factor(winners, levels = c("ZICMP", "ZIP", "ZINB", "ZIGP"))))))
}))
write.csv(win_rate, "compare_win_rate_ext.csv", row.names = FALSE)

cat("\n--- Mean model-selection criteria by regime (lower = better) ---\n")
print(summary_tab)
cat("\n--- Fraction of replicates where each model has the LOWEST WAIC2 ---\n")
print(win_rate)
cat("\nDone. See compare_summary_table_ext.csv and compare_win_rate_ext.csv.\n",
    "If EXTRA_ONLY was TRUE, merge these two new regime rows with your\n",
    "existing compare_summary_table.csv / compare_win_rate.csv by hand,\n",
    "or just rerun with EXTRA_ONLY <- FALSE for all five regimes at once.\n")
