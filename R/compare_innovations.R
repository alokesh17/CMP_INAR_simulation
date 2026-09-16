############################################################################
# Model comparison: ZICMP-INAR(1) (this paper) vs. Poisson / Negative
# Binomial / Generalized Poisson innovations, all fit as ZI-INAR(1) models
# via the ZIHINAR1 package (S. Yang, CRAN), using the SAME EAIC/EBIC/DIC/
# WAIC1/WAIC2 criteria already defined in Section 3 of the draft.
#
# WHY THIS SCRIPT: run_simulation_grid.R checks PARAMETER RECOVERY under the
# true model (Table 1). This script checks MODEL SELECTION -- given data
# truly generated from ZICMP-INAR(1), do information criteria correctly
# prefer CMP over the nested/related alternatives, and by how much, across
# the three dispersion regimes? ZICMP nests neither ZIP, ZINB nor ZIGP
# exactly, but nu=1 in the CMP family reduces to the equidispersed Poisson
# case, so nu=1.0 is the hardest regime for CMP to "win" on parsimony
# grounds (ZIP has one fewer parameter) -- that is the honest, informative
# test, not a softball.
#
# HOW THE CRITERIA ARE MATCHED (verified against the ZIHINAR1 source,
# github.com/cran/ZIHINAR1, R/mod_sel_criteria.R and R/get_mol_sel.R):
#   EAIC  = mean over posterior draws of Stan's own `aic` generated quantity
#           (= -2*ll + 2*npar per draw); EBIC analogous with `bic`.
#   DIC   = -2*logphat + 2*pdic,  pdic = 2*(logphat - mean(ll))
#           where logphat = log-likelihood at the POSTERIOR-MEAN parameter
#           point estimate (their get_loglik(); for CMP this is
#           cmp_transition_loglik() below, an exact R port of the .stan
#           transformed-parameters block).
#   WAIC1/WAIC2 = the standard Gelman et al. (2014) pWAIC1/pWAIC2 forms,
#           computed from the T-1 pointwise (log-)likelihood vectors that
#           every one of these Stan files exposes as `lik`/`log_lik`.
# All four fitted models (ZICMP, ZINB, ZIGP: 4 params; ZIP: 3 params) use
# this same convention, so EAIC/EBIC/DIC/WAIC are directly comparable
# across models without any further adjustment.
#
# INSTALL (once)
#   install.packages("ZIHINAR1")
#   # or, if not yet on CRAN for your R version:
#   # remotes::install_github("fushengyy/ZIHINAR1")
#
# HOW TO USE
#   Put this file next to ZIINAR1-CMP-fast-reparam.stan and run:
#     Rscript compare_innovations.R
#   RUN_MODE controls cost the same way as in run_simulation_grid.R:
#     "quick" -- 1 regime, 2 reps  : confirms everything compiles and runs
#     "pilot" -- all 3 regimes, 10 reps, n=600 : a first honest read
#     "full"  -- all 3 regimes, 30 reps, n=600 : the numbers to report
#   Outputs: compare_summary_table.csv (mean EAIC/EBIC/DIC/WAIC1/WAIC2 per
#   model per regime) and compare_win_rate.csv (fraction of replicates in
#   which each model achieves the lowest WAIC2), plus a checkpointed
#   compare_results_partial.rds after every regime.
#
# NOTE: this is the ORIGINAL 3-regime/30-rep version. run_simulation_grid_
# extended.R / part4_resume.R in this same folder extend the design to 5
# regimes (adding nu=0.2 and nu=2.5) at R=10 -- that is the version whose
# numbers actually appear in the published draft's Table 2/3. Kept here for
# reference and as the simpler starting point if you want to re-derive the
# extended grid yourself.
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

if (RUN_MODE == "quick") {
  NUS <- c(1.0); R <- 2
  ITER <- 800; WARMUP <- 400; CHAINS <- 2
} else if (RUN_MODE == "pilot") {
  NUS <- c(0.5, 1.0, 1.5); R <- 10
  ITER <- 2000; WARMUP <- 1000; CHAINS <- 2
} else {
  NUS <- c(0.5, 1.0, 1.5); R <- 30
  ITER <- 2000; WARMUP <- 1000; CHAINS <- 2
}

# ---------------------------------------------------------------------
# 1. DGP (identical to run_simulation_grid.R)
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
# 2. R port of the CMP transition log-likelihood (exact match to the
#    `transformed parameters` block of ZIINAR1-CMP-fast-reparam.stan),
#    used only to evaluate logphat at a POINT ESTIMATE for DIC -- this is
#    the CMP analogue of ZIHINAR1's internal get_loglik().
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

# EAIC/EBIC/DIC/WAIC1/WAIC2 for the CMP fit -- exact same formulas as
# ZIHINAR1::get_mod_sel(), so all four models are judged identically.
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
# 3. Calibrate lambda per regime (same target mean as the main grid)
# ---------------------------------------------------------------------
lambdas <- setNames(sapply(NUS, solve_lambda_for_mean, target_mu = TARGET_MU), as.character(NUS))
cat("Calibrated lambdas (target CMP mean =", TARGET_MU, "):\n"); print(lambdas)

# ---------------------------------------------------------------------
# 4. Compile the CMP Stan model once; ZIP/ZINB/ZIGP compile lazily inside
#    ZIHINAR1::get_stanfit() the first time each is called (auto_write
#    caches the compiled .rds next to the package's own .stan files).
# ---------------------------------------------------------------------
cat("\nCompiling", CMP_STAN_FILE, "...\n")
cmp_mod <- stan_model(CMP_STAN_FILE)

# ---------------------------------------------------------------------
# 5. Main loop
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
               seed = 2000 + r, refresh = 0, control = list(adapt_delta = 0.9)),
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
  saveRDS(results, "compare_results_partial.rds")
  cat(sprintf("  elapsed so far: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

saveRDS(results, "compare_results_final.rds")

# ---------------------------------------------------------------------
# 6. Summarize: mean criteria per model per regime, and win rate
# ---------------------------------------------------------------------
summary_tab <- do.call(rbind, lapply(names(results), function(k) {
  df <- results[[k]]
  agg <- aggregate(cbind(EAIC, EBIC, DIC, WAIC1, WAIC2) ~ model, df, mean)
  cbind(regime = k, agg)
}))
write.csv(summary_tab, "compare_summary_table.csv", row.names = FALSE)

win_rate <- do.call(rbind, lapply(names(results), function(k) {
  df <- results[[k]]
  reps <- unique(df$rep)
  winners <- sapply(reps, function(rr) {
    sub <- df[df$rep == rr, ]
    sub$model[which.min(sub$WAIC2)]
  })
  data.frame(regime = k, t(prop.table(table(factor(winners, levels = c("ZICMP", "ZIP", "ZINB", "ZIGP"))))))
}))
write.csv(win_rate, "compare_win_rate.csv", row.names = FALSE)

cat("\n--- Mean model-selection criteria by regime (lower = better) ---\n")
print(summary_tab)
cat("\n--- Fraction of replicates where each model has the LOWEST WAIC2 ---\n")
print(win_rate)
cat("\nDone. See compare_summary_table.csv and compare_win_rate.csv.\n")
