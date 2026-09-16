############################################################################
# ZICMP-INAR(1): simulation + plots + predictive (PARTS 1-3 only)
#   PART 1 -- simulation: parameter-recovery grid, 5 dispersion regimes
#             (original nu in {0.5,1.0,1.5} + 2 new extremes nu in {0.2,2.5})
#   PART 2 -- plots: recovery, bias, model-comparison win rate (base R)
#   PART 3 -- predictive analysis: qoi + y_pred forecast plot, this
#             model's real parameters (alpha, lambda, nu, rho), not the
#             "lambda1, lambda2" naming from a different model
#
#   PART 4 (model comparison vs. ZIP/ZINB/ZIGP via ZIHINAR1) has been
#   split OUT into its own standalone script, part4_resume.R -- run
#   that separately (`Rscript part4_resume.R` from a plain Terminal)
#   once Part 1 below has produced cmp_cache.rds. Nothing in this file
#   depends on Part 4; Part 2's win-rate plot below falls back to
#   placeholder numbers if `win_rate` isn't in scope, same as before.
#
# Built around ZIINAR1-CMP-fast-reparam.stan (parameters: alpha, mu_cmp,
# nu, rho; lambda := mu_cmp^nu recovered in transformed parameters;
# generated quantities: y_pred, log_lik, ll, aic, bic).
#
# Nothing fits or plots by itself just from sourcing this file: Part 1's
# grid runs once you set RUN_MODE and execute; Part 2 only plots once
# Part 1 has populated `summary_tab`; Part 3's worked example is gated
# with `if (FALSE)`.
############################################################################

library(rstan)
library(COMPoissonReg)   # for rcmp()
library(parallel)        # base package, no install needed
rstan_options(auto_write = TRUE)

# ---------------------------------------------------------------------
# PARALLELISM: reps within a cell are independent, so they're run
# concurrently across cores via parallel::mclapply (fork-based) --
# Mac/Linux only, since Windows has no fork(). Each individual
# sampling() call is forced to run its own chains SEQUENTIALLY
# (mc.cores=1 below) so cores aren't double-booked: parallelism happens
# ACROSS reps instead, which uses the machine far more fully than the
# old setup (which only ever kept 2 cores busy, one fit at a time).
# ---------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")   # gives mclapply's forked workers statistically
                            # independent RNG streams (R's own recommended
                            # setup for parallel simulation -- see ?mclapply)
set.seed(2026)

HAS_FORK  <- .Platform$OS.type == "unix"   # mclapply needs Mac/Linux (fork())
N_WORKERS <- if (HAS_FORK) max(1, parallel::detectCores() - 1) else 1
if (!HAS_FORK) {
  message("Note: fork-based parallel reps (parallel::mclapply) need Mac/Linux; ",
          "Windows detected here, so reps will run one at a time below. For ",
          "real speedup on Windows, run this script inside WSL2.")
} else {
  message("Parallel reps enabled: ", N_WORKERS, " worker(s) (", parallel::detectCores(),
          " cores detected, 1 held back for the OS).")
}
options(mc.cores = 1)   # forces each sampling() call's own chains to run
                         # sequentially -- see comment block above

STAN_FILE <- "ZIINAR1-CMP-fast-reparam.stan"   # samples (alpha, mu_cmp, nu, rho),
                                                # lambda := mu_cmp^nu in
                                                # transformed parameters

RUN_MODE <- "full_r10" # <- "quick" | "pilot" | "full_r10" | "full"  (drives Part 1)
                        #    "full_r10": every (nu,n) cell, R=10 reps (requested tradeoff)
                        #    "full":     every (nu,n) cell, R=30 reps (original, most expensive)


############################################################################
# PART 1 -- SIMULATION GRID (parameter recovery)
############################################################################

# ---------------------------------------------------------------------
# 1.0 Config
# ---------------------------------------------------------------------
M           <- 300
HYBRID_TOL  <- 1e-6
FF          <- 0        # no forecast needed for the parameter-recovery grid;
                         # PART 3 below uses its own ff > 0 fit
ALPHA_TRUE  <- 0.30
RHO_TRUE    <- 0.30
TARGET_MU   <- 3.0      # common CMP innovation mean held across regimes

if (RUN_MODE == "quick") {
  NUS <- c(2.5); NS <- c(100); R <- 3
  ITER <- 800; WARMUP <- 400; CHAINS <- 2
} else if (RUN_MODE == "pilot") {
  NUS <- c(0.2); NS <- c(600); R <- 15
  ITER <- 2000; WARMUP <- 1000; CHAINS <- 2
} else if (RUN_MODE == "full_r10") {
  # Full (nu, n) grid -- all 5 regimes x all 4 sample sizes -- but with
  # R=10 reps/cell instead of 30, to cut cost roughly 3x while still
  # covering every combination. R=10 also matches the pilot rep count
  # already quoted in the draft's Table 2/3 captions.
  NUS <- c(0.2, 0.5, 1.0, 1.5, 2.5)
  NS  <- c(100, 200, 400, 600)
  R   <- 10
  ITER <- 2000; WARMUP <- 1000; CHAINS <- 2
} else {
  NUS <- c(0.2, 0.5, 1.0, 1.5, 2.5)   # original 3 regimes + 2 new extremes:
                                        # nu=0.2 (strong overdispersion),
                                        # nu=2.5 (strong underdispersion)
  NS  <- c(100, 200, 400, 600)
  R   <- 30
  ITER <- 2000; WARMUP <- 1000; CHAINS <- 2
}
QUICK_TEST <- !(RUN_MODE %in% c("full", "full_r10"))

# ---------------------------------------------------------------------
# 1.1 Simulator
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
# 1.2 Exact CMP moments (kmax raised to 600: nu=0.2's heavier tail needs
#     more terms to converge than the original nu in {0.5,1,1.5} did)
# ---------------------------------------------------------------------
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

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
# 1.2b R port of the CMP transition log-likelihood (exact match to
#      ZIINAR1-CMP-fast-reparam.stan's transformed-parameters block) and
#      the EAIC/EBIC/DIC/WAIC1/WAIC2 model-selection criteria for a CMP
#      fit -- exact same formulas as ZIHINAR1::get_mod_sel(), so all four
#      models (CMP + Part 4's ZIP/ZINB/ZIGP, run separately in
#      part4_resume.R) are judged identically.
#      This grid loop below computes and stashes these numbers at fit
#      time into cmp_cache -- part4_resume.R reuses them instead of
#      re-simulating + re-fitting CMP from scratch.
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
# 1.3 Comparators independent of the Stan fit
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
    r_lo <- resid(0.02); r_hi <- resid(12)
    if (!is.finite(r_lo) || !is.finite(r_hi) || sign(r_lo) == sign(r_hi)) return(NA_real_)
    uniroot(resid, c(0.02, 12), tol = 1e-6)$root
  }, error = function(e) NA_real_)
  if (is.na(nu_hat)) return(c(lambda = NA_real_, nu = NA_real_))
  lam_hat <- tryCatch(solve_lambda_for_mean(nu_hat, muU_target, kmax), error = function(e) NA_real_)
  c(lambda = lam_hat, nu = nu_hat)
}

# ---------------------------------------------------------------------
# 1.4 Calibrate lambda per regime
# ---------------------------------------------------------------------
lambdas <- setNames(sapply(NUS, solve_lambda_for_mean, target_mu = TARGET_MU), as.character(NUS))
cat("Calibrated lambdas (target CMP mean =", TARGET_MU, "):\n"); print(lambdas)
for (nu in NUS) {
  mv <- cmp_moments(lambdas[[as.character(nu)]], nu)
  cat(sprintf("  nu=%.2f  lambda=%.4f  mean=%.3f  var=%.3f  d_U=%.3f\n",
              nu, lambdas[[as.character(nu)]], mv["mean"], mv["var"], mv["var"] / mv["mean"]))
}

# ---------------------------------------------------------------------
# 1.5 Compile once (part4_resume.R recompiles its own copy of `mod`
#     separately -- it runs as its own R session)
# ---------------------------------------------------------------------
cat("\nCompiling", STAN_FILE, "...\n")
mod <- stan_model(STAN_FILE)

# ---------------------------------------------------------------------
# 1.6 Main grid -- fills `results` (raw per-replicate draws) used by
#     both the LaTeX table emitter and the plots in PART 2 below.
#     NOTE: the reparam Stan file's `parameters` block has `mu_cmp`, not
#     `lambda` -- but `lambda` is still available to pull posterior
#     summaries from because it's declared in `transformed parameters`
#     (lambda := mu_cmp^nu), so `pars = c("alpha","lambda","nu","rho")`
#     below works unchanged against ZIINAR1-CMP-fast-reparam.stan.
# ---------------------------------------------------------------------
results <- list()
# cmp_cache stores, per (nu,n) cell and rep, the exact y that was
# simulated and this fit's EAIC/EBIC/DIC/WAIC1/WAIC2 -- part4_resume.R
# reuses these, so the CMP model is never re-simulated + re-fit from
# scratch there (only ZIP/ZINB/ZIGP are fit fresh in that script).
cmp_cache <- list()
t_start <- Sys.time()

# One rep's full unit of work -- called in parallel across reps via
# mclapply below. Returns NULL on any failure (sampling error or the
# fit simply didn't happen); returns list(row=..., cache=...) otherwise.
# Must be self-contained: everything it uses (mod, ALPHA_TRUE, M,
# HYBRID_TOL, FF, CHAINS, ITER, WARMUP, cls_alpha, mom_oracle_lambda_nu,
# get_mod_sel_cmp, RHO_TRUE) is inherited by the forked child processes
# automatically (fork copies the parent's whole memory), so nothing
# needs to be passed in explicitly beyond what varies per call.
fit_one_rep_part1 <- function(r, n, nu_true, lam_true) {
  y <- simul_zinarCMP(n, ALPHA_TRUE, lam_true, nu_true, RHO_TRUE)

  fit <- tryCatch(
    sampling(mod,
             data = list(T = n, y = y, M = M, hybrid_tol = HYBRID_TOL, ff = FF),
             chains = CHAINS, iter = ITER, warmup = WARMUP,
             seed = 1000 + r, refresh = 0,
             control = list(adapt_delta = 0.95)),
    error = function(e) { message("  [rep ", r, "] sampling failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(fit)) return(NULL)

  s <- summary(fit, pars = c("alpha", "lambda", "nu", "rho"))$summary
  rhat_ok <- all(s[, "Rhat"] < 1.05, na.rm = TRUE)

  post <- rstan::extract(fit, pars = c("lambda", "nu"))
  lamnu_cor <- suppressWarnings(cor(post$lambda, post$nu))

  a_cls <- cls_alpha(y)
  mom <- mom_oracle_lambda_nu(y, ALPHA_TRUE, RHO_TRUE)

  # Cache this fit's model-selection criteria for part4_resume.R to
  # reuse -- if this fails for some numerical reason, just skip caching
  # (that script falls back to fitting CMP fresh for this rep).
  crit_cmp <- tryCatch(get_mod_sel_cmp(y, fit, M, HYBRID_TOL), error = function(e) NULL)

  row <- data.frame(
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
  cat(sprintf("  rep %d/%d  Rhat_ok=%s  alpha=%.3f  lambda=%.3f  nu=%.3f  rho=%.3f  cor(lam,nu)=%.3f\n",
              r, R, rhat_ok, s["alpha", "mean"], s["lambda", "mean"],
              s["nu", "mean"], s["rho", "mean"], lamnu_cor))

  list(row = row, cache = list(y = y, crit_cmp = crit_cmp))
}

ok_result <- function(x) !is.null(x) && !inherits(x, "try-error")

for (nu_true in NUS) {
  lam_true <- lambdas[[as.character(nu_true)]]
  for (n in NS) {
    key <- paste0("nu=", nu_true, "_n=", n)
    cat("\n===", key, "=== (", N_WORKERS, "worker(s) in parallel; output below may interleave)\n")

    rep_results <- mclapply(seq_len(R), fit_one_rep_part1, n = n, nu_true = nu_true,
                             lam_true = lam_true, mc.cores = N_WORKERS, mc.preschedule = FALSE)

    rows <- lapply(rep_results, function(x) if (ok_result(x)) x$row else NULL)
    cmp_cache[[key]] <- lapply(rep_results, function(x) if (ok_result(x)) x$cache else NULL)

    results[[key]] <- do.call(rbind, rows[!sapply(rows, is.null)])
    saveRDS(results, "sim_results_partial.rds")
    saveRDS(cmp_cache, "cmp_cache_partial.rds")   # checkpointed every cell now, not just at the end
    cat(sprintf("  elapsed so far: %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }
}
saveRDS(results, "sim_results_final.rds")
saveRDS(cmp_cache, "cmp_cache.rds")
cat("\nTOTAL TIME (min):", as.numeric(difftime(Sys.time(), t_start, units = "mins")), "\n")

# ---------------------------------------------------------------------
# 1.7 Summarize + write CSV/LaTeX (generic over however many regimes
#     are in NUS)
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
      true_val <- switch(par, alpha = 0.30, rho = 0.30, nu = nu, lambda = lambdas[[as.character(nu)]])
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
  writeLines(lines, "sim_table_body.tex")
  cat("\nWrote sim_table_body.tex\n")
}
if (!QUICK_TEST) emit_latex(summary_tab, lambdas, NS, NUS)


############################################################################
# PART 2 -- PLOTS (parameter recovery, bias, model-comparison win rate)
#           Uses `summary_tab` and `NUS`/`NS` from Part 1 directly --
#           no CSV re-read needed since it's the same script/session.
#           Nothing is plotted until you actually run Part 1 above (or
#           point `summary_tab` at a saved sim_summary_table.csv via
#           read.csv() instead).
#
#           The win-rate panel reads `win_rate`, which now only ever
#           comes from part4_resume.R's compare_win_rate.csv (that
#           script writes it once all 20 cells are done). Load it here
#           with: win_rate <- read.csv("compare_win_rate.csv") -- if
#           you skip that, plot_winrate() below falls back to the
#           original 3-regime placeholder numbers, same as always.
############################################################################

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; AQUA <- "#1baf7a"
YELLOW <- "#eda100"; PURPLE <- "#a259d9"
regime_colors <- c(BLUE, ORANGE, AQUA, YELLOW, PURPLE)

build_regime_list <- function(summary_tab, n_vals) {
  tab <- summary_tab
  tab$nu <- as.numeric(sub("nu=([^_]+)_n=.*", "\\1", tab$cell))
  tab$n  <- as.numeric(sub(".*_n=", "", tab$cell))
  nus <- sort(unique(tab$nu))
  out <- list()
  for (i in seq_along(nus)) {
    nu <- nus[i]
    sub_tab <- tab[tab$nu == nu, ]
    sub_tab <- sub_tab[match(n_vals, sub_tab$n), ]
    out[[as.character(nu)]] <- list(
      color = regime_colors[((i - 1) %% length(regime_colors)) + 1],
      true_nu = nu, true_lambda = lambdas[[as.character(nu)]],
      alpha  = list(mean = sub_tab$alpha_MCmean,  sd = sub_tab$alpha_MCsd),
      lambda = list(mean = sub_tab$lambda_MCmean, sd = sub_tab$lambda_MCsd),
      nu     = list(mean = sub_tab$nu_MCmean,     sd = sub_tab$nu_MCsd),
      rho    = list(mean = sub_tab$rho_MCmean,    sd = sub_tab$rho_MCsd)
    )
  }
  out
}

param_labels <- c(alpha = "hat(alpha)", lambda = "hat(lambda)", nu = "hat(nu)", rho = "hat(rho)")
panel_order  <- c("nu", "alpha", "lambda", "rho")

plot_recovery <- function(data_list, n_vals, file_stem, title_suffix = "") {
  regimes <- names(data_list); n_regimes <- length(regimes)
  jitters <- (seq_len(n_regimes) - (n_regimes + 1) / 2) * 0.03 * diff(range(n_vals))
  make_plot <- function() {
    op <- par(mfrow = c(2, 2), mar = c(4, 4.5, 2.5, 1), oma = c(4.5, 0, 2, 0)); on.exit(par(op))
    for (pname in panel_order) {
      all_means <- unlist(lapply(data_list, function(d) d[[pname]]$mean))
      all_sds   <- unlist(lapply(data_list, function(d) d[[pname]]$sd))
      ylim <- range(c(all_means - all_sds, all_means + all_sds), na.rm = TRUE)
      plot(NA, xlim = range(n_vals), ylim = ylim, xaxt = "n",
           xlab = "Sample size n", ylab = paste0(param_labels[pname], " (posterior mean +/- SD)"),
           main = param_labels[pname])
      axis(1, at = n_vals)
      for (i in seq_along(regimes)) {
        reg <- regimes[i]; d <- data_list[[reg]]; x <- n_vals + jitters[i]
        abline(h = if (pname == "nu") d$true_nu else if (!is.na(d$true_lambda) && pname == "lambda") d$true_lambda else NA,
               col = d$color, lty = 2, lwd = 1)
        arrows(x, d[[pname]]$mean - d[[pname]]$sd, x, d[[pname]]$mean + d[[pname]]$sd,
               angle = 90, code = 3, length = 0.03, col = d$color, lwd = 1.2)
        lines(x, d[[pname]]$mean, col = d$color, lwd = 1.4)
        points(x, d[[pname]]$mean, col = d$color, pch = 16, cex = 1)
      }
    }
    legend_labels <- sapply(regimes, function(r) sprintf("nu = %s", r))
    legend_cols <- sapply(data_list, function(d) d$color)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", legend = legend_labels, col = legend_cols, lwd = 1.4, pch = 16,
           horiz = TRUE, bty = "n", inset = c(0, -0.01), xpd = TRUE)
    mtext(paste0("Parameter recovery across dispersion regimes and sample sizes", title_suffix),
          outer = TRUE, cex = 1.1, line = 0.5)
  }
  pdf(paste0(file_stem, ".pdf"), width = 8.5, height = 7); make_plot(); dev.off()
  png(paste0(file_stem, ".png"), width = 2000, height = 1650, res = 220); make_plot(); dev.off()
}

plot_bias <- function(data_list, n_vals, file_stem) {
  regimes <- names(data_list); n_regimes <- length(regimes)
  jitters <- (seq_len(n_regimes) - (n_regimes + 1) / 2) * 0.03 * diff(range(n_vals))
  make_plot <- function() {
    op <- par(mfrow = c(2, 2), mar = c(4, 4.5, 2.5, 1), oma = c(4.5, 0, 2, 0)); on.exit(par(op))
    for (pname in panel_order) {
      biases <- lapply(data_list, function(d) {
        true_val <- if (pname == "nu") d$true_nu else if (pname %in% c("alpha", "rho")) 0.30 else d$true_lambda
        d[[pname]]$mean - true_val
      })
      ylim <- range(unlist(biases), na.rm = TRUE)
      plot(NA, xlim = range(n_vals), ylim = ylim, xaxt = "n",
           xlab = "Sample size n", ylab = paste0("Bias: ", param_labels[pname], " - true value"),
           main = param_labels[pname])
      axis(1, at = n_vals); abline(h = 0, col = "grey40", lwd = 1)
      for (i in seq_along(regimes)) {
        reg <- regimes[i]; d <- data_list[[reg]]; x <- n_vals + jitters[i]
        lines(x, biases[[reg]], col = d$color, lwd = 1.4)
        points(x, biases[[reg]], col = d$color, pch = 16, cex = 1)
      }
    }
    legend_labels <- sapply(regimes, function(r) sprintf("nu = %s", r))
    legend_cols <- sapply(data_list, function(d) d$color)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", legend = legend_labels, col = legend_cols, lwd = 1.4, pch = 16,
           horiz = TRUE, bty = "n", inset = c(0, -0.01), xpd = TRUE)
    mtext("Estimation bias (posterior mean - true value) by dispersion regime and n",
          outer = TRUE, cex = 1.1, line = 0.5)
  }
  pdf(paste0(file_stem, ".pdf"), width = 8.5, height = 7); make_plot(); dev.off()
  png(paste0(file_stem, ".png"), width = 2000, height = 1650, res = 220); make_plot(); dev.off()
}

# Win-rate plot: reads directly from a `win_rate` data.frame (load one
# with win_rate <- read.csv("compare_win_rate.csv") once part4_resume.R
# has finished, then rerun this section) -- falls back to the original
# 3-regime numbers already in the draft
# (Table~\ref{tab:model_comparison_winrate}) if you haven't loaded one,
# so this still produces a figure either way.
plot_winrate <- function(win_rate_df = NULL) {
  models_wr  <- c("ZICMP", "ZIP", "ZINB", "ZIGP")
  model_colors <- c(BLUE, ORANGE, AQUA, YELLOW)
  if (is.null(win_rate_df)) {
    regimes_wr <- c("Overdispersion\n(nu=0.5)", "Equidispersion\n(nu=1.0)", "Underdispersion\n(nu=1.5)")
    winrate <- rbind(c(0.70, 0.00, 0.30, 0.00), c(0.00, 0.90, 0.10, 0.00), c(0.90, 0.10, 0.00, 0.00))
  } else {
    # win_rate_df may be old-style (just a "regime" column, e.g. "nu=0.5",
    # one row per regime) or the new-style (separate "nu"/"n" columns,
    # possibly several n's per nu) -- prefer the "nu" column when present
    # so ordering/labels don't depend on parsing "nu=0.5_n=600" strings.
    nu_num <- if (!is.null(win_rate_df$nu)) as.numeric(win_rate_df$nu) else as.numeric(sub("nu=", "", win_rate_df$regime))
    ord <- order(nu_num)
    win_rate_df <- win_rate_df[ord, ]
    regimes_wr <- if (!is.null(win_rate_df$nu)) sprintf("nu=%s", win_rate_df$nu) else win_rate_df$regime
    winrate <- as.matrix(win_rate_df[, models_wr])
  }
  colnames(winrate) <- models_wr; rownames(winrate) <- regimes_wr
  bp <- barplot(t(winrate), beside = TRUE, col = model_colors, ylim = c(0, 1.05),
                ylab = "Fraction of replicates with lowest WAIC2",
                names.arg = regimes_wr, legend.text = models_wr,
                args.legend = list(x = "top", horiz = TRUE, bty = "n", inset = c(0, -0.08)),
                main = "Which model wins, by dispersion regime")
  for (i in seq_len(nrow(winrate))) for (j in seq_len(ncol(winrate))) {
    v <- winrate[i, j]; if (!is.na(v) && v > 0) text(bp[j, i], v + 0.03, sprintf("%.2f", v), cex = 0.75)
  }
}

# Run Part 2 (only meaningful once Part 1 has populated `summary_tab`):
# Safety: if you're re-pasting this into an R session that already has
# an OLD/stale `win_rate` object lying around (e.g. from an earlier
# partial run), drop it here so this always falls back cleanly to the
# placeholder win-rate numbers instead of crashing on a mismatched
# object shape.
if (exists("win_rate")) rm(win_rate)
if (exists("summary_tab") && nrow(summary_tab) > 0) {
  n_vals <- NS
  data_all <- build_regime_list(summary_tab, n_vals)
  plot_recovery(data_all, n_vals, "fig_recovery")
  plot_bias(data_all, n_vals, "fig_bias")
  pdf("fig_modelcomparison.pdf", width = 7.5, height = 4.2)
  plot_winrate(if (exists("win_rate")) win_rate else NULL)
  dev.off()
  png("fig_modelcomparison.png", width = 1900, height = 1050, res = 220)
  plot_winrate(if (exists("win_rate")) win_rate else NULL)
  dev.off()
  cat("\nWrote fig_recovery.{pdf,png}, fig_bias.{pdf,png}, fig_modelcomparison.{pdf,png}\n")
}


############################################################################
# PART 3 -- PREDICTIVE ANALYSIS
#           qoi + print(fit, pars=qoi) + y_pred-based forecast plot,
#           adapted to this model's real parameters/generated quantities.
############################################################################

qoi <- c("alpha", "lambda", "nu", "rho", "aic", "bic")
# print(fit, pars = qoi)   # once you have a fitted `fit` in scope

plot_predictive <- function(fit, y, ff, probs = c(0.1, 0.9),
                             series_name = "", model_name = "ZICMP", ylim = NULL) {
  n <- length(y)
  stopifnot(ff > 0, ff < n)

  # y_pred has ff+1 entries: y_pred[1] = y[T] (last TRAINING point, an
  # anchor, not a forecast); y_pred[2:(ff+1)] are the ff chained
  # one-step-ahead forecasts, matching y[(n-ff+1):n].
  fitPred   <- summary(fit, pars = "y_pred", probs = probs)$summary
  fitPredM  <- fitPred[-1, "mean"]
  loName    <- grep(paste0(probs[1] * 100, "%"), colnames(fitPred), value = TRUE)[1]
  hiName    <- grep(paste0(probs[2] * 100, "%"), colnames(fitPred), value = TRUE)[1]
  fitPredLo <- fitPred[-1, loName]
  fitPredHi <- fitPred[-1, hiName]

  y_obs <- y[(n - ff + 1):n]
  if (is.null(ylim)) {
    ylim <- range(c(y_obs, fitPredM, fitPredLo, fitPredHi), na.rm = TRUE)
    ylim[2] <- ylim[2] * 1.05
  }

  plot(y_obs, type = "o", pch = 16, ylim = ylim,
       xlab = "Held-out time index", ylab = "Count",
       main = sprintf("%s%s one-step-ahead forecasts (%d%% band)",
                       model_name, if (nzchar(series_name)) paste0(" -- ", series_name) else "",
                       round(100 * (probs[2] - probs[1]))))
  lines(fitPredM, col = "#2a78d6", lwd = 1.8)
  lines(fitPredLo, col = "#2a78d6", lty = 2)
  lines(fitPredHi, col = "#2a78d6", lty = 2)
  legend("topleft", legend = c("Observed", "Posterior predictive mean",
                                sprintf("%d%%-%d%% interval", 100 * probs[1], 100 * probs[2])),
         col = c("black", "#2a78d6", "#2a78d6"), lty = c(1, 1, 2), pch = c(16, NA, NA),
         bty = "n", cex = 0.85)

  invisible(data.frame(t = seq_len(ff), observed = y_obs,
                        pred_mean = fitPredM, pred_lo = fitPredLo, pred_hi = fitPredHi))
}

# Worked example: fit one real (or simulated) series with a genuine
# holdout (ff > 0) and run the predictive plot. Gated with `if (FALSE)`
# so sourcing this script doesn't try to fit anything on its own --
# fill in your own series and flip to TRUE, or just copy the body out.
if (FALSE) {
  y_full <- scan("sexoffences.txt")      # replace with your real loader
  n      <- length(y_full)
  ff     <- round(0.20 * n)              # same 20% holdout as forecast_evaluation.R
  n_tr   <- n - ff
  y_train <- y_full[1:n_tr]

  fit_pred <- sampling(mod,              # reuses the `mod` compiled in Part 1
                        data = list(T = n_tr, y = y_train, M = M,
                                    hybrid_tol = HYBRID_TOL, ff = ff),
                        chains = 4, iter = 2000, warmup = 1000, seed = 1,
                        control = list(adapt_delta = 0.95))

  print(fit_pred, pars = qoi)

  pred_tab    <- plot_predictive(fit_pred, y_full, ff, probs = c(0.1, 0.9),
                                  series_name = "sex offenses")
  pred_tab_95 <- plot_predictive(fit_pred, y_full, ff, probs = c(0.025, 0.975),
                                  series_name = "sex offenses")
}

cat("\nDone. Part 1 fits the recovery grid; Part 2 plots off `summary_tab`",
    "(load win_rate <- read.csv(\"compare_win_rate.csv\") after part4_resume.R",
    "has run and rerun Part 2 to get the real win-rate panel); Part 3's",
    "plot_predictive()/qoi are ready for any fit with ff > 0. Part 4 (model",
    "comparison vs ZIP/ZINB/ZIGP) now lives entirely in part4_resume.R.\n")
