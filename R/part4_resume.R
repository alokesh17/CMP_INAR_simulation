############################################################################
# PART 4 RESUME -- standalone script (safe to run on its own)
#
# WHY THIS EXISTS: the last live run had `cmp_cache` reloaded as an
# EMPTY list (length 0, should have been 20) and `mod` (the compiled
# Stan model) never recompiled after an R session restart. That made
# every rep in 4 cells (nu=1_n=600, nu=1.5_n=100/200/400) silently fail
# ("0 reused, 0 fresh CMP" in ~12 sec each). This script fixes both
# root causes explicitly (reload cmp_cache from the REAL file on disk +
# recompile mod every time it's run) and resumes Part 4 exactly where
# it left off, with a defensive check that drops any cell whose row
# count doesn't look complete (so a half-failed cell is always retried,
# never silently treated as "done").
#
# HOW TO RUN: from a plain Terminal (NOT pasted into the RStudio
# console -- that's what caused the earlier fork stall), in the same
# folder as your .stan/.rds files:
#     Rscript part4_resume.R
# This can be re-run safely any number of times -- already-complete
# cells are skipped every time.
############################################################################

library(rstan)
library(COMPoissonReg)
library(ZIHINAR1)
library(parallel)
rstan_options(auto_write = TRUE)

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)

# --- cores: dialed down from detectCores()-1 on purpose. Bump this up
#     (e.g. to detectCores()-1) once you've confirmed a plain-Terminal
#     Rscript run is stable on your machine; 3 is a conservative start. ---
HAS_FORK  <- .Platform$OS.type == "unix"
N_WORKERS <- if (HAS_FORK) 3 else 1
options(mc.cores = 1)
cat("Workers for this run:", N_WORKERS, "\n")

STAN_FILE <- "ZIINAR1-CMP-fast-reparam.stan"

# ---------------------------------------------------------------------
# Config -- MUST match what generated cmp_cache.rds (RUN_MODE="full_r10")
# ---------------------------------------------------------------------
M          <- 300
HYBRID_TOL <- 1e-6
ALPHA_TRUE <- 0.30
RHO_TRUE   <- 0.30
TARGET_MU  <- 3.0

COMPARE_EXTRA_ONLY <- FALSE
CMP_NUS <- c(0.2, 0.5, 1.0, 1.5, 2.5)
CMP_NS  <- c(100, 200, 400, 600)
CMP_R   <- 10
CMP_ITER <- 2000; CMP_WARMUP <- 1000; CMP_CHAINS <- 2

# ---------------------------------------------------------------------
# Helper functions (verbatim copies from the main script's Part 1,
# sections 1.1 / 1.2 / 1.2b -- fit_one_rep_part4 below needs all of them)
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

ok_result <- function(x) !is.null(x) && !inherits(x, "try-error")

nu_regime_label <- function(nu) {
  if (nu < 0.5) return(sprintf("Strong overdispersion regime ($\\nu=%.1f$)", nu))
  if (nu < 1.0) return(sprintf("Overdispersion regime ($\\nu=%.1f$)", nu))
  if (nu == 1.0) return("Equidispersion regime ($\\nu=1.0$)")
  if (nu <= 1.5) return(sprintf("Underdispersion regime ($\\nu=%.1f$)", nu))
  return(sprintf("Strong underdispersion regime ($\\nu=%.1f$)", nu))
}

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; AQUA <- "#1baf7a"; YELLOW <- "#eda100"

plot_winrate <- function(win_rate_df) {
  models_wr  <- c("ZICMP", "ZIP", "ZINB", "ZIGP")
  model_colors <- c(BLUE, ORANGE, AQUA, YELLOW)
  nu_num <- if (!is.null(win_rate_df$nu)) as.numeric(win_rate_df$nu) else as.numeric(sub("nu=", "", win_rate_df$regime))
  win_rate_df <- win_rate_df[order(nu_num), ]
  regimes_wr <- if (!is.null(win_rate_df$nu)) sprintf("nu=%s", win_rate_df$nu) else win_rate_df$regime
  winrate <- as.matrix(win_rate_df[, models_wr])
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

model_order <- c("ZICMP", "ZIGP", "ZINB", "ZIP")
crit_cols   <- c("EAIC", "EBIC", "DIC", "WAIC1", "WAIC2")

emit_compare_block <- function(sub, label) {
  sub <- sub[match(model_order, sub$model), ]
  best <- sapply(crit_cols, function(cc) which.min(sub[[cc]]))
  lines <- paste0("\\multirow{4}{*}{", label, "}")
  rows <- sapply(seq_len(nrow(sub)), function(i) {
    vals <- sapply(seq_along(crit_cols), function(j) {
      v <- sprintf("%.2f", sub[[crit_cols[j]]][i])
      if (best[j] == i) paste0("\\textbf{", v, "}") else v
    })
    paste0(" & ", sub$model[i], " & ", paste(vals, collapse = " & "), " \\\\")
  })
  c(paste0(lines, rows[1]), rows[-1], "\\midrule")
}

emit_winrate_row <- function(wr_row, label) {
  vals <- sapply(model_order, function(m) sprintf("%.2f", wr_row[[m]]))
  paste0(label, " & ", paste(vals, collapse = " & "), " \\\\")
}

# ---------------------------------------------------------------------
# Recompile mod -- REQUIRED every fresh R session (a saved stanmodel's
# compiled DSO pointer doesn't survive a restart). This is a real
# compile step, expect it to take a minute or two, not be instant.
# ---------------------------------------------------------------------
cat("Compiling", STAN_FILE, "-- this takes a minute or two, not instant...\n")
mod <- stan_model(STAN_FILE)
cat("Done compiling.\n")

# ---------------------------------------------------------------------
# Reload the REAL cmp_cache (all 20 cells were independently verified
# complete in the file you uploaded earlier) -- NOT the empty
# placeholder that a bare `if (!exists("cmp_cache")) cmp_cache <- list()`
# would silently create.
# ---------------------------------------------------------------------
cmp_cache <- readRDS("cmp_cache.rds")
cat("Loaded cmp_cache:", length(cmp_cache), "cells (expect 20).\n")
if (length(cmp_cache) != 20) {
  warning("cmp_cache has ", length(cmp_cache), " cells, not the expected 20 -- ",
          "double-check you're running this from the folder with the real cmp_cache.rds.")
}

# ---------------------------------------------------------------------
# Reload compare_results and defensively drop any cell that doesn't
# look complete (exactly CMP_R * 4 rows = 10 reps x 4 models). This
# covers the 4 known-broken cells AND the one that was mid-run
# (nu=1.5_n=600) when the earlier session was interrupted, without
# needing to hardcode which keys are bad.
# ---------------------------------------------------------------------
compare_results <- if (file.exists("compare_results_partial.rds")) {
  readRDS("compare_results_partial.rds")
} else list()

expected_rows <- CMP_R * 4
incomplete <- names(compare_results)[sapply(compare_results, function(df) is.null(df) || nrow(df) != expected_rows)]
if (length(incomplete) > 0) {
  cat("Dropping incomplete/broken cells before resuming:", paste(incomplete, collapse = ", "), "\n")
  compare_results[incomplete] <- NULL
}
cat("Cells already complete and kept:", length(compare_results), "of 20 --",
    paste(names(compare_results), collapse = ", "), "\n")

cmp_lambdas <- setNames(sapply(CMP_NUS, solve_lambda_for_mean, target_mu = TARGET_MU),
                         as.character(CMP_NUS))

# ---------------------------------------------------------------------
# Same per-rep worker as the main script's Part 4, section 4.3
# ---------------------------------------------------------------------
fit_one_rep_part4 <- function(r, nu_true, n_val, lam_true, key) {
  cached_cell <- cmp_cache[[key]]
  cached_rep <- if (!is.null(cached_cell) && r <= length(cached_cell)) cached_cell[[r]] else NULL
  reuse_cmp <- !is.null(cached_rep) && !is.null(cached_rep$crit_cmp)

  if (reuse_cmp) {
    y <- cached_rep$y
    crit_cmp <- cached_rep$crit_cmp
  } else {
    y <- simul_zinarCMP(n_val, ALPHA_TRUE, lam_true, nu_true, RHO_TRUE)
    fit_cmp <- tryCatch(
      sampling(mod, data = list(T = n_val, y = y, M = M,
                                 hybrid_tol = HYBRID_TOL, ff = 0),
               chains = CMP_CHAINS, iter = CMP_ITER, warmup = CMP_WARMUP,
               seed = 2000 + r, refresh = 0, control = list(adapt_delta = 0.95)),
      error = function(e) { message("  [rep ", r, "] ZICMP fit failed: ", conditionMessage(e)); NULL })
    if (is.null(fit_cmp)) return(NULL)
    crit_cmp <- get_mod_sel_cmp(y, fit_cmp, M, HYBRID_TOL)
  }

  fit_poi <- tryCatch(
    ZIHINAR1::get_stanfit(mod_type = "zi", distri = "poi", y = y,
                           n_pred = 0, chains = CMP_CHAINS, iter = CMP_ITER,
                           warmup = CMP_WARMUP, seed = 2000 + r),
    error = function(e) { message("  [rep ", r, "] ZIP fit failed: ", conditionMessage(e)); NULL })

  fit_nb <- tryCatch(
    ZIHINAR1::get_stanfit(mod_type = "zi", distri = "nb", y = y,
                           n_pred = 0, chains = CMP_CHAINS, iter = CMP_ITER,
                           warmup = CMP_WARMUP, seed = 2000 + r),
    error = function(e) { message("  [rep ", r, "] ZINB fit failed: ", conditionMessage(e)); NULL })

  fit_gp <- tryCatch(
    ZIHINAR1::get_stanfit(mod_type = "zi", distri = "gp", y = y,
                           n_pred = 0, chains = CMP_CHAINS, iter = CMP_ITER,
                           warmup = CMP_WARMUP, seed = 2000 + r),
    error = function(e) { message("  [rep ", r, "] ZIGP fit failed: ", conditionMessage(e)); NULL })

  if (is.null(fit_poi) || is.null(fit_nb) || is.null(fit_gp)) return(NULL)

  crit_poi <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "poi", stan_fit = fit_poi)
  crit_nb  <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "nb",  stan_fit = fit_nb)
  crit_gp  <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "gp",  stan_fit = fit_gp)

  row <- rbind(
    cbind(model = "ZICMP", rep = r, crit_cmp),
    cbind(model = "ZIP",   rep = r, crit_poi),
    cbind(model = "ZINB",  rep = r, crit_nb),
    cbind(model = "ZIGP",  rep = r, crit_gp)
  )
  cat(sprintf("  rep %d/%d  [%s]  EAIC[cmp,poi,nb,gp] = %.1f, %.1f, %.1f, %.1f  |  WAIC2[cmp,poi,nb,gp] = %.1f, %.1f, %.1f, %.1f\n",
              r, CMP_R, if (reuse_cmp) "CMP reused from Part 1" else "CMP fit fresh",
              crit_cmp$EAIC, crit_poi$EAIC, crit_nb$EAIC, crit_gp$EAIC,
              crit_cmp$WAIC2, crit_poi$WAIC2, crit_nb$WAIC2, crit_gp$WAIC2))

  list(row = row, reused = reuse_cmp)
}

# ---------------------------------------------------------------------
# Resume the grid
# ---------------------------------------------------------------------
t0 <- Sys.time()
total_reused <- 0L; total_fresh <- 0L

for (nu_true in CMP_NUS) {
  lam_true <- cmp_lambdas[[as.character(nu_true)]]
  for (n_val in CMP_NS) {
    key <- paste0("nu=", nu_true, "_n=", n_val)
    if (key %in% names(compare_results)) {
      cat("\n===", key, "=== already done, skipping\n")
      next
    }
    cat("\n===", key, "=== (", N_WORKERS, "worker(s) in parallel; output below may interleave)\n")

    rep_results <- mclapply(seq_len(CMP_R), fit_one_rep_part4, nu_true = nu_true, n_val = n_val,
                             lam_true = lam_true, key = key, mc.cores = N_WORKERS, mc.preschedule = FALSE)

    rows <- lapply(rep_results, function(x) if (ok_result(x)) x$row else NULL)
    cell_reused <- sum(sapply(rep_results, function(x) ok_result(x) && isTRUE(x$reused)))
    cell_fresh  <- sum(sapply(rep_results, function(x) ok_result(x) && !isTRUE(x$reused)))
    total_reused <- total_reused + cell_reused
    total_fresh  <- total_fresh + cell_fresh

    combined <- do.call(rbind, rows[!sapply(rows, is.null)])
    if (!is.null(combined) && nrow(combined) == expected_rows) {
      compare_results[[key]] <- combined
    } else {
      cat("  WARNING:", key, "did not complete cleanly (",
          if (is.null(combined)) 0 else nrow(combined), "/", expected_rows,
          "rows) -- left OUT of compare_results so it will be retried next run.\n")
    }
    saveRDS(compare_results, "compare_results_partial.rds")
    cat(sprintf("  elapsed so far: %.1f min  (this cell: %d reused, %d fresh CMP)\n",
                as.numeric(difftime(Sys.time(), t0, units = "mins")), cell_reused, cell_fresh))
  }
}
saveRDS(compare_results, "compare_results_final.rds")
cat(sprintf("\nCMP reused from Part 1's cache: %d reps total  |  CMP fit fresh here: %d reps total\n",
            total_reused, total_fresh))
cat("Cells complete:", length(compare_results), "of 20\n")

# ---------------------------------------------------------------------
# Only emit final summaries/LaTeX once ALL 20 cells are present -- if
# you're stopping partway through, just rerun this script later and it
# will pick up here once the grid is finally complete.
# ---------------------------------------------------------------------
if (length(compare_results) == 20) {

  compare_summary_tab <- do.call(rbind, lapply(names(compare_results), function(k) {
    df <- compare_results[[k]]
    agg <- aggregate(cbind(EAIC, EBIC, DIC, WAIC1, WAIC2) ~ model, df, mean)
    parts <- regmatches(k, regexec("^nu=([^_]+)_n=(.+)$", k))[[1]]
    cbind(regime = k, nu = parts[2], n = parts[3], agg)
  }))
  write.csv(compare_summary_tab, "compare_summary_table.csv", row.names = FALSE)

  win_rate <- do.call(rbind, lapply(names(compare_results), function(k) {
    df <- compare_results[[k]]
    reps <- unique(df$rep)
    winners <- sapply(reps, function(rr) {
      sub <- df[df$rep == rr, ]
      sub$model[which.min(sub$WAIC2)]
    })
    parts <- regmatches(k, regexec("^nu=([^_]+)_n=(.+)$", k))[[1]]
    wr_tab <- prop.table(table(factor(winners, levels = c("ZICMP", "ZIP", "ZINB", "ZIGP"))))
    data.frame(regime = k, nu = parts[2], n = parts[3],
               as.list(unclass(wr_tab)), check.names = FALSE)
  }))
  write.csv(win_rate, "compare_win_rate.csv", row.names = FALSE)

  cat("\n--- Mean model-selection criteria by regime (lower = better) ---\n")
  print(compare_summary_tab)
  cat("\n--- Fraction of replicates where each model has the LOWEST WAIC2 ---\n")
  print(win_rate)

  tab600 <- compare_summary_tab[compare_summary_tab$n == "600", ]
  wr600  <- win_rate[win_rate$n == "600", ]
  lines_tab <- unlist(lapply(sort(unique(as.numeric(tab600$nu))), function(nu) {
    emit_compare_block(tab600[tab600$nu == as.character(nu), ], nu_regime_label(nu))
  }))
  writeLines(lines_tab, "compare_table_n600_body.tex")
  lines_wr <- sapply(sort(unique(as.numeric(wr600$nu))), function(nu) {
    emit_winrate_row(wr600[wr600$nu == as.character(nu), ], nu_regime_label(nu))
  })
  writeLines(lines_wr, "compare_winrate_n600_body.tex")

  lines_tab_full <- unlist(lapply(names(compare_results), function(k) {
    parts <- regmatches(k, regexec("^nu=([^_]+)_n=(.+)$", k))[[1]]
    nu_val <- as.numeric(parts[2]); n_val <- parts[3]
    emit_compare_block(compare_summary_tab[compare_summary_tab$regime == k, ],
                        paste0(nu_regime_label(nu_val), ", $n=", n_val, "$"))
  }))
  writeLines(lines_tab_full, "compare_table_full_body.tex")
  lines_wr_full <- sapply(names(compare_results), function(k) {
    parts <- regmatches(k, regexec("^nu=([^_]+)_n=(.+)$", k))[[1]]
    nu_val <- as.numeric(parts[2]); n_val <- parts[3]
    emit_winrate_row(win_rate[win_rate$regime == k, ],
                      paste0(nu_regime_label(nu_val), ", $n=", n_val, "$"))
  })
  writeLines(lines_wr_full, "compare_winrate_full_body.tex")

  pdf("fig_modelcomparison.pdf", width = 7.5, height = 4.2); plot_winrate(wr600); dev.off()
  png("fig_modelcomparison.png", width = 1900, height = 1050, res = 220); plot_winrate(wr600); dev.off()

  cat("\nPART 4 FULLY COMPLETE. Wrote compare_summary_table.csv, compare_win_rate.csv,",
      "compare_table_n600_body.tex, compare_winrate_n600_body.tex,",
      "compare_table_full_body.tex, compare_winrate_full_body.tex,",
      "fig_modelcomparison.{pdf,png}.\n")
} else {
  cat("\nNot all 20 cells done yet (", length(compare_results),
      "/20) -- rerun `Rscript part4_resume.R` to continue; summaries/LaTeX/figure",
      "are written automatically once the grid is complete.\n")
}
