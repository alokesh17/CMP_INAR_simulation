############################################################################
# Real Data Application driver for the ZICMP-INAR(1) paper (Section 5)
#
# DATA SOURCE -- please read this before running.
# The draft previously claimed the two series were "monthly drug offenses"
# and "monthly sex offenses" for census tract 2206, available in an R
# package called `tsinteger`. That is not quite right, and this script uses
# the corrected, verified version:
#
#   - There is no CRAN package called `tsinteger`. There IS a GitHub-only
#     package whose internal Package: field really is "tsinteger"
#     (Bourguignon, Santos & Patriota), at
#     https://github.com/projecttsinteger/tsintegerpackage -- never
#     submitted to CRAN.
#   - That package ships two real Pittsburgh crime series (144 monthly
#     counts, Jan 1990-Dec 2001), both traceable to the crime-data archive
#     at forecastingprinciples.com (the data underlying the 2005 NIJ-funded
#     Pittsburgh/Rochester crime-forecasting project, ICPSR study 4545):
#       * sexoffences -- sex offences, 21st police car beat
#       * violence    -- family violence, 11th police car beat
#     There is NO drug-offenses series in this package. The geography is
#     police car beats, not census tracts.
#   - The old "sex offenses" placeholder numbers already in the draft
#     (mean 0.59, variance 1.03) are in fact the real numbers for
#     `sexoffences` -- confirmed below -- so only the surrounding text
#     needs correcting, not that table. The old "drug offenses" numbers
#     (mean 2.11, variance 12.91) do not correspond to anything in this
#     package and have been dropped; `violence` (family violence, 11th car
#     beat) is used as the second series instead.
#
# WHAT THIS SCRIPT DOES
#   1. Downloads the two .rda series directly from the package's GitHub repo
#      (no CRAN dependency needed just to get the data).
#   2. Reports real descriptive statistics for both series.
#   3. Fits the ZICMP-INAR(1) model (ZIINAR1-CMP-fast-reparam.stan) to each
#      series via Stan/NUTS, exactly as in run_simulation_grid.R.
#   4. Fits ZIP/ZINB/ZIGP to each series via the ZIHINAR1 package, exactly
#      as in compare_innovations.R, for the "external point of comparison"
#      already promised in the draft text.
#   5. Writes CSVs and ready-to-paste LaTeX table bodies for both the
#      posterior-estimate table (format of Table \ref{tab:realdata_placeholder})
#      and a model-comparison table (format of Table \ref{tab:model_comparison}).
#
# HOW TO USE
#   1. Put this file plus ZIINAR1-CMP-fast-reparam.stan in the same directory
#      (the same one used for run_simulation_grid.R / compare_innovations.R).
#   2. Make sure ZIHINAR1 is installed (see compare_innovations.R's header
#      for the expint/actuar binary-install workaround if source compilation
#      fails).
#   3. Rscript real_data_analysis.R
############################################################################

library(rstan)
library(ZIHINAR1)
rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))

set.seed(2026)

STAN_FILE  <- "ZIINAR1-CMP-fast-reparam.stan"
M          <- 300
HYBRID_TOL <- 1e-6
CHAINS     <- 4
ITER       <- 4000
WARMUP     <- 2000

DATA_URLS <- list(
  sexoffences = "https://raw.githubusercontent.com/projecttsinteger/tsintegerpackage/master/data/sexoffences.rda",
  violence    = "https://raw.githubusercontent.com/projecttsinteger/tsintegerpackage/master/data/violence.rda"
)
SERIES_LABEL <- c(
  sexoffences = "Sex offenses (21st police car beat)",
  violence    = "Family violence (11th police car beat)"
)

# ---------------------------------------------------------------------
# 1. Download and load the two series
# ---------------------------------------------------------------------
dir.create("real_data", showWarnings = FALSE)
series <- list()
for (nm in names(DATA_URLS)) {
  destfile <- file.path("real_data", paste0(nm, ".rda"))
  if (!file.exists(destfile)) download.file(DATA_URLS[[nm]], destfile, mode = "wb", quiet = TRUE)
  e <- new.env()
  load(destfile, envir = e)
  series[[nm]] <- as.numeric(get(nm, envir = e))
  stopifnot(length(series[[nm]]) == 144)
}

# ---------------------------------------------------------------------
# 2. Descriptive statistics
# ---------------------------------------------------------------------
describe <- function(y) {
  data.frame(n = length(y), mean = mean(y), variance = var(y),
             vmr = var(y) / mean(y), pct_zero = 100 * mean(y == 0),
             max = max(y), acf1 = acf(y, plot = FALSE)$acf[2, 1, 1])
}
desc_tab <- do.call(rbind, lapply(names(series), function(nm) cbind(series = nm, describe(series[[nm]]))))
write.csv(desc_tab, "real_data_descriptives.csv", row.names = FALSE)
cat("\n--- Descriptive statistics (also written to real_data_descriptives.csv) ---\n")
print(desc_tab)

# ---------------------------------------------------------------------
# 3. Compile the CMP Stan model once
# ---------------------------------------------------------------------
cat("\nCompiling", STAN_FILE, "...\n")
cmp_mod <- stan_model(STAN_FILE)

# ---------------------------------------------------------------------
# 4. Fit ZICMP + ZIP/ZINB/ZIGP to each series
# ---------------------------------------------------------------------
source_functions_from_compare_innovations <- function() {
  # Minimal re-implementation of get_mod_sel_cmp() from compare_innovations.R
  # so this script is self-contained; keep in sync if that file changes.
  logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
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
                                          ph["nu", "mean"], ph["rho", "mean"], M, hybrid_tol))
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
  list(get_mod_sel_cmp = get_mod_sel_cmp)
}
helpers <- source_functions_from_compare_innovations()

est_rows <- list()
crit_rows <- list()

for (nm in names(series)) {
  y <- series[[nm]]
  cat("\n===", SERIES_LABEL[[nm]], "(n =", length(y), ") ===\n")

  fit_cmp <- sampling(cmp_mod, data = list(T = length(y), y = y, M = M,
                                            hybrid_tol = HYBRID_TOL, ff = 0),
                       chains = CHAINS, iter = ITER, warmup = WARMUP,
                       seed = 4242, control = list(adapt_delta = 0.95, max_treedepth = 12))

  s <- summary(fit_cmp, pars = c("alpha", "lambda", "nu", "rho"),
               probs = c(0.025, 0.5, 0.975))$summary
  est_rows[[nm]] <- data.frame(
    series = nm, parameter = rownames(s),
    mean = s[, "mean"], sd = s[, "sd"], median = s[, "50%"],
    q025 = s[, "2.5%"], q975 = s[, "97.5%"], rhat = s[, "Rhat"]
  )
  cat("\nZICMP-INAR(1) posterior summary:\n"); print(round(s, 4))

  fit_poi <- ZIHINAR1::get_stanfit(mod_type = "zi", distri = "poi", y = y,
                                    n_pred = 0, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = 4242)
  fit_nb  <- ZIHINAR1::get_stanfit(mod_type = "zi", distri = "nb",  y = y,
                                    n_pred = 0, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = 4242)
  fit_gp  <- ZIHINAR1::get_stanfit(mod_type = "zi", distri = "gp",  y = y,
                                    n_pred = 0, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = 4242)

  crit_cmp <- helpers$get_mod_sel_cmp(y, fit_cmp, M, HYBRID_TOL)
  crit_poi <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "poi", stan_fit = fit_poi)
  crit_nb  <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "nb",  stan_fit = fit_nb)
  crit_gp  <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "gp",  stan_fit = fit_gp)

  crit_rows[[nm]] <- rbind(
    cbind(series = nm, model = "ZICMP", crit_cmp),
    cbind(series = nm, model = "ZIP",   crit_poi),
    cbind(series = nm, model = "ZINB",  crit_nb),
    cbind(series = nm, model = "ZIGP",  crit_gp)
  )
  cat("\nModel comparison:\n"); print(crit_rows[[nm]])

  saveRDS(list(fit_cmp = fit_cmp, fit_poi = fit_poi, fit_nb = fit_nb, fit_gp = fit_gp),
          file.path("real_data", paste0("fits_", nm, ".rds")))
}

est_tab  <- do.call(rbind, est_rows)
crit_tab <- do.call(rbind, crit_rows)
write.csv(est_tab,  "real_data_posterior_estimates.csv", row.names = FALSE)
write.csv(crit_tab, "real_data_model_comparison.csv", row.names = FALSE)

# ---------------------------------------------------------------------
# 5. Emit LaTeX table bodies ready to paste into Section 5
# ---------------------------------------------------------------------
emit_est_table <- function(tab, series_name) {
  sub <- tab[tab$series == series_name, ]
  ord <- c("alpha", "lambda", "nu", "rho")
  symb <- c(alpha = "\\alpha", lambda = "\\lambda", nu = "\\nu", rho = "\\rho")
  lines <- sapply(ord, function(p) {
    r <- sub[sub$parameter == p, ]
    sprintf("$%s$ & %.3f & %.3f & %.3f & %.3f & %.3f \\\\",
            symb[[p]], r$mean, r$sd, r$median, r$q025, r$q975)
  })
  lines
}

emit_crit_table <- function(tab, series_name) {
  sub <- tab[tab$series == series_name, ]
  best <- sapply(c("EAIC", "EBIC", "DIC", "WAIC1", "WAIC2"), function(cc) which.min(sub[[cc]]))
  fmtcell <- function(val, is_best) if (is_best) sprintf("\\textbf{%.2f}", val) else sprintf("%.2f", val)
  lines <- character(nrow(sub))
  for (i in seq_len(nrow(sub))) {
    lines[i] <- sprintf(" & %s & %s & %s & %s & %s & %s \\\\",
                         sub$model[i],
                         fmtcell(sub$EAIC[i],  best["EAIC"]  == i),
                         fmtcell(sub$EBIC[i],  best["EBIC"]  == i),
                         fmtcell(sub$DIC[i],   best["DIC"]   == i),
                         fmtcell(sub$WAIC1[i], best["WAIC1"] == i),
                         fmtcell(sub$WAIC2[i], best["WAIC2"] == i))
  }
  lines
}

tex_lines <- c(
  "% --- Sex offenses (21st police car beat) ---",
  "% Posterior estimate table body (Table \\ref{tab:realdata_sexoffences}):",
  emit_est_table(est_tab, "sexoffences"),
  "% Model comparison table body:",
  emit_crit_table(crit_tab, "sexoffences"),
  "",
  "% --- Family violence (11th police car beat) ---",
  "% Posterior estimate table body (Table \\ref{tab:realdata_violence}):",
  emit_est_table(est_tab, "violence"),
  "% Model comparison table body:",
  emit_crit_table(crit_tab, "violence")
)
writeLines(tex_lines, "real_data_table_bodies.tex")
cat("\nWrote real_data_table_bodies.tex -- paste these rows into Section 5 of the draft.\n")
cat("Done.\n")
