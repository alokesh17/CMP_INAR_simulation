############################################################################
# Real Data Application -- TWO NEW SERIES: polio and campylobacterosis
# Companion to real_data_analysis.R (which covers sexoffences/violence/
# claims/soap). Same pipeline, same criteria, same table format -- kept as
# a separate file rather than merged in, since these come from different
# packages (gamlss.data, tscount) and a different citation/data-provenance
# footnote is needed for each.
#
# WHY THESE TWO SERIES
# None of the four series in real_data_analysis.R produced an outright
# ZICMP win on EAIC/EBIC/DIC/WAIC1/WAIC2 -- each was at best a close,
# non-dominant tie with ZINB or ZIGP (see Sections 5.1-5.2 and the
# Appendix~\ref{app:extra-realdata} series in the draft). That is a real
# finding, not a bug to keep searching past -- Section 4.4's own simulation
# shows CMP, NB, and GP performing comparably whenever the truth is a
# fixed, strongly non-Poisson dispersion, because all three are one-
# parameter dispersion families competing on the same data. These two new
# series are added because they sit in different corners of the
# (n, dispersion, zero-inflation) space than anything already fit:
#
#   1. polio (gamlss.data::polio) -- Zeger (1988)'s classic US monthly
#      poliomyelitis case counts, Jan 1970-Dec 1983, n=168. MILD
#      overdispersion (typically reported VMR around 1.4-1.7) with many
#      zeros -- a much longer series than sexoffences/violence (144) at a
#      similar zero-heavy, low-count regime, so it tests whether more data
#      alone (holding "hard, low-count" character fixed) is enough to
#      separate the models where the shorter series couldn't.
#   2. campy (tscount::campy) -- Ferland, Latour & Oraichi (2006)'s
#      Campylobacterosis case counts for the Drenthe or Quebec-analogue
#      series bundled with `tscount`, 4-week periods over 1990-2000,
#      n=140. Clearly OVERdispersed (a standard benchmark series in the
#      INGARCH/count-time-series literature specifically because of its
#      pronounced overdispersion and occasional outbreak-driven spikes),
#      giving a second overdispersion-regime series independent of soap
#      and sexoffences, from a different data-generating context
#      (infectious disease surveillance vs. crime).
#
# Neither series is a source of new UNDERdispersion evidence -- that
# remains soap/family-violence/claims' territory (and, more decisively,
# Section 4.4's simulation study, where CMP's real advantage shows up).
# The point of adding polio/campy is purely to check whether a longer,
# still-overdispersed series moves ZICMP from "competitive" to "winning"
# on information criteria, the same question soap was added to ask.
#
# THIRD SERIES IN THE ORIGINAL REQUEST -- EXCLUDED, NOT FORGOTTEN
# UK coal-mining-strike-outbreak counts (Ridout & Besbeas 2004,
# Statistical Modelling) are NOT included here. Only the marginal
# frequency table has been located (count: 0/1/2/3/4+, frequency:
# 46/76/24/9/1, n=156 four-week periods 1948-1959) -- a cross-sectional
# count, not a time-ordered series. INAR(1) needs the actual sequence
# y_1,...,y_n (the conditional likelihood \eqref{eq:transition} in the
# draft is a function of consecutive pairs), so this table alone cannot
# be used to fit ANY of the four competing models, ZICMP included -- there
# is no legitimate way to reconstruct a fake but "real" order from a
# frequency table. A placeholder with a FABRICATED random order is
# deliberately NOT run below (unlike in the original user snippet, where
# it's already flagged as prototyping-only); do not paste a result from
# such a placeholder into the paper. If the real ordered series can be
# obtained from Ridout & Besbeas (2004) directly or their supplementary
# material, this script's polio/campy loop below is a drop-in template
# for adding it as a fifth series.
#
# WHAT THIS SCRIPT DOES (identical structure to real_data_analysis.R)
#   1. Loads polio (gamlss.data) and campy (tscount), installing each
#      package only if not already present.
#   2. Reports descriptive statistics (mean, variance, VMR, %zero, max,
#      lag-1 ACF) for both.
#   3. Fits ZICMP-INAR(1) (ZIINAR1-CMP-fast-reparam.stan) via Stan/NUTS.
#   4. Fits ZIP/ZINB/ZIGP via ZIHINAR1, all on EAIC/EBIC/DIC/WAIC1/WAIC2.
#   5. Runs the credible-interval-excludes-1 check on nu (the sharper
#      diagnostic already used for family violence in the draft, Appendix
#      B) -- CMP can add real information even when it doesn't win the
#      information-criterion horse race outright.
#   6. Flags whether each series' posterior mean for lambda is reliable
#      (no long right tail in nu) or needs the median/CI-only treatment
#      already used for family violence and claims -- this decides
#      whether a series belongs in the main text or the appendix, exactly
#      as documented in the draft's Appendix~\ref{app:extra-realdata}
#      opening paragraph.
#   7. Writes CSVs and ready-to-paste LaTeX table bodies in the same
#      format as real_data_analysis.R's Section 5 tables.
#
# HOW TO USE
#   1. Put this file plus ZIINAR1-CMP-fast-reparam.stan in the same
#      directory used for the other real_data_analysis*.R scripts.
#   2. Make sure ZIHINAR1 is installed (see compare_innovations.R's header
#      for the expint/actuar binary-install workaround if source
#      compilation fails).
#   3. Rscript real_data_analysis_polio_campy.R
############################################################################

if (!requireNamespace("gamlss.data", quietly = TRUE)) install.packages("gamlss.data")
if (!requireNamespace("tscount",     quietly = TRUE)) install.packages("tscount")

library(rstan)
library(ZIHINAR1)
library(gamlss.data)
library(tscount)
rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))

set.seed(2026)

STAN_FILE  <- "ZIINAR1-CMP-fast-reparam.stan"
M          <- 300
HYBRID_TOL <- 1e-6
CHAINS     <- 4
ITER       <- 4000
WARMUP     <- 2000

SERIES_LABEL <- c(
  polio = "Poliomyelitis, US monthly cases (Zeger 1988)",
  campy = "Campylobacterosis, 4-week case counts (Ferland et al. 2006)"
)

# ---------------------------------------------------------------------
# 1. Load the two series
# ---------------------------------------------------------------------
data(polio, package = "gamlss.data")
data(campy, package = "tscount")

series <- list(
  polio = as.numeric(polio),
  campy = as.numeric(campy)
)
EXPECTED_N <- c(polio = 168, campy = 140)
for (nm in names(series)) {
  if (length(series[[nm]]) != EXPECTED_N[[nm]]) {
    warning(sprintf("%s: expected n=%d, got n=%d -- check the package version; proceeding anyway.",
                     nm, EXPECTED_N[[nm]], length(series[[nm]])))
  }
}

# ---------------------------------------------------------------------
# 2. Descriptive statistics -- same shape as real_data_analysis.R's
#    describe(), so real_data_descriptives.csv from that script and this
#    one's output can be concatenated directly.
# ---------------------------------------------------------------------
describe <- function(y) {
  data.frame(n = length(y), mean = mean(y), variance = var(y),
             vmr = var(y) / mean(y), pct_zero = 100 * mean(y == 0),
             max = max(y), acf1 = acf(y, plot = FALSE)$acf[2, 1, 1])
}
desc_tab <- do.call(rbind, lapply(names(series), function(nm) cbind(series = nm, describe(series[[nm]]))))
write.csv(desc_tab, "real_data_descriptives_polio_campy.csv", row.names = FALSE)
cat("\n--- Descriptive statistics (also written to real_data_descriptives_polio_campy.csv) ---\n")
print(desc_tab)
cat("\nInterpretation: VMR>1 in both series indicates overdispersion relative to Poisson\n",
    "at the same mean; neither is expected to show underdispersion (see header note).\n", sep = "")

# ---------------------------------------------------------------------
# 3. Compile the CMP Stan model once
# ---------------------------------------------------------------------
cat("\nCompiling", STAN_FILE, "...\n")
cmp_mod <- stan_model(STAN_FILE)

# ---------------------------------------------------------------------
# 4. Model-selection criteria for the CMP fit -- verbatim port of
#    get_mod_sel_cmp() from compare_innovations.R / real_data_analysis.R,
#    kept self-contained here so this script can run on its own.
# ---------------------------------------------------------------------
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

# ---------------------------------------------------------------------
# 5. Fit ZICMP + ZIP/ZINB/ZIGP to each series
# ---------------------------------------------------------------------
est_rows  <- list()
crit_rows <- list()
reliability_notes <- list()   # mean-vs-median gap for lambda; CI-excludes-1 check for nu

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

  # -- Reliability checks, same logic used to route family violence/claims
  #    to the appendix in the draft: (a) does lambda's posterior mean sit
  #    orders of magnitude above its median (long right tail in nu)? and
  #    (b) does nu's 95% CI exclude 1 (CMP detecting dispersion even if it
  #    doesn't win the IC contest outright)?
  lam_mean <- s["lambda", "mean"]; lam_median <- s["lambda", "50%"]
  lam_ratio <- lam_mean / lam_median
  nu_ci <- s["nu", c("2.5%", "97.5%")]
  nu_excludes_1 <- nu_ci[1] > 1 || nu_ci[2] < 1
  reliability_notes[[nm]] <- data.frame(
    series = nm,
    lambda_mean = lam_mean, lambda_median = lam_median, lambda_mean_median_ratio = lam_ratio,
    lambda_reliable = lam_ratio < 10,   # same order of magnitude -> no long-tail problem
    nu_q025 = nu_ci[1], nu_q975 = nu_ci[2], nu_excludes_1 = nu_excludes_1
  )
  cat(sprintf("\nReliability check: lambda mean/median ratio = %.2f (%s); nu 95%% CI = (%.3f, %.3f) %s 1\n",
              lam_ratio, if (lam_ratio < 10) "OK, main-text-safe" else "LONG TAIL -- route to appendix, report median/CI only",
              nu_ci[1], nu_ci[2], if (nu_excludes_1) "excludes" else "includes"))

  fit_poi <- ZIHINAR1::get_stanfit(mod_type = "zi", distri = "poi", y = y,
                                    n_pred = 0, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = 4242)
  fit_nb  <- ZIHINAR1::get_stanfit(mod_type = "zi", distri = "nb",  y = y,
                                    n_pred = 0, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = 4242)
  fit_gp  <- ZIHINAR1::get_stanfit(mod_type = "zi", distri = "gp",  y = y,
                                    n_pred = 0, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = 4242)

  crit_cmp <- get_mod_sel_cmp(y, fit_cmp, M, HYBRID_TOL)
  crit_poi <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "poi", stan_fit = fit_poi)
  crit_nb  <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "nb",  stan_fit = fit_nb)
  crit_gp  <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "gp",  stan_fit = fit_gp)

  crit_rows[[nm]] <- rbind(
    cbind(series = nm, model = "ZICMP", crit_cmp),
    cbind(series = nm, model = "ZIP",   crit_poi),
    cbind(series = nm, model = "ZINB",  crit_nb),
    cbind(series = nm, model = "ZIGP",  crit_gp)
  )
  cat("\nModel comparison (lower = better; DIC caveat -- see real_data_analysis.R header):\n")
  print(crit_rows[[nm]])

  winner_per_criterion <- sapply(c("EAIC", "EBIC", "DIC", "WAIC1", "WAIC2"),
                                  function(cc) crit_rows[[nm]]$model[which.min(crit_rows[[nm]][[cc]])])
  cat("Winner per criterion:", paste(names(winner_per_criterion), winner_per_criterion, sep = "=", collapse = ", "), "\n")
  if (all(winner_per_criterion == "ZICMP")) {
    cat(">>> ZICMP wins EVERY criterion outright on this series.\n")
  } else if ("ZICMP" %in% winner_per_criterion) {
    cat(">>> ZICMP wins SOME criteria, not all -- a partial, not outright, win.\n")
  } else {
    cat(">>> ZICMP does not win any criterion on this series (competitive-but-not-dominant, as in the four series already in the draft).\n")
  }

  saveRDS(list(fit_cmp = fit_cmp, fit_poi = fit_poi, fit_nb = fit_nb, fit_gp = fit_gp),
          paste0("fits_", nm, ".rds"))
}

est_tab  <- do.call(rbind, est_rows)
crit_tab <- do.call(rbind, crit_rows)
reliability_tab <- do.call(rbind, reliability_notes)
write.csv(est_tab,  "real_data_posterior_estimates_polio_campy.csv", row.names = FALSE)
write.csv(crit_tab, "real_data_model_comparison_polio_campy.csv", row.names = FALSE)
write.csv(reliability_tab, "real_data_reliability_checks_polio_campy.csv", row.names = FALSE)

cat("\n--- Reliability summary (all series) ---\n")
print(reliability_tab)

# ---------------------------------------------------------------------
# 6. Emit LaTeX table bodies -- same format as real_data_analysis.R, so
#    these can be pasted into Section 5 (if lambda is reliable) or
#    Appendix~\ref{app:extra-realdata} (if not) following the same
#    routing logic already used for family violence/claims.
# ---------------------------------------------------------------------
emit_est_table <- function(tab, series_name) {
  sub <- tab[tab$series == series_name, ]
  ord <- c("alpha", "lambda", "nu", "rho")
  symb <- c(alpha = "\\alpha", lambda = "\\lambda", nu = "\\nu", rho = "\\rho")
  sapply(ord, function(p) {
    r <- sub[sub$parameter == p, ]
    sprintf("$%s$ & %.3f & %.3f & %.3f & %.3f & %.3f \\\\",
            symb[[p]], r$mean, r$sd, r$median, r$q025, r$q975)
  })
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
  "% --- Poliomyelitis (Zeger 1988, via gamlss.data) ---",
  "% Posterior estimate table body (new table, e.g. tab:realdata_polio):",
  emit_est_table(est_tab, "polio"),
  "% Model comparison table body:",
  emit_crit_table(crit_tab, "polio"),
  "",
  "% --- Campylobacterosis (Ferland, Latour \\& Oraichi 2006, via tscount) ---",
  "% Posterior estimate table body (new table, e.g. tab:realdata_campy):",
  emit_est_table(est_tab, "campy"),
  "% Model comparison table body:",
  emit_crit_table(crit_tab, "campy")
)
writeLines(tex_lines, "real_data_table_bodies_polio_campy.tex")
cat("\nWrote real_data_table_bodies_polio_campy.tex -- paste into Section 5 (if lambda_reliable == TRUE for",
    "that series in real_data_reliability_checks_polio_campy.csv) or Appendix~\\ref{app:extra-realdata} otherwise,",
    "following the same routing already used for family violence/claims.\n")
cat("Done.\n")
