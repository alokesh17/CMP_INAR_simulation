############################################################################
# Out-of-sample forecasting evaluation for Section 3.5 / 5.5 of the
# ZICMP-INAR(1) paper -- addresses predicted-reviewer point #16:
#   "The paper doesn't appear to have a real forecasting experiment ...
#    train on first 80%; predict next 20%; compare log predictive score;
#    MAE/RMSE; interval coverage; zero-event prediction."
#
# DESIGN (matches the methodology text of Section 3.5 in cmp_inar_draft.tex
# exactly -- keep the two in sync if either changes):
#   For each series of length n, hold out the final ceiling(0.2*n)
#   observations. Fit each of the four models ONCE to the first
#   n_tr = n - ceiling(0.2*n) observations only. Then, for every held-out
#   time point y_t (t = n_tr+1, ..., n), form the ONE-STEP-AHEAD posterior
#   predictive distribution of y_t given the TRUE observed y_{t-1} -- never
#   a simulated/compounded previous value -- using the posterior draws from
#   the training-only fit:
#       pi_hat(y_t = k | y_{t-1}) ~= (1/S) sum_s pi(k | y_{t-1}, theta^(s))
#   which is exactly the recursion already used for the in-sample posterior
#   predictive checks of Section 3.4, conditioned on data the model did not
#   see during estimation. This is computed EXACTLY (no simulation noise)
#   by convolving, for every posterior draw, a Binomial(y_{t-1}, alpha)
#   thinning distribution with the model's own zero-inflated innovation
#   pmf, truncated at K_MAX (chosen well above the largest observed count),
#   and averaging the resulting discrete pmf across draws.
#
# METRICS (Table tab:forecast_eval / Table~13 in the standalone excerpt)
#   LPS      = mean_t log(pi_hat(y_t | y_{t-1}))                (higher better)
#   MAE/RMSE = built from the predictive mean sum_k k*pi_hat(k | y_{t-1})
#   Coverage = fraction of y_t inside the equal-tailed 95% predictive
#              interval read off the cumulative predictive pmf
#   BS0      = mean_t ( 1(y_t=0) - pi_hat(0 | y_{t-1}) )^2       (lower better)
#
# IMPORTANT -- READ BEFORE TRUSTING THE ZINB / ZIGP NUMBERS
# The Poisson transition pmf is exact by construction (no free parameterization
# choice). The negative-binomial and generalized-Poisson transition pmfs below
# use the parameterizations standard in the INAR literature:
#   NB : mean mu, size r  (R's dnbinom(mu=, size=); Var = mu + mu^2/r)
#   GP : Consul-Jain mean mu, shape xi in (-1,1) (Var = mu/(1-xi)^2; xi<0
#        gives restricted underdispersion, xi>0 overdispersion, xi=0 Poisson)
# These are believed to match ZIHINAR1's internal parameterization but have
# NOT been re-derived from that package's Stan source in this environment
# (no package/network access here -- this script must be run in yours); the
# exact NB/GP dispersion parameter is a double unknown -- both its NAME
# (already confirmed NOT "r"/fixed across installs) and its MEANING (e.g.
# dnbinom's size directly, or its reciprocal) -- so nb_candidates()/
# gp_candidates() enumerate every plausible (name, meaning) pair from the
# fit's own parameter list and resolve_by_selfcheck() picks whichever one
# actually reproduces Stan's own likelihood, rather than assuming one.
# Concretely: for the TRAINING fit already needed for forecasting, the
# hand-rolled transition log-likelihood is computed AVERAGED OVER POSTERIOR
# DRAWS and compared to Stan's own `ll` generated quantity (also a draw
# average) for that same fit -- deliberately avoiding a point-estimate-vs-
# average comparison, which is confounded by Jensen's inequality (see the
# comment above check_against_stan() for why). If no candidate agrees with
# Stan's likelihood within LL_CHECK_TOL nats, that model's forecast row is
# filled with NA instead of an unverified number -- check the console
# output (it prints which candidate, if any, was resolved) before pasting
# anything into the table.
#
# DATA
# All four series load exactly as in real_data_analysis.R (verified GitHub
# source, projecttsinteger/tsintegerpackage): sexoffences, violence, and
# soap are each a plain numeric vector; claims.rda is a 120x11 matrix (one
# column per claimant) and column "claims6" (VMR=0.836, the genuinely
# underdispersed one used for Table~tab:realdata_claims) is selected below,
# exactly as real_data_analysis.R does.
#
# HOW TO USE
#   Put this file next to ZIINAR1-CMP-fast-reparam.stan, then:
#     Rscript forecast_evaluation.R
# Requires: rstan, ZIHINAR1, matrixStats (same as real_data_analysis.R /
# compare_innovations.R).
############################################################################

library(rstan)
library(ZIHINAR1)
rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))
set.seed(2026)

STAN_FILE    <- "ZIINAR1-CMP-fast-reparam.stan"
M            <- 300
HYBRID_TOL   <- 1e-6
CHAINS       <- 4
ITER         <- 4000
WARMUP       <- 2000
LL_CHECK_TOL <- 5.0     # max |hand-rolled - Stan ll| (nats, summed over series) before a model is trusted
K_MAX        <- 60      # truncation for the predictive pmf; raise if any series has max(y) > ~K_MAX/3
S_USE        <- 500     # posterior draws used to average the predictive pmf (subsampled for speed)

# ---------------------------------------------------------------------
# 1. Data: all four series, loaded exactly as in real_data_analysis.R
#    (same GitHub repo, same claims6 column selection).
# ---------------------------------------------------------------------
dir.create("real_data", showWarnings = FALSE)

DATA_URLS <- list(
  sexoffences = "https://raw.githubusercontent.com/projecttsinteger/tsintegerpackage/master/data/sexoffences.rda",
  violence    = "https://raw.githubusercontent.com/projecttsinteger/tsintegerpackage/master/data/violence.rda",
  claims      = "https://raw.githubusercontent.com/projecttsinteger/tsintegerpackage/master/data/claims.rda",
  soap        = "https://raw.githubusercontent.com/projecttsinteger/tsintegerpackage/master/data/soap.rda"
)
EXPECTED_N <- c(sexoffences = 144, violence = 144, claims = 120, soap = 242)

series <- list()
for (nm in names(DATA_URLS)) {
  destfile <- file.path("real_data", paste0(nm, ".rda"))
  if (!file.exists(destfile)) download.file(DATA_URLS[[nm]], destfile, mode = "wb", quiet = TRUE)
  e <- new.env(); load(destfile, envir = e)
  obj <- get(ls(e)[1], envir = e)
  series[[nm]] <- if (nm == "claims") as.numeric(obj[, "claims6"]) else as.numeric(obj)
  stopifnot(length(series[[nm]]) == EXPECTED_N[[nm]])
}
# Sanity check against the descriptive stats already reported in the draft.
stopifnot(abs(mean(series$claims) - 0.917) < 0.01, abs(var(series$claims) - 0.766) < 0.01)
stopifnot(abs(mean(series$soap) - 5.442) < 0.01, abs(var(series$soap) - 15.401) < 0.01)

SERIES_LABEL <- c(sexoffences = "Sex offenses", violence = "Family violence",
                   claims = "Claims", soap = "Soap")
SERIES_ORDER <- c("sexoffences", "violence", "claims", "soap")

# ---------------------------------------------------------------------
# 2. Shared machinery: zero-inflated-innovation transition pmf, generic
#    over the choice of base innovation log-pmf (CMP / Poisson / NB / GP).
#    The convolution logic (thinning-binomial * ZI-innovation) is the exact
#    same structure as cmp_transition_loglik() in real_data_analysis.R /
#    compare_innovations.R, generalized to return a full pmf vector (for
#    the forecast evaluation) rather than only the log-density at observed y.
# ---------------------------------------------------------------------
logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

# innov_logpmf(k, ...): vectorized over k = 0:K_MAX, returns log f(k; params)
cmp_innov_logpmf <- function(k, lambda, nu, M = M, hybrid_tol = HYBRID_TOL) {
  lgam <- lgamma(k + 1)
  log_lambda <- log(lambda)
  test <- exp(-log_lambda / nu)
  log_Z <- if (test < hybrid_tol) {
    nu * exp(log_lambda / nu) - ((nu - 1) / (2 * nu)) * log_lambda -
      ((nu - 1) / 2) * log(2 * pi) - 0.5 * log(nu)
  } else {
    r <- 0:M; logsumexp(r * log_lambda - nu * lgamma(r + 1))
  }
  k * log_lambda - nu * lgam - log_Z
}
poisson_innov_logpmf <- function(k, lambda) dpois(k, lambda, log = TRUE)
nb_innov_logpmf      <- function(k, mu, size) dnbinom(k, size = size, mu = mu, log = TRUE)
gp_innov_logpmf <- function(k, mu, xi) {
  # Consul-Jain (1973) generalized Poisson, mean mu, shape xi in (-1,1);
  # theta := mu*(1-xi) so that E[Y] = theta/(1-xi) = mu. Var = mu/(1-xi)^2.
  theta <- mu * (1 - xi)
  log(theta) + (k - 1) * log(pmax(theta + k * xi, 1e-12)) - lgamma(k + 1) - theta - k * xi
}

# Zero-inflate a base innovation log-pmf vector (k = 0:K) and convolve with
# Binomial(y_prev, alpha) thinning to get the full transition pmf vector
# pi(k | y_prev) for k = 0:K.
zi_transition_pmf <- function(y_prev, alpha, rho, base_logpmf_0K) {
  K <- length(base_logpmf_0K) - 1
  zi_pmf <- exp(base_logpmf_0K) * (1 - rho)
  zi_pmf[1] <- zi_pmf[1] + rho              # k = 0 term gets the extra zero-inflation mass
  bin_pmf <- dbinom(0:min(y_prev, K), y_prev, alpha)
  out <- numeric(K + 1)
  for (j in seq_along(bin_pmf) - 1) {
    idx <- (j):K
    out[idx + 1] <- out[idx + 1] + bin_pmf[j + 1] * zi_pmf[idx - j + 1]
  }
  out
}

# Average the transition pmf over S_use posterior draws.
avg_predictive_pmf <- function(y_prev, draws, model, K = K_MAX, m_cmp = M, tol_cmp = HYBRID_TOL) {
  S <- nrow(draws)
  idx <- if (S > S_USE) sample.int(S, S_USE) else seq_len(S)
  acc <- numeric(K + 1)
  for (s in idx) {
    d <- draws[s, ]
    base_logpmf <- switch(model,
      ZICMP = cmp_innov_logpmf(0:K, d[["lambda"]], d[["nu"]], m_cmp, tol_cmp),
      ZIP   = poisson_innov_logpmf(0:K, d[["lambda"]]),
      ZINB  = nb_innov_logpmf(0:K, d[["mu"]], d[["size"]]),
      ZIGP  = gp_innov_logpmf(0:K, d[["mu"]], d[["xi"]]))
    acc <- acc + zi_transition_pmf(y_prev, d[["alpha"]], d[["rho"]], base_logpmf)
  }
  acc / length(idx)
}

# ---------------------------------------------------------------------
# 3. In-sample self-check: does the hand-rolled transition log-likelihood,
#    AVERAGED OVER POSTERIOR DRAWS, match Stan's own `ll` generated
#    quantity (also a posterior-draw average) for the same fit?
#
#    NOTE ON DESIGN: this deliberately compares two draw-averaged
#    quantities, not "hand-rolled loglik at the posterior mean" vs
#    "Stan's mean(ll)". Those are NOT the same thing whenever the
#    posterior isn't a point mass -- log p(y | E[theta]) != E[log p(y |
#    theta)] (Jensen's inequality), and the gap between them is exactly
#    what DIC's effective-parameter correction (pdic = 2*(logphat -
#    mean(ll))) is built to measure. A point-estimate-vs-average check
#    would flag essentially every correctly-implemented model as a
#    "mismatch" whenever pdic is non-trivial, which is not a useful
#    trustworthiness signal. Averaging the hand-rolled loglik over the
#    same posterior draws removes that confound: both sides then
#    estimate E_theta[loglik(theta) | training data], so a real gap
#    reflects a genuine formula/parameterization mismatch, not curvature
#    of the log-likelihood surface.
# ---------------------------------------------------------------------
transition_loglik_one_draw <- function(y, alpha, rho, model, mu = NULL, nu = NULL, disp = NULL,
                                        K = K_MAX, m_cmp = M, tol_cmp = HYBRID_TOL) {
  base_logpmf <- switch(model,
    ZICMP = cmp_innov_logpmf(0:K, mu, nu, m_cmp, tol_cmp),
    ZIP   = poisson_innov_logpmf(0:K, mu),
    ZINB  = nb_innov_logpmf(0:K, mu, disp),
    ZIGP  = gp_innov_logpmf(0:K, mu, disp))
  Tt <- length(y); ll <- 0
  for (t in 2:Tt) {
    pmf <- zi_transition_pmf(y[t - 1], alpha, rho, base_logpmf)
    pmf <- pmf / sum(pmf)   # correct for K_MAX truncation so this isn't confused with a formula error
    yt <- y[t]
    ll <- ll + log(if (yt <= K) max(pmf[yt + 1], 1e-300) else 1e-300)
  }
  ll
}

check_against_stan <- function(y, stan_fit, model, ext, S_check = 300) {
  draws <- ext$draws
  S <- nrow(draws)
  idx <- if (S > S_check) sample.int(S, S_check) else seq_len(S)
  ll_vals <- vapply(idx, function(s) {
    dr <- draws[s, ]
    switch(model,
      ZICMP = transition_loglik_one_draw(y, dr[["alpha"]], dr[["rho"]], model, mu = dr[["lambda"]], nu = dr[["nu"]]),
      ZIP   = transition_loglik_one_draw(y, dr[["alpha"]], dr[["rho"]], model, mu = dr[["lambda"]]),
      ZINB  = transition_loglik_one_draw(y, dr[["alpha"]], dr[["rho"]], model, mu = dr[["mu"]], disp = dr[["size"]]),
      ZIGP  = transition_loglik_one_draw(y, dr[["alpha"]], dr[["rho"]], model, mu = dr[["mu"]], disp = dr[["xi"]]))
  }, numeric(1))
  ll_hand_avg <- mean(ll_vals)
  ll_stan <- mean(rstan::extract(stan_fit, pars = "ll")[[1]])
  diff <- abs(ll_hand_avg - ll_stan)
  ok <- diff <= LL_CHECK_TOL
  cat(sprintf("  [%s] hand-rolled E[ll] (%d draws) = %.2f, Stan mean(ll) = %.2f, |diff| = %.2f -> %s\n",
              model, length(idx), ll_hand_avg, ll_stan, diff, if (ok) "OK" else "MISMATCH -- see header comment"))
  ok
}

# Pull posterior draws in the layout avg_predictive_pmf() /
# check_against_stan() expect, per model. alpha, rho, and the
# innovation-mean parameter "lambda" are confirmed present under those
# exact names in ZIHINAR1's zi/poi, zi/nb, zi/gp stanfits (the ZIP
# self-check ran without a missing-parameter error).
#
# The DISPERSION parameter for NB/GP is a double unknown, not just a
# naming question: even once we know WHICH parameter it is (e.g. "phi"),
# there are multiple standard conventions for WHAT it means -- dnbinom's
# "size" directly, or its reciprocal (a "1/size" heterogeneity index);
# for GP, the Consul-Jain shape xi in (-1,1) directly, or a variance-to-
# mean-ratio phi with xi = 1 - 1/sqrt(phi) (the more common convention in
# applied GP regression). Rather than guess a single (name, meaning) pair
# again, build every plausible candidate below and let check_against_stan()
# -- already verified immune to the Jensen's-gap confound -- pick whichever
# one actually reproduces Stan's own likelihood.
extract_cmp <- function(fit) {
  d <- as.data.frame(rstan::extract(fit, pars = c("alpha", "lambda", "nu", "rho")))
  list(draws = d)
}
extract_poi <- function(fit) {
  d <- as.data.frame(rstan::extract(fit, pars = c("alpha", "lambda", "rho")))
  list(draws = d)
}

base_pars_renamed <- function(fit, extra_nm) {
  d <- as.data.frame(rstan::extract(fit, pars = c("alpha", "lambda", extra_nm, "rho")))
  names(d)[names(d) == "lambda"] <- "mu"
  d
}
nb_candidates <- function(fit) {
  avail <- fit@sim$pars_oi
  cand_names <- intersect(c("phi", "r", "size", "kappa", "theta", "disp", "nb_size", "od"), avail)
  if (length(cand_names) == 0)
    stop("No NB dispersion parameter found among this fit's parameters: ", paste(avail, collapse = ", "))
  out <- list()
  for (nm in cand_names) {
    raw <- base_pars_renamed(fit, nm)
    d_direct <- raw; names(d_direct)[names(d_direct) == nm] <- "size"                              # nm IS dnbinom's size
    d_inverse <- raw; d_inverse[[nm]] <- 1 / pmax(d_inverse[[nm]], 1e-8); names(d_inverse)[names(d_inverse) == nm] <- "size"  # nm is 1/size
    out[[paste0(nm, ":direct")]]  <- d_direct
    out[[paste0(nm, ":inverse")]] <- d_inverse
  }
  out
}
gp_candidates <- function(fit) {
  avail <- fit@sim$pars_oi
  cand_names <- intersect(c("phi", "xi", "disp", "kappa", "theta", "gp_disp"), avail)
  if (length(cand_names) == 0)
    stop("No GP dispersion parameter found among this fit's parameters: ", paste(avail, collapse = ", "))
  out <- list()
  for (nm in cand_names) {
    raw <- base_pars_renamed(fit, nm)
    d_asis <- raw; names(d_asis)[names(d_asis) == nm] <- "xi"                                          # nm IS Consul-Jain xi in (-1,1)
    d_vmr  <- raw; d_vmr[[nm]] <- 1 - 1 / sqrt(pmax(d_vmr[[nm]], 1e-6)); names(d_vmr)[names(d_vmr) == nm] <- "xi"  # nm is Var/mean ratio, phi=1/(1-xi)^2
    out[[paste0(nm, ":asis")]] <- d_asis
    out[[paste0(nm, ":vmr")]]  <- d_vmr
  }
  out
}

# Try each candidate (name, meaning) pair in turn; keep the first one whose
# hand-rolled likelihood matches Stan's own (within LL_CHECK_TOL).
resolve_by_selfcheck <- function(y, fit, model, candidates) {
  for (cand_name in names(candidates)) {
    ext <- list(draws = candidates[[cand_name]])
    cat("  [", model, "] trying '", cand_name, "':\n", sep = "")
    if (check_against_stan(y, fit, model, ext)) {
      cat("  [", model, "] resolved dispersion parameterization: '", cand_name, "'\n", sep = "")
      return(list(ok = TRUE, ext = ext, chosen = cand_name))
    }
  }
  cat("  [", model, "] no candidate parameterization matched Stan's likelihood -- row will be NA.\n", sep = "")
  list(ok = FALSE, ext = list(draws = candidates[[1]]), chosen = NA)
}

# ---------------------------------------------------------------------
# 4. Main loop: fit on training data, self-check, forecast the holdout.
# ---------------------------------------------------------------------
cat("Compiling", STAN_FILE, "...\n")
cmp_mod <- stan_model(STAN_FILE)

metric_rows <- list()

for (nm in SERIES_ORDER) {
  y <- series[[nm]]
  n <- length(y)
  n_test <- ceiling(0.2 * n)
  n_tr <- n - n_test
  y_tr <- y[1:n_tr]
  test_idx <- (n_tr + 1):n
  cat("\n===", SERIES_LABEL[[nm]], "(n =", n, ", n_tr =", n_tr, ", n_test =", n_test, ") ===\n")

  fit_cmp <- sampling(cmp_mod, data = list(T = n_tr, y = y_tr, M = M, hybrid_tol = HYBRID_TOL, ff = 0),
                       chains = CHAINS, iter = ITER, warmup = WARMUP, seed = 4242,
                       control = list(adapt_delta = 0.95, max_treedepth = 12))
  fit_poi <- ZIHINAR1::get_stanfit(mod_type = "zi", distri = "poi", y = y_tr, n_pred = 0,
                                    chains = CHAINS, iter = ITER, warmup = WARMUP, seed = 4242)
  fit_nb  <- ZIHINAR1::get_stanfit(mod_type = "zi", distri = "nb",  y = y_tr, n_pred = 0,
                                    chains = CHAINS, iter = ITER, warmup = WARMUP, seed = 4242)
  fit_gp  <- ZIHINAR1::get_stanfit(mod_type = "zi", distri = "gp",  y = y_tr, n_pred = 0,
                                    chains = CHAINS, iter = ITER, warmup = WARMUP, seed = 4242)

  cat("  Self-check (posterior-draw-averaged hand-rolled transition loglik vs Stan mean(ll)):\n")
  ext_cmp <- extract_cmp(fit_cmp); ok_cmp <- check_against_stan(y_tr, fit_cmp, "ZICMP", ext_cmp)
  ext_poi <- extract_poi(fit_poi); ok_poi <- check_against_stan(y_tr, fit_poi, "ZIP",   ext_poi)
  res_nb  <- resolve_by_selfcheck(y_tr, fit_nb, "ZINB", nb_candidates(fit_nb)); ext_nb <- res_nb$ext; ok_nb <- res_nb$ok
  res_gp  <- resolve_by_selfcheck(y_tr, fit_gp, "ZIGP", gp_candidates(fit_gp)); ext_gp <- res_gp$ext; ok_gp <- res_gp$ok

  models <- list(ZICMP = list(ext = ext_cmp, ok = ok_cmp), ZIP = list(ext = ext_poi, ok = ok_poi),
                 ZINB  = list(ext = ext_nb,  ok = ok_nb),  ZIGP = list(ext = ext_gp,  ok = ok_gp))

  for (model in names(models)) {
    mres <- models[[model]]
    if (!mres$ok) {
      cat("  Skipping", model, "forecast (self-check failed) -- row will be NA.\n")
      metric_rows[[paste(nm, model)]] <- data.frame(series = nm, model = model,
                                                      LPS = NA, MAE = NA, RMSE = NA, Coverage = NA, BS0 = NA)
      next
    }
    logdens <- numeric(n_test); pred_mean <- numeric(n_test)
    covered <- logical(n_test); p_zero <- numeric(n_test)
    for (i in seq_along(test_idx)) {
      t <- test_idx[i]
      pmf <- avg_predictive_pmf(y[t - 1], mres$ext$draws, model)
      pmf <- pmax(pmf, 1e-300); pmf <- pmf / sum(pmf)
      yt <- y[t]
      logdens[i] <- log(if (yt <= K_MAX) pmf[yt + 1] else 1e-300)
      pred_mean[i] <- sum((0:K_MAX) * pmf)
      cdf <- cumsum(pmf)
      lo <- which(cdf >= 0.025)[1] - 1; hi <- which(cdf >= 0.975)[1] - 1
      covered[i] <- (yt >= lo) && (yt <= hi)
      p_zero[i] <- pmf[1]
    }
    metric_rows[[paste(nm, model)]] <- data.frame(
      series = nm, model = model,
      LPS = mean(logdens),
      MAE = mean(abs(y[test_idx] - pred_mean)),
      RMSE = sqrt(mean((y[test_idx] - pred_mean)^2)),
      Coverage = mean(covered),
      BS0 = mean((as.numeric(y[test_idx] == 0) - p_zero)^2))
  }

  saveRDS(list(fit_cmp = fit_cmp, fit_poi = fit_poi, fit_nb = fit_nb, fit_gp = fit_gp),
          file.path("real_data", paste0("forecast_fits_", nm, ".rds")))
}

forecast_tab <- do.call(rbind, metric_rows)
write.csv(forecast_tab, "forecast_evaluation_results.csv", row.names = FALSE)
cat("\n--- Held-out forecasting metrics (also written to forecast_evaluation_results.csv) ---\n")
print(forecast_tab, digits = 3)

# ---------------------------------------------------------------------
# 5. Emit the LaTeX table body, ready to paste into Table~tab:forecast_eval
#    (replacing the placeholder rows in cmp_inar_draft.tex and
#    cmp-inar-realdata-section.tex; also delete the %% Pending comment
#    above that table once real numbers are in).
# ---------------------------------------------------------------------
n_test_by_series <- sapply(series[SERIES_ORDER], function(y) ceiling(0.2 * length(y)))
fmt <- function(x, d = 3) if (is.na(x)) "--" else formatC(x, digits = d, format = "f")

best_flags <- function(sub) {
  # boldface the best value per column within a series block
  list(LPS = sub$LPS == max(sub$LPS, na.rm = TRUE),
       MAE = sub$MAE == min(sub$MAE, na.rm = TRUE),
       RMSE = sub$RMSE == min(sub$RMSE, na.rm = TRUE),
       Coverage = abs(sub$Coverage - 0.95) == min(abs(sub$Coverage - 0.95), na.rm = TRUE),
       BS0 = sub$BS0 == min(sub$BS0, na.rm = TRUE))
}
bf <- function(val_str, is_best) if (isTRUE(is_best)) paste0("\\textbf{", val_str, "}") else val_str

lines <- character(0)
for (nm in SERIES_ORDER) {
  sub <- forecast_tab[forecast_tab$series == nm, ]
  sub <- sub[match(c("ZICMP", "ZIP", "ZINB", "ZIGP"), sub$model), ]
  flags <- best_flags(sub)
  cap <- sprintf("\\multirow{4}{*}{%s ($n_{\\mathrm{test}}=%d$)}", SERIES_LABEL[[nm]], n_test_by_series[[nm]])
  if (nrow(sub) > 0) lines <- c(lines, cap)
  for (i in seq_len(nrow(sub))) {
    row <- sub[i, ]
    lines <- c(lines, sprintf(" & %s & %s & %s & %s & %s & %s \\\\",
      row$model,
      bf(fmt(row$LPS), flags$LPS[i]), bf(fmt(row$MAE), flags$MAE[i]),
      bf(fmt(row$RMSE), flags$RMSE[i]), bf(fmt(row$Coverage, 2), flags$Coverage[i]),
      bf(fmt(row$BS0), flags$BS0[i])))
  }
  lines <- c(lines, "\\midrule")
}
lines[length(lines)] <- "\\bottomrule"  # last block doesn't want a trailing midrule
writeLines(lines, "forecast_eval_table_body.tex")
cat("\nWrote forecast_eval_table_body.tex -- paste its contents into Table~tab:forecast_eval\n",
    "in place of the placeholder rows (both cmp_inar_draft.tex and\n",
    "cmp-inar-realdata-section.tex), then delete the '%% Pending' comment above it.\n")
