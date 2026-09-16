############################################################################
# Predictive analysis for the ZICMP-INAR(1) model: posterior summaries of
# theta = (alpha, lambda, nu, rho) plus a posterior-predictive forecast
# plot of observed vs. predicted counts, in the style requested:
#
#   qoi <- c("lambda1", "lambda2", "alpha", "rho", "bic", "aic")
#   print(fit_stanGP, pars = qoi)
#   fitPred  <- summary(fit_stanGP, pars = "y_pred", probs = c(0.1, 0.9))$summary
#   fitPredM <- fitPred[, 1]
#   fitPred05 <- fitPred[, 4]
#   fitPred95 <- fitPred[, 5]
#   plot(y1[(n - ff + 1):n], ylim = c(0, 20))
#   lines(fitPredM); lines(fitPred05); lines(fitPred95)
#
# adapted to the ZICMP-INAR(1) model's ACTUAL parameters and generated
# quantities, which are not "lambda1, lambda2" (that naming is from a
# different model) but alpha, lambda, nu, rho, with aic/bic/y_pred already
# defined in the `generated quantities` block of both
# ZIINAR1-CMP-fast.stan and ZIINAR1-CMP-fast-reparam.stan:
#
#   parameters:          alpha, lambda (or mu_cmp with lambda := mu_cmp^nu
#                        in the -reparam version), nu, rho
#   generated quantities: y_pred[1:(ff+1)]  -- y_pred[1] = y[T] (the last
#                            TRAINING observation, an anchor, not a
#                            forecast); y_pred[2:(ff+1)] are the actual
#                            ff-step-ahead chained forecasts
#                         log_lik[1:T], ll, aic, bic
#
# This script does NOT fit anything itself -- it assumes you already have
# a `stanfit` object (from rstan::sampling() on one of the two .stan files
# above, called with data$ff = ff > 0 so y_pred actually forecasts
# something) and shows how to pull the quantities of interest out of it
# and plot them. No plots are rendered by this file being sourced without
# a fit object in scope; wire in your own `fit` and `y` before running.
############################################################################

library(rstan)

# ---------------------------------------------------------------------
# 0. Quantities of interest for the print() summary.
#    "lambda1, lambda2" in the requested style don't exist in this model;
#    the ZICMP-INAR(1) posterior has exactly four structural parameters
#    (alpha, lambda, nu, rho) plus the two generated-quantities model-
#    selection numbers (aic, bic) computed inside the Stan file itself.
# ---------------------------------------------------------------------
qoi <- c("alpha", "lambda", "nu", "rho", "aic", "bic")

# Example (uncomment once `fit` is in scope):
# print(fit, pars = qoi)

# ---------------------------------------------------------------------
# 1. Core function: pull the posterior-predictive forecast summary out of
#    a fitted model and overlay it against the held-out observed counts.
#
#    Arguments:
#      fit    -- a stanfit object from sampling() on ZIINAR1-CMP-fast.stan
#                or ZIINAR1-CMP-fast-reparam.stan, fit with data$ff = ff
#      y      -- the FULL observed series (length n) that was split into
#                training (first n - ff points, used as data$y/data$T) and
#                the ff held-out points forecast by y_pred
#      ff     -- number of one-step-ahead-chained forecasts requested from
#                Stan (must match data$ff used when fitting)
#      probs  -- predictive-interval quantiles to extract; c(0.1, 0.9) to
#                match the requested style (an 80% band), or c(0.025,
#                0.975) for the 95% band used elsewhere in the paper
#      series_name, model_name -- only used for the plot title
#
#    Returns (invisibly) a data.frame with the observed, mean-forecast,
#    and lower/upper predictive-interval columns, so the same numbers can
#    feed a table as well as the plot.
# ---------------------------------------------------------------------
plot_predictive <- function(fit, y, ff, probs = c(0.1, 0.9),
                             series_name = "", model_name = "ZICMP",
                             ylim = NULL) {
  n <- length(y)
  stopifnot(ff > 0, ff < n)

  # y_pred has ff+1 entries: y_pred[1] = y[T] (last TRAINING point, an
  # anchor for the recursion, not itself a forecast); y_pred[2:(ff+1)]
  # are the ff chained one-step-ahead forecasts, corresponding to
  # y[(n-ff+1):n] in the full series.
  fitPred   <- summary(fit, pars = "y_pred", probs = probs)$summary
  fitPredM  <- fitPred[-1, "mean"]                       # drop the anchor row
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
       main = sprintf("%s%s one-step-ahead forecasts%s",
                       model_name, if (nzchar(series_name)) paste0(" -- ", series_name) else "",
                       sprintf(" (%d%% band)", round(100 * (probs[2] - probs[1])))))
  lines(fitPredM, col = "#2a78d6", lwd = 1.8)
  lines(fitPredLo, col = "#2a78d6", lty = 2)
  lines(fitPredHi, col = "#2a78d6", lty = 2)
  legend("topleft", legend = c("Observed", "Posterior predictive mean",
                                sprintf("%d%%-%d%% interval", 100 * probs[1], 100 * probs[2])),
         col = c("black", "#2a78d6", "#2a78d6"), lty = c(1, 1, 2), pch = c(16, NA, NA),
         bty = "n", cex = 0.85)

  invisible(data.frame(
    t = seq_len(ff), observed = y_obs,
    pred_mean = fitPredM, pred_lo = fitPredLo, pred_hi = fitPredHi
  ))
}

# ---------------------------------------------------------------------
# 2. Worked example, matching the training/held-out split already used in
#    forecast_evaluation.R (final 20% of each series held out). Fill in
#    the paths/objects for whichever series and model you're plotting.
# ---------------------------------------------------------------------
if (FALSE) {   # guard so sourcing this file does nothing until you edit it

  # -- Load one real series (example: sex offenses) --------------------
  y_full <- scan("sexoffences.txt")     # replace with your actual loader
  n      <- length(y_full)
  ff     <- round(0.20 * n)             # same 20% holdout as forecast_evaluation.R
  n_tr   <- n - ff

  y_train <- y_full[1:n_tr]

  # -- Fit (example: reparameterized Stan file) -------------------------
  mod <- stan_model("ZIINAR1-CMP-fast-reparam.stan")
  fit <- sampling(mod,
                   data = list(T = n_tr, y = y_train, M = 300,
                               hybrid_tol = 1e-6, ff = ff),
                   chains = 4, iter = 2000, warmup = 1000, seed = 1,
                   control = list(adapt_delta = 0.95))

  # -- Posterior summary of theta + AIC/BIC, in the requested style -----
  print(fit, pars = qoi)

  # -- Predictive plot, 80% band as in the requested style --------------
  pred_tab <- plot_predictive(fit, y_full, ff, probs = c(0.1, 0.9),
                               series_name = "sex offenses", model_name = "ZICMP")

  # -- Or the 95% band used elsewhere in the paper's tables --------------
  pred_tab95 <- plot_predictive(fit, y_full, ff, probs = c(0.025, 0.975),
                                 series_name = "sex offenses", model_name = "ZICMP")

  # -- Repeat for the other three series by changing y_full/ff/fit; a ---
  #    2x2 layout puts all four series' forecast plots on one page:
  # par(mfrow = c(2, 2))
  # plot_predictive(fit_sexoffences, y_sexoffences, ff_sexoffences, series_name = "Sex offenses")
  # plot_predictive(fit_violence,    y_violence,    ff_violence,    series_name = "Family violence")
  # plot_predictive(fit_claims,      y_claims,      ff_claims,      series_name = "Claims")
  # plot_predictive(fit_soap,        y_soap,        ff_soap,        series_name = "Soap")
}
