############################################################################
# R port of make_sim_figures.py, for readers who want to stay in R end to
# end (simulation -> figures) rather than switching to Python for plotting.
# Produces the same three figures as make_sim_figures.py:
#   fig_recovery.pdf/png        -- Section 4.2, Table 1 (parameter recovery)
#   fig_modelcomparison.pdf/png -- Section 4.4, Table 3 (win-rate by regime)
#   fig_bias.pdf/png            -- Section 4.2, estimation bias
# plus a fourth figure this version adds for the two EXTRA dispersion
# regimes (strong over-/underdispersion) requested alongside the original
# three:
#   fig_recovery_extended.pdf/png -- all five regimes (nu = 0.2, 0.5, 1.0,
#                                     1.5, 2.5) on the same panels, so the
#                                     new extremes can be read directly
#                                     against the original three.
#
# DATA SOURCE: by default this script reads the actual simulation output
# (sim_summary_table.csv from run_simulation_grid.R, and, if present,
# sim_summary_table_ext.csv from run_simulation_grid_extended.R) rather
# than transcribing numbers by hand, so it stays correct if you rerun the
# grid. If those CSVs are not found in the working directory, it falls
# back to the exact numbers already in the draft's Table 1 /
# Table~\ref{tab:model_comparison_winrate}, so the three original figures
# can still be reproduced without rerunning Stan.
#
# No simulation is run here -- this is purely the plotting step. As
# requested, the plots themselves are NOT rendered by this script call;
# running it will produce the PDF/PNG files, but that final decision
# (which numbers, which regimes to include/exclude) is left to you --
# edit EXTRA_NUS / the fallback tables below as needed before running.
#
# Requires only base R graphics (no ggplot2 dependency), so it runs
# anywhere run_simulation_grid.R already runs.
############################################################################

# ---------------------------------------------------------------------
# 0. Palette (matches make_sim_figures.py exactly)
# ---------------------------------------------------------------------
BLUE   <- "#2a78d6"
ORANGE <- "#eb6834"
AQUA   <- "#1baf7a"
YELLOW <- "#eda100"
GREY   <- "#52514e"

n_vals <- c(100, 200, 400, 600)

# ---------------------------------------------------------------------
# 1. Load data: prefer the real simulation-grid CSVs; fall back to the
#    numbers already transcribed into the draft if the CSVs aren't there.
# ---------------------------------------------------------------------
fallback_table <- list(
  "0.5" = list(color = BLUE, true_nu = 0.5, true_lambda = 1.56,
    alpha = list(mean = c(0.291, 0.312, 0.300, 0.301), sd = c(0.063, 0.043, 0.034, 0.022)),
    lambda = list(mean = c(2.435, 2.052, 1.943, 1.719), sd = c(0.803, 0.770, 0.480, 0.314)),
    nu = list(mean = c(0.698, 0.620, 0.605, 0.548), sd = c(0.205, 0.212, 0.148, 0.116)),
    rho = list(mean = c(0.313, 0.296, 0.314, 0.297), sd = c(0.075, 0.060, 0.063, 0.043))),
  "1" = list(color = ORANGE, true_nu = 1.0, true_lambda = 3.00,
    alpha = list(mean = c(0.323, 0.283, 0.295, 0.297), sd = c(0.067, 0.053, 0.036, 0.024)),
    lambda = list(mean = c(5.077, 4.314, 3.410, 3.480), sd = c(1.924, 1.474, 1.057, 1.114)),
    nu = list(mean = c(1.236, 1.150, 1.030, 1.056), sd = c(0.299, 0.240, 0.191, 0.208)),
    rho = list(mean = c(0.340, 0.307, 0.296, 0.293), sd = c(0.078, 0.059, 0.055, 0.039))),
  "1.5" = list(color = AQUA, true_nu = 1.5, true_lambda = 5.66,
    alpha = list(mean = c(0.313, 0.292, 0.307, 0.306), sd = c(0.070, 0.038, 0.032, 0.029)),
    lambda = list(mean = c(5.957, 6.670, 7.370, 6.430), sd = c(2.327, 2.315, 2.658, 1.607)),
    nu = list(mean = c(1.427, 1.514, 1.617, 1.547), sd = c(0.324, 0.210, 0.256, 0.187)),
    rho = list(mean = c(0.308, 0.284, 0.305, 0.297), sd = c(0.074, 0.049, 0.036, 0.032)))
)

load_regime_table <- function(csv_path, lambda_lookup = NULL) {
  # Turns a sim_summary_table[_ext].csv (as written by run_simulation_grid*.R)
  # into the same list-of-regimes structure as fallback_table, so the
  # plotting code below can treat real and fallback data identically.
  tab <- read.csv(csv_path, stringsAsFactors = FALSE)
  # cell column looks like "nu=0.2_n=100"
  tab$nu <- as.numeric(sub("nu=([^_]+)_n=.*", "\\1", tab$cell))
  tab$n  <- as.numeric(sub(".*_n=", "", tab$cell))
  nus <- sort(unique(tab$nu))
  out <- list()
  for (nu in nus) {
    sub_tab <- tab[tab$nu == nu, ]
    sub_tab <- sub_tab[match(n_vals, sub_tab$n), ]  # order by n_vals, NA if missing
    out[[as.character(nu)]] <- list(
      true_nu = nu,
      true_lambda = if (!is.null(lambda_lookup)) lambda_lookup[[as.character(nu)]] else NA,
      alpha  = list(mean = sub_tab$alpha_MCmean,  sd = sub_tab$alpha_MCsd),
      lambda = list(mean = sub_tab$lambda_MCmean, sd = sub_tab$lambda_MCsd),
      nu     = list(mean = sub_tab$nu_MCmean,     sd = sub_tab$nu_MCsd),
      rho    = list(mean = sub_tab$rho_MCmean,    sd = sub_tab$rho_MCsd)
    )
  }
  out
}

regime_colors <- c(BLUE, ORANGE, AQUA, YELLOW, "#a259d9")  # 5th color for a 5th regime

if (file.exists("sim_summary_table.csv")) {
  data3 <- load_regime_table("sim_summary_table.csv")
  for (nm in names(data3)) data3[[nm]]$color <- regime_colors[match(nm, sort(names(data3)))]
} else {
  message("sim_summary_table.csv not found -- using the numbers already transcribed into the draft.")
  data3 <- fallback_table
}

data_ext <- NULL
if (file.exists("sim_summary_table_ext.csv")) {
  data_ext <- load_regime_table("sim_summary_table_ext.csv")
}

param_labels <- c(alpha = "hat(alpha)", lambda = "hat(lambda)", nu = "hat(nu)", rho = "hat(rho)")
panel_order  <- c("nu", "alpha", "lambda", "rho")   # nu first: the parameter of interest

# ---------------------------------------------------------------------
# 2. Figure 1: parameter recovery (mean +/- SD vs n, one panel per parameter)
# ---------------------------------------------------------------------
plot_recovery <- function(data_list, file_stem, jitter_frac = 0.03, title_suffix = "") {
  regimes <- names(data_list)
  n_regimes <- length(regimes)
  jitters <- (seq_len(n_regimes) - (n_regimes + 1) / 2) * jitter_frac * diff(range(n_vals))

  make_plot <- function() {
    op <- par(mfrow = c(2, 2), mar = c(4, 4.5, 2.5, 1), oma = c(4.5, 0, 2, 0))
    on.exit(par(op))
    for (pname in panel_order) {
      all_means <- unlist(lapply(data_list, function(d) d[[pname]]$mean))
      all_sds   <- unlist(lapply(data_list, function(d) d[[pname]]$sd))
      ylim <- range(c(all_means - all_sds, all_means + all_sds), na.rm = TRUE)
      plot(NA, xlim = range(n_vals), ylim = ylim, xaxt = "n",
           xlab = "Sample size n", ylab = paste0(param_labels[pname], " (posterior mean +/- SD)"),
           main = param_labels[pname])
      axis(1, at = n_vals)
      for (i in seq_along(regimes)) {
        reg <- regimes[i]; d <- data_list[[reg]]
        x <- n_vals + jitters[i]
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

plot_recovery(data3, "fig_recovery")
if (!is.null(data_ext)) {
  data5 <- c(data3, data_ext)
  data5 <- data5[order(as.numeric(names(data5)))]
  plot_recovery(data5, "fig_recovery_extended", title_suffix = " (5 regimes: adds strong over/underdispersion)")
}

# ---------------------------------------------------------------------
# 3. Figure 2: estimation bias (posterior mean - true value)
# ---------------------------------------------------------------------
plot_bias <- function(data_list, file_stem) {
  regimes <- names(data_list)
  n_regimes <- length(regimes)
  jitters <- (seq_len(n_regimes) - (n_regimes + 1) / 2) * 0.03 * diff(range(n_vals))

  make_plot <- function() {
    op <- par(mfrow = c(2, 2), mar = c(4, 4.5, 2.5, 1), oma = c(4.5, 0, 2, 0))
    on.exit(par(op))
    for (pname in panel_order) {
      biases <- lapply(data_list, function(d) {
        true_val <- if (pname == "nu") d$true_nu else if (pname %in% c("alpha", "rho")) 0.30 else d$true_lambda
        d[[pname]]$mean - true_val
      })
      ylim <- range(unlist(biases), na.rm = TRUE)
      plot(NA, xlim = range(n_vals), ylim = ylim, xaxt = "n",
           xlab = "Sample size n", ylab = paste0("Bias: ", param_labels[pname], " - true value"),
           main = param_labels[pname])
      axis(1, at = n_vals)
      abline(h = 0, col = "grey40", lwd = 1)
      for (i in seq_along(regimes)) {
        reg <- regimes[i]; d <- data_list[[reg]]
        x <- n_vals + jitters[i]
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

plot_bias(data3, "fig_bias")
if (!is.null(data_ext)) {
  data5 <- c(data3, data_ext)
  data5 <- data5[order(as.numeric(names(data5)))]
  plot_bias(data5, "fig_bias_extended")
}

# ---------------------------------------------------------------------
# 4. Figure 3: model-comparison win rate (Table~\ref{tab:model_comparison_winrate})
#    Numbers transcribed from the draft (R=10 pilot); replace with your
#    own win-rate table if you rerun compare_innovations.R / compare_pois_nb.R
#    for the two new regimes.
# ---------------------------------------------------------------------
regimes_wr <- c("Overdispersion\n(nu=0.5)", "Equidispersion\n(nu=1.0)", "Underdispersion\n(nu=1.5)")
models_wr  <- c("ZICMP", "ZIP", "ZINB", "ZIGP")
model_colors <- c(BLUE, ORANGE, AQUA, YELLOW)
winrate <- rbind(
  c(0.70, 0.00, 0.30, 0.00),
  c(0.00, 0.90, 0.10, 0.00),
  c(0.90, 0.10, 0.00, 0.00)
)
colnames(winrate) <- models_wr
rownames(winrate) <- regimes_wr

plot_winrate <- function() {
  bp <- barplot(t(winrate), beside = TRUE, col = model_colors, ylim = c(0, 1.05),
                ylab = "Fraction of replicates with lowest WAIC2 (R=10, n=600)",
                names.arg = regimes_wr, legend.text = models_wr,
                args.legend = list(x = "top", horiz = TRUE, bty = "n", inset = c(0, -0.08)),
                main = "Which model wins, by dispersion regime")
  for (i in seq_len(nrow(winrate))) {
    for (j in seq_len(ncol(winrate))) {
      v <- winrate[i, j]
      if (v > 0) text(bp[j, i], v + 0.03, sprintf("%.2f", v), cex = 0.75)
    }
  }
}

pdf("fig_modelcomparison.pdf", width = 7.5, height = 4.2); plot_winrate(); dev.off()
png("fig_modelcomparison.png", width = 1900, height = 1050, res = 220); plot_winrate(); dev.off()

cat("Wrote fig_recovery.{pdf,png}, fig_bias.{pdf,png}, fig_modelcomparison.{pdf,png}",
    "\n(and, if sim_summary_table_ext.csv was found, fig_recovery_extended.{pdf,png},",
    "fig_bias_extended.{pdf,png} covering all five dispersion regimes).\n")
