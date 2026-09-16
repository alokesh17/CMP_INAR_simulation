# CMP-INAR: Zero-Inflated Conway-Maxwell-Poisson INAR(1)

Bayesian estimation and model comparison for a zero-inflated integer-valued
autoregressive process of order 1 (INAR(1)) with Conway-Maxwell-Poisson (CMP)
innovations. The CMP innovation lets a single model span underdispersion,
equidispersion, and overdispersion through one shape parameter `nu`, instead
of committing in advance to Poisson (`nu=1` only), Negative Binomial
(overdispersion only), or Generalized Poisson.

The process:

```
Y_t = alpha o Y_{t-1} + epsilon_t,      epsilon_t ~ ZICMP(lambda, nu, rho)
```

where `alpha o Y_{t-1}` is binomial thinning and `ZICMP(lambda, nu, rho)` is a
zero-inflated CMP innovation with zero-inflation probability `rho`.

Estimation is Bayesian via Stan, using a hybrid asymptotic/truncated-series
approximation to the intractable CMP normalizing constant. Model comparison
against Poisson, Negative Binomial, and Generalized Poisson innovations
(fit via the `ZIHINAR1` package) uses EAIC, EBIC, DIC, WAIC1, and WAIC2,
computed identically across all four models.

Authors: Alokesh Manna, Dipak K. Dey, Víctor Lachos.

The manuscript itself is kept in a separate, not-yet-public location and is
not part of this repository; this repo is the reproducible code (Stan model
+ R pipeline) behind it.

## Repository layout

```
CMP_INAR_simulation/
├── stan/
│   ├── ZIINAR1-CMP-fast-reparam.stan   <- THE model every script uses
│   └── ZIINAR1-CMP-fast.stan           <- legacy natural parameterization (reference only)
├── R/
│   ├── run_simulation_grid.R           <- Table 1: parameter recovery, original 3-regime grid
│   ├── run_simulation_grid_extended.R  <- Table 1, extended to 5 regimes (nu = 0.2,0.5,1,1.5,2.5)
│   ├── compare_innovations.R           <- Table 2/3: model-selection simulation, original 3-regime
│   ├── compare_innovations_extended.R  <- Table 2/3, extended 5-regime x 4-n grid setup
│   ├── part4_resume.R                  <- standalone resume/completion script for the extended grid
│   ├── ZICMP_INAR_simulation_plot_predictive.R  <- merged Parts 1-3 driver + plotting + predictive example
│   ├── make_sim_figures.R / .py        <- figure generation from simulation outputs
│   ├── real_data_analysis.R            <- Section 5: sexoffences, violence, claims, soap
│   ├── real_data_analysis_polio_campy.R <- companion real-data script: polio, campylobacterosis
│   ├── forecast_evaluation.R           <- Table 4: out-of-sample forecasting evaluation
│   └── predictive_analysis.R           <- posterior predictive checks (Section 3.4 worked example)
└── docs/
    └── review_response_summary.md      <- notes on responses to review/feedback
```

## Getting started

### 1. Download the repository

Clone it with git (recommended, so you can pull future updates):

```
git clone https://github.com/<your-username>/CMP_INAR_simulation.git
cd CMP_INAR_simulation
```

Or, if you just want a one-off copy without git: on the repo's GitHub page,
click **Code -> Download ZIP**, then unzip it and `cd` into the folder.

### 2. Install R and the required packages

Requires R (>= 4.0), a C++ toolchain (needed by `rstan` to compile Stan
models), and:

```r
install.packages(c("rstan", "COMPoissonReg", "matrixStats",
                    "gamlss.data", "tscount"))
# ZIHINAR1 (Poisson/NB/GP competitor fits + matching model-selection criteria):
install.packages("ZIHINAR1")
# or, if not yet on CRAN for your R version:
# remotes::install_github("fushengyy/ZIHINAR1")
```

On macOS, `rstan` additionally needs Xcode command line tools
(`xcode-select --install`) and, for best performance, a configured
`~/.R/Makevars` — see the RStan getting-started guide if `install.packages("rstan")`
fails to compile:
https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

### 3. Point the scripts at the Stan model

Every script in `R/` compiles `ZIINAR1-CMP-fast-reparam.stan` via a
`STAN_FILE` (or `CMP_STAN_FILE`) variable near the top. The simplest setup is
to run scripts from inside `R/` with a copy of (or symlink to) the Stan file
alongside them:

```
cd R
ln -s ../stan/ZIINAR1-CMP-fast-reparam.stan .
```

or just edit that variable at the top of whichever script you're running to
point at `"../stan/ZIINAR1-CMP-fast-reparam.stan"`.

### 4. Run something small first

Before any full run, open `R/run_simulation_grid.R` and set
`RUN_MODE <- "quick"` (already the safest default in most scripts), then:

```
Rscript run_simulation_grid.R
```

This should compile the Stan model (can take a minute or two the first
time) and finish a tiny 1-cell, 3-replicate run in well under a minute. If
that works, the environment is set up correctly and you can move on to the
full pipeline below.

## What to run, in order

1. **Parameter recovery (Table 1).**
   `R/run_simulation_grid.R` (or the 5-regime `run_simulation_grid_extended.R`
   for the version actually reported in the draft). Set `RUN_MODE` to
   `"quick"` first to confirm the Stan model compiles, then `"pilot"`, then
   `"full"`. Checkpoints to `sim_results_partial.rds` every `(nu, n)` cell, so
   a killed run resumes without repeating finished cells. Emits
   `sim_summary_table.csv` and `sim_table_body.tex`.

2. **Model-selection simulation (Table 2/3, win-rate).**
   `R/compare_innovations.R` for the original 3-regime design, or
   `R/compare_innovations_extended.R` + `R/part4_resume.R` for the 5-regime
   x 4-sample-size grid whose numbers appear in the current draft. `part4_resume.R`
   is a standalone, self-contained script meant to be run fresh (it recompiles
   the Stan model and reloads `cmp_cache.rds` itself) to finish or resume an
   interrupted extended grid. Emits `compare_summary_table.csv`,
   `compare_win_rate.csv`, and ready-to-paste LaTeX table bodies.

3. **Real-data application.**
   `R/real_data_analysis.R` fits ZICMP/ZIP/ZINB/ZIGP to four real series
   (sexoffences, family violence, workers'-compensation claims, soap sales)
   and emits posterior estimate and model-comparison LaTeX table bodies.
   `R/real_data_analysis_polio_campy.R` does the same for two additional
   series (US poliomyelitis, Quebec campylobacterosis), including the
   lambda-reliability and nu-CI-excludes-1 diagnostics used to decide whether
   a series's results are numerically trustworthy.

4. **Out-of-sample forecasting.**
   `R/forecast_evaluation.R`. Trains each model on the first 80% of each
   series, forecasts the held-out 20% one step ahead, and reports LPS, MAE,
   RMSE, coverage, and zero-event Brier score. Includes a self-check that
   drops any model/series cell where a hand-rolled log-likelihood can't be
   verified against Stan's own posterior-averaged log-likelihood within
   tolerance, rather than reporting an unverified forecast number.

5. **Figures.**
   `R/make_sim_figures.R` (or the Python port `make_sim_figures.py`) builds
   the win-rate and comparison figures from the CSVs produced in steps 1-2.

Every script prints its progress to the console and checkpoints intermediate
results to `.rds`/`.csv` files in whatever directory you run it from, so
long grid runs can be interrupted and picked back up (see the Notes below).

## Notes

- `ZIINAR1-CMP-fast-reparam.stan` samples `(mu_cmp, nu)` with
  `lambda := mu_cmp^nu`, rather than sampling `(lambda, nu)` directly. The
  natural parameterization showed a strong posterior ridge
  (`cor(lambda, nu) ~ 0.9`), a known CMP identifiability issue, not a bug;
  the reparameterization is close to an orthogonal (location, dispersion)
  parameterization and is what every script in this repo uses.
- Long-running grid scripts checkpoint to `.rds` files after every cell so
  they can be killed and resumed; generated `.rds`/`.csv`/log files are not
  tracked in this repo (see `.gitignore`) since they're reproducible from the
  scripts and can be large.
- This repository does not yet include a LICENSE file. Add one (or ask your
  co-authors/institution what's appropriate) before making the repo public,
  since the default with no LICENSE is "all rights reserved."
- This repo is currently **private**. If you make it public later, double
  check `git log --stat` first to confirm no unpublished paper draft or
  results you don't want public ended up in the history.
