/**
 * ZINAR(1) with Conway-Maxwell Poisson (CMP) innovations
 * Optimized version: precomputed log_Z, log_lambda, lgamma table
 * Based on COMPoissonReg (Raim & Sellers, 2022)
 * https://www.census.gov/content/dam/Census/library/working-papers/2022/adrm/RRC2022-01.pdf
 *
 * LEGACY / REFERENCE ONLY: this is the natural (lambda, nu) parameterization.
 * It showed a strong posterior ridge (cor(lambda, nu) ~ 0.9) and is kept
 * here for reference; every script in this repo uses
 * ZIINAR1-CMP-fast-reparam.stan instead. See that file's header for why.
 */
functions {
  real cmp_log_Z_asymp(real lambda, real nu) {
    real log_lambda = log(lambda);
    return nu * exp(log_lambda / nu)
           - ((nu - 1.0) / (2.0 * nu)) * log_lambda
           - ((nu - 1.0) / 2.0) * log(2.0 * pi())
           - 0.5 * log(nu);
  }
  real cmp_log_Z_trunc(real log_lambda, real nu, int M,
                       vector lgam) {
    vector[M + 1] log_terms;
    for (r in 0:M) {
      log_terms[r + 1] = r * log_lambda - nu * lgam[r + 1];
    }
    return log_sum_exp(log_terms);
  }
  real cmp_log_Z(real lambda, real nu, int M,
                 real hybrid_tol, vector lgam) {
    real log_lambda = log(lambda);
    real test       = exp(-log_lambda / nu);
    if (test < hybrid_tol) {
      return cmp_log_Z_asymp(lambda, nu);
    } else {
      return cmp_log_Z_trunc(log_lambda, nu, M, lgam);
    }
  }
  int cmp_rng(real lambda, real nu, int M) {
    vector[M + 1] log_terms;
    real log_Z;
    real u;
    real cdf;
    real log_lambda;
    if (lambda <= 0) reject("lambda must be positive.");
    if (nu     <= 0) reject("nu must be positive.");
    if (M      <  1) reject("M must be at least 1.");
    log_lambda = log(lambda);
    for (r in 0:M) {
      log_terms[r + 1] = r * log_lambda - nu * lgamma(r + 1);
    }
    log_Z = log_sum_exp(log_terms);
    u     = uniform_rng(0, 1);
    cdf   = 0;
    for (y in 0:M) {
      cdf += exp(log_terms[y + 1] - log_Z);
      if (u <= cdf) return y;
    }
    reject("CMP RNG failed: truncated CDF did not reach 1. Increase M.");
    return -1;
  }
}
data {
  int<lower=0>        T;
  array[T] int<lower=0> y;
  int<lower=1>        M;
  real<lower=0>       hybrid_tol;
  int<lower=0>        ff;
}
transformed data {
  // ── Precompute lgamma table ONCE before sampling ──────────────
  // covers 0, 1, ..., max(y)+1  (need +1 for safety)
  int y_max = max(y);
  vector[M + 2] lgam;          // covers 0:M for log_Z truncation
  for (k in 0:(M + 1)) {
    lgam[k + 1] = lgamma(k + 1);
  }
  // ── Precompute min(y[t-1], y[t]) for t=2:T ───────────────────
  array[T - 1] int pp;
  for (t in 2:T) {
    pp[t - 1] = min(y[t - 1], y[t]);
  }
}
parameters {
  real<lower=0, upper=1> alpha;
  real<lower=0>          lambda;
  real<lower=0>          nu;
  real<lower=0, upper=1> rho;
}
transformed parameters {
  vector[T] log_mu;
  // ── Precompute ONCE per gradient evaluation ───────────────────
  real log_Z    = cmp_log_Z(lambda, nu, M, hybrid_tol, lgam);
  real log_lam  = log(lambda);
  real log_rho  = log(rho);
  real log1mrho = log1m(rho);
  log_mu[1] = log(y[1] + 1e-10);   // initialise (not used in likelihood)
  for (t in 2:T) {
    int  p      = pp[t - 1];        // precomputed min(y[t-1], y[t])
    int  yt     = y[t];
    int  yt1    = y[t - 1];
    // ── j = 0 term ─────────────────────────────────────────────
    real lbin0  = binomial_lpmf(0 | yt1, alpha);
    real lcmp0  = yt * log_lam - nu * lgam[yt + 1] - log_Z;
    real lterm0;
    if (yt == 0)
      lterm0 = lbin0 + log_sum_exp(log_rho, log1mrho + lcmp0);
    else
      lterm0 = lbin0 + log1mrho + lcmp0;
    // ── j = 1:p terms ──────────────────────────────────────────
    if (p == 0) {
      log_mu[t] = lterm0;
    } else {
      vector[p + 1] lterms;
      lterms[1] = lterm0;
      for (j in 1:p) {
        real lbinj  = binomial_lpmf(j | yt1, alpha);
        int  diff   = yt - j;
        real lcmpj  = diff * log_lam - nu * lgam[diff + 1] - log_Z;
        if (yt == j)
          lterms[j + 1] = lbinj + log_sum_exp(log_rho,
                                               log1mrho + lcmpj);
        else
          lterms[j + 1] = lbinj + log1mrho + lcmpj;
      }
      log_mu[t] = log_sum_exp(lterms);
    }
  }
}
model {
  alpha  ~ uniform(0, 1);
  rho    ~ uniform(0, 1);
  //lambda ~ lognormal(0, 2);
  //nu     ~ lognormal(0, 2);
  lambda ~ student_t(5, 0, 5);
  nu ~ student_t(5, 0, 5);
  // ── Vectorized log likelihood ──────────────────────────────
  target += sum(log_mu[2:T]);
}
generated quantities {
  // ── Posterior predictive forecast ─────────────────────────
  array[ff + 1] int y_pred;
  y_pred[1] = y[T];
  for (t in 2:(ff + 1)) {
    y_pred[t] = binomial_rng(y_pred[t - 1], alpha);
    int aa     = bernoulli_rng(rho);
    if (aa == 0)
      y_pred[t] += cmp_rng(lambda, nu, M);
  }
  // ── Pointwise log likelihood ───────────────────────────────
  vector[T] log_lik;
  log_lik[1] = 0;
  for (t in 2:T)
    log_lik[t] = log_mu[t];
  // ── Model selection criteria ───────────────────────────────
  real ll  = sum(log_lik[2:T]);
  real aic = -2 * ll + 2  * 4;          // 4 parameters
  real bic = -2 * ll + 4  * log(T - 1);
}
