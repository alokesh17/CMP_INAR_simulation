/**
 * ZINAR(1) with Conway-Maxwell Poisson (CMP) innovations
 * Reparameterized version: samples (mu_cmp, nu) instead of (lambda, nu),
 * with lambda := mu_cmp^nu recovered in transformed parameters.
 *
 * WHY: the natural (lambda, nu) parameterization exhibits a strong ridge
 * (empirically cor(lambda, nu) ~ 0.9 in posterior draws here), a known
 * pathology of CMP models generally. mu_cmp ~ lambda^{1/nu} is
 * approximately the CMP mean/location, so (mu_cmp, nu) is much closer to
 * an orthogonal (location, dispersion) parameterization.
 *
 * Everything else (functions, data, transformed data, generated quantities)
 * is unchanged from ZIINAR1-CMP-fast.stan.
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
  int y_max = max(y);
  vector[M + 2] lgam;
  for (k in 0:(M + 1)) {
    lgam[k + 1] = lgamma(k + 1);
  }
  array[T - 1] int pp;
  for (t in 2:T) {
    pp[t - 1] = min(y[t - 1], y[t]);
  }
}
parameters {
  real<lower=0, upper=1> alpha;
  real<lower=0>          mu_cmp;   // ~ CMP mean scale; lambda = mu_cmp^nu
  real<lower=0>          nu;
  real<lower=0, upper=1> rho;
}
transformed parameters {
  real lambda   = mu_cmp ^ nu;     // recover the natural CMP rate parameter
  vector[T] log_mu;
  real log_Z    = cmp_log_Z(lambda, nu, M, hybrid_tol, lgam);
  real log_lam  = log(lambda);
  real log_rho  = log(rho);
  real log1mrho = log1m(rho);
  log_mu[1] = log(y[1] + 1e-10);
  for (t in 2:T) {
    int  p      = pp[t - 1];
    int  yt     = y[t];
    int  yt1    = y[t - 1];
    real lbin0  = binomial_lpmf(0 | yt1, alpha);
    real lcmp0  = yt * log_lam - nu * lgam[yt + 1] - log_Z;
    real lterm0;
    if (yt == 0)
      lterm0 = lbin0 + log_sum_exp(log_rho, log1mrho + lcmp0);
    else
      lterm0 = lbin0 + log1mrho + lcmp0;
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
  mu_cmp ~ student_t(5, 0, 5);   // prior now on the (approx) CMP mean directly
  nu     ~ student_t(5, 0, 5);
  target += sum(log_mu[2:T]);
}
generated quantities {
  array[ff + 1] int y_pred;
  y_pred[1] = y[T];
  for (t in 2:(ff + 1)) {
    y_pred[t] = binomial_rng(y_pred[t - 1], alpha);
    int aa     = bernoulli_rng(rho);
    if (aa == 0)
      y_pred[t] += cmp_rng(lambda, nu, M);
  }
  vector[T] log_lik;
  log_lik[1] = 0;
  for (t in 2:T)
    log_lik[t] = log_mu[t];
  real ll  = sum(log_lik[2:T]);
  real aic = -2 * ll + 2  * 4;
  real bic = -2 * ll + 4  * log(T - 1);
}
