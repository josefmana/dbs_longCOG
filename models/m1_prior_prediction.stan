// generated with brms 2.17.0
functions {
  /* compute the logm1 link
   * Args:
   *   p: a positive scalar
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real logm1(real y) {
     return log(y - 1);
   }
  /* compute the inverse of the logm1 link
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a positive scalar
   */
   real expp1(real y) {
     return exp(y) + 1;
   }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=-1,upper=2> cens[N];  // indicates censoring
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // dispersion parameter
  real<lower=1> nu;  // degrees of freedom or shape
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N_1] r_1_2;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_1_2 = (sd_1[2] * (z_1[2]));
  lprior += normal_lpdf(b[1] | -0.2, 0.1);
  lprior += normal_lpdf(b[2] | 0, 0.1);
  lprior += normal_lpdf(b[3] | 0, 0.1);
  lprior += normal_lpdf(b[4] | 0, 0.1);
  lprior += normal_lpdf(b[5] | 0, 0.1);
  lprior += normal_lpdf(b[6] | 0, 0.1);
  lprior += normal_lpdf(b[7] | 0, 0.1);
  lprior += normal_lpdf(b[8] | 0, 0.1);
  lprior += normal_lpdf(b[9] | 0, 0.1);
  lprior += normal_lpdf(b[10] | 0, 0.1);
  lprior += normal_lpdf(b[11] | 0, 0.1);
  lprior += normal_lpdf(b[12] | 0, 0.1);
  lprior += normal_lpdf(b[13] | 0, 0.1);
  lprior += normal_lpdf(b[14] | 0, 0.1);
  lprior += normal_lpdf(b[15] | 0, 0.1);
  lprior += normal_lpdf(Intercept | 0.3, 0.1);
  lprior += exponential_lpdf(sigma | 1);
  lprior += gamma_lpdf(nu | 2, 0.1)
    - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(sd_1[1] | 0, 0.1)
    - 1 * normal_lccdf(0 | 0, 0.1);
  lprior += normal_lpdf(sd_1[2] | 0, 0.1)
    - 1 * normal_lccdf(0 | 0, 0.1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n];
    }
    for (n in 1:N) {
    // special treatment of censored data
      if (cens[n] == 0) {
        target += student_t_lpdf(Y[n] | nu, mu[n], sigma);
      } else if (cens[n] == 1) {
        target += student_t_lccdf(Y[n] | nu, mu[n], sigma);
      } else if (cens[n] == -1) {
        target += student_t_lcdf(Y[n] | nu, mu[n], sigma);
      }
    }
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_1[2]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
