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
  // data for splines
  int Ks;  // number of linear effects
  matrix[N, Ks] Xs;  // design matrix for the linear effects
  // data for spline t2(time)
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
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
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  vector[Ks] bs;  // spline coefficients
  // parameters for spline t2(time)
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  real<lower=0> sds_1_1;  // standard deviations of spline coefficients
  real<lower=0> sigma;  // dispersion parameter
  real<lower=1> nu;  // degrees of freedom or shape
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1;
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N_1] r_1_2;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  // compute actual spline coefficients
  s_1_1 = sds_1_1 * zs_1_1;
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_1_2 = (sd_1[2] * (z_1[2]));
  lprior += student_t_lpdf(Intercept | 3, 0.3, 2.5);
  lprior += student_t_lpdf(sds_1_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu | 2, 0.1)
    - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 2 * student_t_lccdf(0 | 3, 0, 2.5);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + rep_vector(0.0, N) + Xs * bs + Zs_1_1 * s_1_1;
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
  target += std_normal_lpdf(zs_1_1);
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_1[2]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}
