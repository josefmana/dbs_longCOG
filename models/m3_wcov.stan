"x"
"1" "// generated with brms 2.17.0
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
  int<lower=1> N_drs;  // number of observations
  vector[N_drs] Y_drs;  // response variable
  int<lower=-1,upper=2> cens_drs[N_drs];  // indicates censoring
  int<lower=1> K_drs;  // number of population-level effects
  matrix[N_drs, K_drs] X_drs;  // population-level design matrix
  int<lower=1> Ksp_drs;  // number of special effects terms
  int<lower=1> N_bdi;  // number of observations
  vector[N_bdi] Y_bdi;  // response variable
  int<lower=0> Nmi_bdi;  // number of missings
  int<lower=1> Jmi_bdi[Nmi_bdi];  // positions of missings
  int<lower=1> K_bdi;  // number of population-level effects
  matrix[N_bdi, K_bdi] X_bdi;  // population-level design matrix
  int<lower=1> Ksp_bdi;  // number of special effects terms
  int<lower=1> N_led;  // number of observations
  vector[N_led] Y_led;  // response variable
  int<lower=0> Nmi_led;  // number of missings
  int<lower=1> Jmi_led[Nmi_led];  // positions of missings
  // data for splines
  int Ks_led;  // number of linear effects
  matrix[N_led, Ks_led] Xs_led;  // design matrix for the linear effects
  // data for spline t2(time)
  int nb_led_1;  // number of bases
  int knots_led_1[nb_led_1];  // number of knots
  // basis function matrices
  matrix[N_led, knots_led_1[1]] Zs_led_1_1;
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1_drs[N_drs];  // grouping indicator per observation
  // group-level predictor values
  vector[N_drs] Z_1_drs_1;
  vector[N_drs] Z_1_drs_2;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2_bdi[N_bdi];  // grouping indicator per observation
  // group-level predictor values
  vector[N_bdi] Z_2_bdi_1;
  vector[N_bdi] Z_2_bdi_2;
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  int<lower=1> J_3_led[N_led];  // grouping indicator per observation
  // group-level predictor values
  vector[N_led] Z_3_led_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc_drs = K_drs - 1;
  matrix[N_drs, Kc_drs] Xc_drs;  // centered version of X_drs without an intercept
  vector[Kc_drs] means_X_drs;  // column means of X_drs before centering
  int Kc_bdi = K_bdi - 1;
  matrix[N_bdi, Kc_bdi] Xc_bdi;  // centered version of X_bdi without an intercept
  vector[Kc_bdi] means_X_bdi;  // column means of X_bdi before centering
  for (i in 2:K_drs) {
    means_X_drs[i - 1] = mean(X_drs[, i]);
    Xc_drs[, i - 1] = X_drs[, i] - means_X_drs[i - 1];
  }
  for (i in 2:K_bdi) {
    means_X_bdi[i - 1] = mean(X_bdi[, i]);
    Xc_bdi[, i - 1] = X_bdi[, i] - means_X_bdi[i - 1];
  }
}
parameters {
  vector[Kc_drs] b_drs;  // population-level effects
  real Intercept_drs;  // temporary intercept for centered predictors
  vector[Ksp_drs] bsp_drs;  // special effects coefficients
  real<lower=0> sigma_drs;  // dispersion parameter
  real<lower=1> nu_drs;  // degrees of freedom or shape
  vector[Nmi_bdi] Ymi_bdi;  // estimated missings
  vector[Kc_bdi] b_bdi;  // population-level effects
  real Intercept_bdi;  // temporary intercept for centered predictors
  vector[Ksp_bdi] bsp_bdi;  // special effects coefficients
  real<lower=0> sigma_bdi;  // dispersion parameter
  vector[Nmi_led] Ymi_led;  // estimated missings
  real Intercept_led;  // temporary intercept for centered predictors
  vector[Ks_led] bs_led;  // spline coefficients
  // parameters for spline t2(time)
  // standarized spline coefficients
  vector[knots_led_1[1]] zs_led_1_1;
  real<lower=0> sds_led_1_1;  // standard deviations of spline coefficients
  real<lower=0> sigma_led;  // dispersion parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // standardized group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  vector[N_3] z_3[M_3];  // standardized group-level effects
}
transformed parameters {
  // actual spline coefficients
  vector[knots_led_1[1]] s_led_1_1;
  vector[N_1] r_1_drs_1;  // actual group-level effects
  vector[N_1] r_1_drs_2;  // actual group-level effects
  vector[N_2] r_2_bdi_1;  // actual group-level effects
  vector[N_2] r_2_bdi_2;  // actual group-level effects
  vector[N_3] r_3_led_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  // compute actual spline coefficients
  s_led_1_1 = sds_led_1_1 * zs_led_1_1;
  r_1_drs_1 = (sd_1[1] * (z_1[1]));
  r_1_drs_2 = (sd_1[2] * (z_1[2]));
  r_2_bdi_1 = (sd_2[1] * (z_2[1]));
  r_2_bdi_2 = (sd_2[2] * (z_2[2]));
  r_3_led_1 = (sd_3[1] * (z_3[1]));
  lprior += normal_lpdf(b_drs[1] | 0, 0.1);
  lprior += normal_lpdf(b_drs[2] | -0.2, 0.1);
  lprior += normal_lpdf(b_drs[3] | 0, 0.1);
  lprior += normal_lpdf(b_drs[4] | 0, 0.1);
  lprior += normal_lpdf(b_drs[5] | 0, 0.1);
  lprior += normal_lpdf(b_drs[6] | 0, 0.1);
  lprior += normal_lpdf(b_drs[7] | 0, 0.1);
  lprior += normal_lpdf(b_drs[8] | 0, 0.1);
  lprior += normal_lpdf(b_drs[9] | 0, 0.1);
  lprior += normal_lpdf(b_drs[10] | 0, 0.1);
  lprior += normal_lpdf(b_drs[11] | 0, 0.1);
  lprior += normal_lpdf(b_drs[12] | 0, 0.1);
  lprior += normal_lpdf(b_drs[13] | 0, 0.1);
  lprior += normal_lpdf(b_drs[14] | 0, 0.1);
  lprior += normal_lpdf(b_drs[15] | 0, 0.1);
  lprior += normal_lpdf(b_drs[16] | 0, 0.1);
  lprior += normal_lpdf(Intercept_drs | 0.3, 0.1);
  lprior += normal_lpdf(bsp_drs[1] | 0, 0.1);
  lprior += normal_lpdf(bsp_drs[2] | 0, 0.1);
  lprior += exponential_lpdf(sigma_drs | 1);
  lprior += gamma_lpdf(nu_drs | 2, 0.1)
    - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(b_bdi[1] | 0, 0.5);
  lprior += normal_lpdf(b_bdi[2] | 0, 0.5);
  lprior += normal_lpdf(Intercept_bdi | 0.6, 0.5);
  lprior += normal_lpdf(bsp_bdi[1] | 0, 0.5);
  lprior += exponential_lpdf(sigma_bdi | 1);
  lprior += normal_lpdf(Intercept_led | 0, 100);
  lprior += normal_lpdf(bs_led[1] | 0, 100);
  lprior += student_t_lpdf(sds_led_1_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += exponential_lpdf(sigma_led | 1);
  lprior += normal_lpdf(sd_1[1] | 0, 0.1)
    - 1 * normal_lccdf(0 | 0, 0.1);
  lprior += normal_lpdf(sd_1[2] | 0, 0.1)
    - 1 * normal_lccdf(0 | 0, 0.1);
  lprior += normal_lpdf(sd_2[1] | 0, 0.5)
    - 1 * normal_lccdf(0 | 0, 0.5);
  lprior += normal_lpdf(sd_2[2] | 0, 0.5)
    - 1 * normal_lccdf(0 | 0, 0.5);
  lprior += normal_lpdf(sd_3[1] | 0, 0.5)
    - 1 * normal_lccdf(0 | 0, 0.5);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // vector combining observed and missing responses
    vector[N_bdi] Yl_bdi = Y_bdi;
    // vector combining observed and missing responses
    vector[N_led] Yl_led = Y_led;
    // initialize linear predictor term
    vector[N_drs] mu_drs = Intercept_drs + Xc_drs * b_drs;
    // initialize linear predictor term
    vector[N_bdi] mu_bdi = Intercept_bdi + Xc_bdi * b_bdi;
    // initialize linear predictor term
    vector[N_led] mu_led = Intercept_led + rep_vector(0.0, N_led) + Xs_led * bs_led + Zs_led_1_1 * s_led_1_1;
    Yl_bdi[Jmi_bdi] = Ymi_bdi;
    Yl_led[Jmi_led] = Ymi_led;
    for (n in 1:N_drs) {
      // add more terms to the linear predictor
      mu_drs[n] += (bsp_drs[1]) * Yl_bdi[n] + (bsp_drs[2]) * Yl_led[n] + r_1_drs_1[J_1_drs[n]] * Z_1_drs_1[n] + r_1_drs_2[J_1_drs[n]] * Z_1_drs_2[n];
    }
    for (n in 1:N_bdi) {
      // add more terms to the linear predictor
      mu_bdi[n] += (bsp_bdi[1]) * Yl_led[n] + r_2_bdi_1[J_2_bdi[n]] * Z_2_bdi_1[n] + r_2_bdi_2[J_2_bdi[n]] * Z_2_bdi_2[n];
    }
    for (n in 1:N_led) {
      // add more terms to the linear predictor
      mu_led[n] += r_3_led_1[J_3_led[n]] * Z_3_led_1[n];
    }
    for (n in 1:N_drs) {
    // special treatment of censored data
      if (cens_drs[n] == 0) {
        target += student_t_lpdf(Y_drs[n] | nu_drs, mu_drs[n], sigma_drs);
      } else if (cens_drs[n] == 1) {
        target += student_t_lccdf(Y_drs[n] | nu_drs, mu_drs[n], sigma_drs);
      } else if (cens_drs[n] == -1) {
        target += student_t_lcdf(Y_drs[n] | nu_drs, mu_drs[n], sigma_drs);
      }
    }
    target += normal_lpdf(Yl_bdi | mu_bdi, sigma_bdi);
    target += normal_lpdf(Yl_led | mu_led, sigma_led);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zs_led_1_1);
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_1[2]);
  target += std_normal_lpdf(z_2[1]);
  target += std_normal_lpdf(z_2[2]);
  target += std_normal_lpdf(z_3[1]);
}
generated quantities {
  // actual population-level intercept
  real b_drs_Intercept = Intercept_drs - dot_product(means_X_drs, b_drs);
  // actual population-level intercept
  real b_bdi_Intercept = Intercept_bdi - dot_product(means_X_bdi, b_bdi);
  // actual population-level intercept
  real b_led_Intercept = Intercept_led;
}
"
