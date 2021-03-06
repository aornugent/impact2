// M2 - Joint species tobit with varying interactions.

// Species respond to environmental gradients of fertility and rainfall the
// same across all treatements and interact with each other differently
// in each experimental treatment.

// author: Andrew O'Reilly-Nugent
// email:  andrew.oreilly-nugent@canberra.edu.au
// date:   22/01/18

data{
  int<lower=1> N;
  int<lower=0> N_censored;
  int<lower=1> S;
  int<lower=1> K;
  int<lower=1> E;
  int<lower=1> N_plots;
  int<lower=1> N_sites;

  vector[K] X[N];                             //environmental covariates
  vector<lower=0>[S] y_observed[N];           //left censored data

  int<lower=1, upper=E> treatment[N];         //treatment index
  int<lower=1, upper=N_plots> plot[N];        //plot index
  int<lower=1, upper=N_sites> site[N_plots];  //site index
  int<lower=1> shape_prior;                   //LKJ prior shape parameter
}
parameters{
  matrix[S, K] B;                             //scaled species coefficients
  vector[K] mu_K;
  vector<lower=0>[K] sigma_B;

  vector[N_plots] B_plot;                     //scaled random effects
  real<lower=0> sigma_plot;

  vector[N_sites] B_site;
  real<lower=0> sigma_site;

  cholesky_factor_corr[S] L_Omega[E];         //correlation cholesky factor
  vector<lower=0>[S] sigma[E];                   //observation variances

  real<upper=0> y_censored[N_censored];       //censored values as parameters
}
transformed parameters{
  cholesky_factor_cov[S] L_Sigma[E];          //covariance cholesky factor

  // Convert correlation to covariance
  for(i in 1:E){
    L_Sigma[i] = diag_pre_multiply(sigma[i], L_Omega[i]);
  }
}
model{
  vector[S] y_latent[N];                      //uncensored observations
  int pos = 1;

  // Uncensor data
  for(i in 1:N) {
    for (j in 1:S) {
      if (y_observed[i, j] > 0) {
        y_latent[i, j] = y_observed[i, j];
      }
      else {
        y_latent[i, j] = y_censored[pos];
        pos = pos + 1;
      }
    }

    // Likelihood
    y_latent[i] ~ multi_normal_cholesky(B * X[i] + B_plot[plot[i]], L_Sigma[treatment[i]]);
  }

  // Priors
  // Coefficients have central means and variance.
  for(j in 1:K){
    B[ , j] ~ student_t(4, mu_K[j], sigma_B[j]);
  }

  mu_K ~ normal(0, 1);
  sigma_B ~ normal(0, 1);

  // Quadrats are nested within plots, within sites.
  B_plot ~ normal(B_site[site], sigma_plot);
  sigma_plot ~ cauchy(0, 1);

  B_site ~ normal(0, sigma_site);
  sigma_site ~ cauchy(0, 1);

  // Correlations and covarinces change between treatments.
  for(i in 1:E){
    L_Omega[i] ~ lkj_corr_cholesky(shape_prior);
    sigma[i] ~ cauchy(0, 1);
  }
}
generated quantities{
  vector[S] y_pred[N];                        //predicted fit
  corr_matrix[S] Omega[E];                    //correlation matrix
  cov_matrix[S] Sigma[E];                     //covariance matrix

  for(i in 1:N){
    y_pred[i] = multi_normal_cholesky_rng(B * X[i] + B_plot[plot[i]], L_Sigma[treatment[i]]);
  }

  for(i in 1:E){
    Omega[i] = multiply_lower_tri_self_transpose(L_Omega[i]);
    Sigma[i] = quad_form_diag(Omega[i], sigma[i]);
  }
}
