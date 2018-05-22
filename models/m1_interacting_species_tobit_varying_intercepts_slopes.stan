// M4 - Joint species tobit with varying intercepts and slopes..

// Species respond to environmental gradients of fertility and rainfall the
// differently in each treatement, and interact with each other the same in
// each treatment.

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
   int<lower=1, upper=N_plots> plot[N];  //plot index
  int<lower=1, upper=N_sites> site[N_plots];  //site index
  int<lower=1> shape_prior;                   //LKJ prior shape parameter
}
parameters{
  matrix[S, K] B[E];                          //scaled species coefficients
  vector[K] mu_B[E];
  vector<lower=0>[K] sigma_B[E];
  vector[K] mu_K;

  vector[N_plots] B_plot;                    //scaled random effects
  real<lower=0> sigma_plot;

  vector[N_sites] B_site;
  real<lower=0> sigma_site;

  cholesky_factor_corr[S] L_Omega;            //correlation cholesky factor
  vector<lower=0>[S] sigma;                   //observation variances

  real<upper=0> y_censored[N_censored];       //censored values as parameters
}
transformed parameters{
  cholesky_factor_cov[S] L_Sigma;             //covariance cholesky factor

  // Convert correlation to covariance
  L_Sigma = diag_pre_multiply(sigma, L_Omega);
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
    y_latent[i] ~ multi_normal_cholesky(B[treatment[i]] * X[i] + B_plot[plot[i]], L_Sigma);
  }

  // Priors
  // Coefficients have central means and variances in each treatment.
  // Priors
  for(i in 1:E){
    for(j in 1:K){
      B[i, , j] ~ student_t(4, mu_B[i, j], sigma_B[i, j]);

      mu_B[i, j] ~ normal(mu_K[j], 1);
      sigma_B[i, j] ~ normal(0, 1);
    }
  }

  mu_K ~ normal(0, 1);

  // Plots nested within sites.
  B_plot ~ normal(B_site[site], sigma_plot);
  sigma_plot ~ cauchy(0, 1);

  B_site ~ normal(0, sigma_site);
  sigma_site ~ cauchy(0, 1);

  // Correlations and covariances are constant.
  L_Omega ~ lkj_corr_cholesky(shape_prior);
  sigma ~ cauchy(0, 1);
}
generated quantities{
  vector[S] y_pred[N];                        //predicted fit
  corr_matrix[S] Omega;                       //correlation matrix
  cov_matrix[S] Sigma;                        //covariance matrix

  for(i in 1:N){
    y_pred[i] = multi_normal_cholesky_rng(B[treatment[i]] * X[i] + B_plot[plot[i]], L_Sigma);
  }

  Omega = multiply_lower_tri_self_transpose(L_Omega);
  Sigma = quad_form_diag(Omega, sigma);
}
