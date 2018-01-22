// M1 - Independent species tobit

// Species respond to environmental gradients of fertility and rainfall,
// but not to experimental treatment and do not interact.

// author: Andrew O'Reilly-Nugent
// email:  andrew.oreilly-nugent@canberra.edu.au
// date:   22/01/18

data{
  int<lower=1> N;
  int<lower=0> N_censored;
  int<lower=1> S;
  int<lower=1> K;
  int<lower=1> N_plots;
  int<lower=1> N_sites;

  vector<lower=0>[S] y_observed[N];           //left censored data
  vector[K] X[N];                             //environmental covariates
  int<lower=1, upper=N_plots> plot[N];        //plot index
  int<lower=1, upper=N_sites> site[N_plots];  //site index
  int<lower=1> shape_prior;                   //LKJ prior shape parameter
}
parameters{
  matrix[S, K] B;                             //scaled species coefficients
  vector[K] mu_B;
  vector<lower=0>[K] sigma_B;

  vector[N_plots] B_plot;                     //scaled random effects
  real<lower=0> sigma_plot;

  vector[N_sites] B_site;
  real<lower=0> sigma_site;

  vector<lower=0>[S] sigma;                   //observation variances

  real<upper=0> y_censored[N_censored];       //censored values as parameters
}
transformed parameters{
  vector[S] mu[N];                            //linear predictor of cover

  // Generate linear predictors
  for(i in 1:N){
    for(j in 1:S){
      mu[i, j] = (B[j] * X[i]) + B_plot[plot[i]];
    }
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

      // Likelihood
      y_latent[i, j] ~ normal(mu[i, j], sigma[j]);
    }
  }

  // Priors
  // Coefficients have central means and variance.
  for(j in 1:K){
    B[ , j] ~ student_t(2, mu_B[j], sigma_B[j]);
  }

  mu_B ~ normal(0, 20);
  sigma_B ~ normal(0, 20);

  // Plots are nested within sites.
  B_plot ~ normal(B_site[site], sigma_plot);
  sigma_plot ~ cauchy(0, 3);

  B_site ~ normal(0, sigma_site);
  sigma_site ~ cauchy(0, 3);

  sigma ~ cauchy(0, 5);
}
generated quantities{
  vector[S] y_pred[N];                        //predicted fit

  for(i in 1:N){
    for(j in 1:S){
      y_pred[i, j] = normal_rng(mu[i, j], sigma[j]);
    }
  }
}
