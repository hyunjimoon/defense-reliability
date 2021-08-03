
data {
  int n_obs;
  int n_covariates;
  int y[n_obs];
  vector[n_obs] ye;
  vector[n_covariates] x[n_obs];
  real alpha_mu_prior;
  real alpha_sd_prior;
  real rho_mu_prior;
  real rho_sd_prior;
}

transformed data {
  real delta = 1e-8;
}

parameters {
  real<lower = 0> alpha;
  real<lower = 0> rho;
  vector[n_obs] eta;
}

transformed parameters {
   vector[n_obs] theta;
   {
     matrix[n_obs, n_obs] L_Sigma;
     matrix[n_obs, n_obs] Sigma;
     Sigma = cov_exp_quad(x, alpha, rho);
     for (n in 1:n_obs) Sigma[n, n] = Sigma[n,n] + delta;
     L_Sigma = cholesky_decompose(Sigma);
     theta = L_Sigma * eta;
   }
}

model {
  rho ~ normal(rho_mu_prior, rho_sd_prior);
  alpha ~ normal(alpha_mu_prior, alpha_sd_prior);

  eta ~ normal(0, 1);
  y ~ poisson_log(log(ye) + theta);
}

generated quantities {
  int y_pred[n_obs] = poisson_log_rng(log(ye) + theta);
}