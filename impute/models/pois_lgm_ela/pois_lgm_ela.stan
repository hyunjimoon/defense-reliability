
functions {
  matrix K_functor (vector phi,
                    vector[] x,
                    real[] delta, int[] delta_int) {
    int n_obs = delta_int[1];
    matrix[n_obs, n_obs] K = cov_exp_quad(x, phi[1], phi[2]);
    for (i in 1:n_obs) K[i, i] += 1e-8;
    return K;
  }
}

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
  real tol = 1e-6;
  vector[n_obs] theta_0 = rep_vector(0, n_obs);
  int n_phi = 2;
  int n_samples[n_obs] = rep_array(1, n_obs);
  real delta[0];
  int delta_int[1] = {n_obs};
}

parameters {
  real<lower = 0> alpha;
  real<lower = 0> rho;
}

transformed parameters {
  vector[n_phi] phi;
  phi[1] = alpha;
  phi[2] = rho; 
}

model {
  rho ~ normal(rho_mu_prior, rho_sd_prior);
  alpha ~ normal(alpha_mu_prior, alpha_sd_prior);

  target += laplace_marginal_poisson(y, n_samples, ye, K_functor,
                                     phi, x, delta, delta_int, theta_0);
}

generated quantities {
  vector[n_obs] theta
    = laplace_approx_poisson_rng(y, n_samples, ye, K_functor,
                                 phi, x, delta, delta_int, theta_0);
  int y_pred[n_obs] = poisson_log_rng(log(ye) + theta);
}
