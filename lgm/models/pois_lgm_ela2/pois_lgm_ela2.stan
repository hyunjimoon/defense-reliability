
functions {
  matrix K_functor (vector phi, vector[] x, real[] delta, int[] delta_int) {
    int n_ages = delta_int[1];
    int n_ships = delta_int[2];
    int s = delta_int[3];
    int ship_engine_ind[n_ships] = delta_int[4:(3+n_ships)];
    real ages[n_ages];
    for (t in 1:n_ages) ages[t] = t;
    vector [5] alpha_engine = phi[1:5];
    vector [5]rho_engine = phi[6:10];
    vector [99] alpha_ship = phi[11:109];
    vector [99]rho_ship = phi[110:208];
    matrix[n_ages, n_ages] K = cov_exp_quad(ages, alpha_engine[ship_engine_ind[s]], rho_engine[ship_engine_ind[s]]) + cov_exp_quad(ages, alpha_ship[s], rho_ship[s]) ;
    for (i in 1:n_ages) K[i, i] += 1e-8;
    return K; 
  }
}

data {
  int n_coordinates;
  int n_ages;
  int n_ships;
  int y [n_ages, n_ships];
  vector[n_ages] ye;
  int<lower=1> ship_engine_ind[n_ships];
  real alpha_mu_prior;
  real alpha_sd_prior;
  real rho_mu_prior;
  real rho_sd_prior;
}

transformed data {
  real tol = 1e-6;
  vector[n_ages] theta_0 = rep_vector(0, n_ages);
  int n_phi = 208;
  int n_samples[n_ages] = rep_array(1, n_ages);
  real delta[0];
  vector[n_ages] x[n_coordinates]; // arrays of vector
  for (t in 1:n_ages) x[1][t] = t;
}

parameters {
  vector<lower = 0> [5] alpha_engine;
  vector<lower = 0> [5] rho_engine;
  vector<lower = 0> [99] alpha_ship;
  vector<lower = 0> [99] rho_ship;
}

transformed parameters {
  vector[n_phi] phi;
  phi[1:5] = alpha_engine;
  phi[6:10] = rho_engine;
  phi[11:109] = alpha_ship;
  phi[110:208] =  rho_ship;
}

model {
  rho_engine ~ normal(rho_mu_prior, rho_sd_prior);
  alpha_engine ~ normal(alpha_mu_prior, alpha_sd_prior);
  rho_ship ~ normal(rho_mu_prior, rho_sd_prior);
  alpha_ship ~ normal(alpha_mu_prior, alpha_sd_prior);

  for (s in 1:n_ships){
    int y_s [n_ages] = y[,s];
    int delta_int[3+n_ships];
    delta_int[1] = n_ages;
    delta_int[2] = n_ships;
    delta_int[3] = s;
    delta_int[4:3+n_ships] = ship_engine_ind;
    target += laplace_marginal_poisson(y_s, n_samples, ye, K_functor,
                                   phi,x, delta, delta_int, theta_0);
  }
}

generated quantities {
  //matrix[n_ships, n_ages] theta;
  vector [n_ages] theta;
  matrix[n_ships, n_ages] y_pred;
  for (s in 1:n_ships){
    int y_s [n_ages] = y[,s];;
    int delta_int[3+n_ships];
    delta_int[1] = n_ages;
    delta_int[2] = n_ships;
    delta_int[3] = s;
    delta_int[4:3+n_ships] = ship_engine_ind;
    theta = laplace_approx_poisson_rng(y_s, n_samples, ye, K_functor,
                                 phi, x, delta, delta_int, theta_0);
    y_pred[s, ] = to_row_vector(poisson_log_rng(log(ye) + theta));                          
  }
}