data {
  int<lower=1> N; //653
  int<lower=1> N_engines; //5
  int<lower=1> N_ships; //99
  int<lower=1> N_ages_obs; // 31
  int<lower=1> N_ages;// 31
  int<lower=1> ship_engine_ind[N_ships];
  int<lower=1,upper=99> ship_ind[N];
  int<lower=1> age_ind[N]; 
  vector[N] y; //
}
transformed data {
  real ages[N_ages];
  vector[8] counts;
  for (t in 1:N_ages)
    ages[t] = t;
  for (i in 1:8)
    counts[i] = 2;
}
parameters {
  matrix[N_ages,N_engines] GP_engine_std;
  matrix[N_ages,N_ships] GP_ship_std;
  vector[N_ages_obs] age_std;
  vector[N_ships] ship_std;
  vector[N_engines] engine_std;
  real<lower=0> tot_var;
  simplex[8] prop_var;
  real mu;
  real<lower=0> length_GP_engine_long;
  real<lower=0> length_GP_ship_long;
  real<lower=0> length_GP_engine_short;
  real<lower=0> length_GP_ship_short;
}
transformed parameters {
  matrix[N_ages,N_engines] GP_engine;
  matrix[N_ages,N_ships] GP_ship;

  vector[N_ages_obs] age_re;
  vector[N_ships] ship_re;
  vector[N_engines] engine_re;
  vector[8] vars;

  real sigma_age;
  real sigma_engine;
  real sigma_ship; // vector[99] sigma_ship;

  real sigma_error_ship_2;

  real sigma_GP_engine_long;
  real sigma_GP_ship_long;
  real sigma_GP_engine_short;
  real sigma_GP_ship_short;

  vars = 8 * prop_var * tot_var;
  sigma_age = sqrt(vars[1]);
  sigma_engine = sqrt(vars[2]);
  sigma_ship = sqrt(vars[3]); // for 1:99
  sigma_GP_engine_long = sqrt(vars[4]);
  sigma_GP_ship_long = sqrt(vars[5]);
  sigma_GP_engine_short = sqrt(vars[6]);
  sigma_GP_ship_short = sqrt(vars[7]);
  sigma_error_ship_2 = sqrt(vars[8]);

  engine_re = sigma_engine * engine_std;
  age_re = sigma_age * age_std;
  ship_re = sigma_ship * ship_std; //ship_re = sigma_ship[ship_engine_ind] .* ship_std;
  
  {
    matrix[N_ages, N_ages] cov_engine; 
    matrix[N_ages, N_ages] cov_ship; 
    matrix[N_ages, N_ages] L_cov_engine; 
    matrix[N_ages, N_ages] L_cov_ship; 

    cov_engine = cov_exp_quad(ages, sigma_GP_engine_long, 
                                  length_GP_engine_long)
               + cov_exp_quad(ages, sigma_GP_engine_short, 
                                  length_GP_engine_short);
    cov_ship = cov_exp_quad(ages, sigma_GP_ship_long, 
                                  length_GP_ship_long)
               + cov_exp_quad(ages, sigma_GP_ship_short, 
                                  length_GP_ship_short);
    for (age in 1:N_ages) {
      cov_engine[age, age] = cov_engine[age, age] + 1e-6;
      cov_ship[age, age] = cov_ship[age, age] + 1e-6;
    }

    L_cov_engine = cholesky_decompose(cov_engine);
    L_cov_ship = cholesky_decompose(cov_ship);
    GP_engine = L_cov_engine * GP_engine_std;
    GP_ship = L_cov_ship * GP_ship_std;
  }
}
model {
  vector[N] obs_mu;

  for (n in 1:N) {
    obs_mu[n] = mu + age_re[age_ind[n]] 
              + ship_re[ship_ind[n]] 
              + engine_re[ship_engine_ind[ship_ind[n]]]
              + GP_engine[age_ind[n],ship_engine_ind[ship_ind[n]]]
              + GP_ship[age_ind[n],ship_ind[n]];
  }
  y ~ normal(obs_mu, sigma_error_ship_2); 

  to_vector(GP_engine_std) ~ normal(0, 1);
  to_vector(GP_ship_std) ~ normal(0, 1);
  age_std ~ normal(0, 1);
  ship_std ~ normal(0, 1);
  engine_std ~ normal(0, 1);
  mu ~ normal(.5, .5);
  tot_var ~ gamma(3, 3);
  prop_var ~ dirichlet(counts); // proper?
  length_GP_engine_long ~ weibull(30,8);
  length_GP_ship_long ~ weibull(30,8);
  length_GP_engine_short ~ weibull(30,3);
  length_GP_ship_short ~ weibull(30,3);
}
generated quantities {
  matrix[N_ages,N_ships] y_new;
  matrix[N_ages,N_ships] y_new_pred;

  {
    for (ship in 1:N_ships) {
      for (t in 1:N_ages) {
          y_new[t, ship] = ship_re[ship] 
                         + engine_re[ship_engine_ind[ship]]
                         + GP_ship[t,ship]
                         + GP_engine[t,ship_engine_ind[ship]]
                         + (mu + age_re[t]);
        
        y_new_pred[t,ship] = normal_rng(y_new[t,ship],
                                         sigma_error_ship_2);
      }
    }
  }
}