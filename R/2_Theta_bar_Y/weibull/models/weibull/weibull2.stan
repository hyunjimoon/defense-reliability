data {
  int<lower=1> N; //653
  int<lower=1> N_engines; //5
  int<lower=1> N_ships; //99
  int<lower=1> N_ages;// 31
  int<lower=1> ship_engine_ind[N_ships];
  int<lower=1> N_missing[N_ships];
  real<lower=0> mean_failure[N_ships];
  real<lower=0,upper=31> mean_time[N_ships];
  real<lower=0,upper=1> y[31, 99];
  int<lower=0,upper=1> NA_ind[31, 99];
}



parameters {
  
    // weibull parameter
    // first layer
    real<lower=0> mu_alpha_bar;
    real<lower=0> mu_sigma_bar;

    // second layer
    real<lower=0> alpha_bar[N_engines];
    real<lower=0> sigma_bar[N_engines];
    real<lower=0> sd_alpha_bar;
    real<lower=0> sd_sigma_bar;
    
    //third layer
    real<lower=0> alpha[N_ships];
    real<lower=0> sigma[N_ships];
    real<lower=0> sd_alpha;
    real<lower=0> sd_sigma; 
    
    
    //failure, time parameter(for each engine)
    real<lower=0> mu_failure[N_engines];
    real<lower=0,upper=31> mu_time[N_engines];
    
    
}

model{
  
  //prior
  mu_alpha_bar ~ normal(0,2);
  mu_sigma_bar ~ normal(0,2);
  sd_alpha_bar ~ exponential(1);
  sd_sigma_bar ~ exponential(1);
  sd_alpha ~ gamma(10,10);
  sd_sigma ~ gamma(10,10);
  
  //layer
    for (e in 1:N_engines) {
        alpha_bar[e] ~ normal(mu_alpha_bar, sd_alpha_bar);
        sigma_bar[e] ~ normal(mu_sigma_bar, sd_sigma_bar);
    }
    
    for (n in 1:N_ships) {
        alpha[n] ~ normal(alpha_bar[ship_engine_ind[n]], sd_alpha);
        sigma[n] ~ normal(sigma_bar[ship_engine_ind[n]], sd_sigma);
    }
    
  
  //model
  for (n in 1:N_ships) {
    for (t in 1:N_ages){
      if (NA_ind[t,n]==0){
            target += weibull_lcdf( y[t,n] | alpha[n], sigma[n]);
    }
  }
  }
  
  //prior
  mu_failure ~ exponential(mean(mean_failure));
  mu_time ~ exponential(mean(mean_time));
  
    //model
    for (n in 1:N_ships) {
        mean_failure[n] ~ normal(mu_failure[ship_engine_ind[n]], 1);
        mean_time[n] ~ normal(mu_time[ship_engine_ind[n]], 1);
    }

}

generated quantities {
    matrix[N_ages,N_ships] y_new_pred;
    for (ship in 1:N_ships) {
      for (t in 1:N_ages) {
        y_new_pred[t,ship] = weibull_cdf(t,alpha[ship],sigma[ship]);
    }
}
}


