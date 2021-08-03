data {
  int<lower=1> N; //28013
  int<lower=1> N_engines; //5
  int<lower=1> N_ships; //99
  int<lower=1> N_ages;// 31
  int<lower=1> ship_engine_ind[N_ships];
  int<lower=1> ship_n[N_ships];
  int<lower=1> ship_index[N];
  real<lower=0,upper=31> time[N];
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
    
    
}

model{
  
  //prior
  mu_alpha_bar ~ normal(0,2);
  mu_sigma_bar ~ normal(0,2);
  sd_alpha_bar ~ exponential(1);
  sd_sigma_bar ~ exponential(1);
  
  //layer
    for (e in 1:N_engines) {
        alpha_bar[e] ~ normal(mu_alpha_bar, sd_alpha_bar);
        sigma_bar[e] ~ normal(mu_sigma_bar, sd_sigma_bar);
    }
  
  //model
  for (i in 1:N){
    time[i] ~ weibull(alpha_bar[ship_engine_ind[ship_index[i]]],sigma_bar[ship_engine_ind[ship_index[i]]]);
  }
}


generated quantities {
    real<lower=0> time_pred[N];
    real<lower=0> pred;
    int ind=0;
    for (ship in 1:N_ships) {
      for (i in 1:ship_n[ship]){
        ind+=1;
        pred= round(weibull_rng(alpha_bar[ship_engine_ind[ship]],sigma_bar[ship_engine_ind[ship]]));
        time_pred[ind]=pred;
      }
    }
}


