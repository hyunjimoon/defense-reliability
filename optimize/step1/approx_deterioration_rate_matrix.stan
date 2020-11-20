data {
  int<lower=0> N; // numbmer of observations
  int<lower=1>n_state; // number of states
  vector[5] state_obs[N];
  int time_obs[N];
  int max_allowed_state;
  int repair_state;
  int<lower=0, upper=n_state> initial_state;
}

transformed data {
  vector[n_state] initial;
  matrix[n_state, n_state] M;

  for(i in 1:n_state){
    initial[i] = 0;
    if(i == initial_state){
      initial[i] = 1;
    }
  }

  for(i in 1:n_state){
    for(j in 1:n_state){
      M[i, j] = 0;
    }
  }
  for(i in 1:max_allowed_state){
    M[i, i] = 1;
  }
  for(i in (max_allowed_state+1):n_state){
    M[i, repair_state] = 1;
  }
}


parameters {
  real<lower=0> rate[n_state-1];
  //matrix[n_state,n_state] D_rate;
  // simplex[5] state_1;
  // simplex[4] state_2;
  // simplex[3] state_3;
  // simplex[2] state_4;
}

transformed parameters {
  matrix[n_state, n_state] D_init = diag_matrix(rep_vector(1, n_state));
  matrix[n_state, n_state] D_rate_c = rep_matrix(0, n_state, n_state);
  matrix[n_state, n_state] D_pow[max(time_obs)];
  matrix[n_state, n_state] DM_pow[max(time_obs)];
  for (i in 1:(n_state-1)){
    D_rate_c[i,i] = -rate[i];
    D_rate_c[i,i+1] = rate[i];
  }
  D_rate_c[n_state,n_state] =1;

  // M should be multiplied in rate
  for (i in 1:max(time_obs)){
    D_pow[i] = scale_matrix_exp_multiply(i,  D_rate_c, D_init);
    DM_pow[i] = scale_matrix_exp_multiply(i,  (D_rate_c * M), D_init);
  }
  
}

model {
  for(i in 1:N){
    if (time_obs[i] ==1){
      target += -log((D_pow[time_obs[i]]*initial) * state_obs[i]);
      } 
      else{
        target += -dot_product(log(DM_pow[time_obs[i]-1] * D_pow[time_obs[i]] * initial), state_obs[i]);
    }
}
}
