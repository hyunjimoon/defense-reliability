data {
  int<lower=0> N; // numbmer of observations
  int<lower=1>n_state; // number of states
  vector[n_state] state_obs[N];
  int time_obs[N];
  int max_allowed_state;
  int repair_state;
  int pm_state;
  int cm_state;
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
      if (i==j && i<= max_allowed_state){
        M[i, i] = 1;
      }
      else if (j==repair_state && i>max_allowed_state ){
        M[i, j] = 1;
      }
      else M[i, j] = 0;
    }
  }
}


parameters {
  real<lower=0> rate[n_state-1];
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
  D_rate_c[n_state,n_state] =0;
  
  D_pow[1] = scale_matrix_exp_multiply(1,  D_rate_c, D_init);
  DM_pow[1] = M * scale_matrix_exp_multiply(1,  D_rate_c, D_init);

  // M should be multiplied in rate
  for (i in 2:max(time_obs)){
    D_pow[i] = matrix_exp(D_rate_c)*D_pow[i-1];
    DM_pow[i] = M * matrix_exp(D_rate_c) * DM_pow[i-1];
  }
}

model {
  for(i in 1:N){
    if (time_obs[i] ==1){
      target += -(D_pow[time_obs[i]]*initial - state_obs[i])'*(D_pow[time_obs[i]]*initial - state_obs[i]); //how to prevent DM_pow[0]?
    }
    else{
      target += -(DM_pow[time_obs[i]-1] * D_pow[time_obs[i]] * initial - state_obs[i])'*(DM_pow[time_obs[i]-1] * D_pow[time_obs[i]] * initial - state_obs[i]);
    }
}
}
