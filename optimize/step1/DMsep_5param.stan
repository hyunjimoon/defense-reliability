functions {
  matrix generate_p_t(int n_state, real[] rate, real t){
    matrix[n_state, n_state] return_matrix = rep_matrix(0, n_state, n_state);

    return_matrix[1, 1] = exp(-(rate[1] + rate[2]) * t);
    return_matrix[1, 2] = rate[1] * exp(-rate[3] * t) * (1 - exp(-(rate[1] + rate[2] - rate[3]) * t)) / (rate[1] + rate[2] - rate[3]);
    return_matrix[1, 3] = 1 - return_matrix[1, 1] - return_matrix[1, 2];

    return_matrix[2, 2] = exp(-rate[3] * t);
    return_matrix[2, 3] = 1 - return_matrix[2, 2];

    return_matrix[3, 3] = 1;
    return(return_matrix);
  }
}

data {
  int<lower=0> N; // numbmer of observations
  int<lower=1>n_state; // number of states
  vector[n_state] state_obs[N];
  int time_obs[N];
  int<lower=0, upper=n_state> initial_state;
}

transformed data {
  vector[n_state] initial;

  for(i in 1:n_state){
    initial[i] = 0;
    if(i == initial_state){
      initial[i] = 1;
    }
  }
}

parameters {
  real<lower=0> rate[3];
  real<lower=0, upper=1> p;
}

transformed parameters {
  matrix[n_state, n_state] D_pow[max(time_obs)];
  matrix[n_state, n_state] DM_pow[max(time_obs)];
  matrix[n_state, n_state] TPM[max(time_obs)];
  matrix[n_state, n_state] M;
  
  
  for(i in 1:n_state){
    for(j in 1:n_state){
      if ((i==j && i< 3)){
        M[i, i] = 1;
      }
      else if (i==3 && j ==1){
        M[i, j] = p;
      }
      else if (i==3 && j ==2){
        M[i, j] = 1-p;
      }
      else M[i, j] = 0;
    }
  }
  
  // TPM[1,1] = exp(-(rate[1]+ rate[2]));
  // TPM[1,2] = rate[1] * exp(-rate[3]) * (1-exp(-(rate[1]+ rate[2] - rate[3]))) / (rate[1]+ rate[2] - rate[3]);
  // TPM[1,3] = 1 - TPM[1,1] - TPM[1,2];
  // TPM[2,1] = 0;
  // TPM[2,2] = exp(-rate[3]);
  // TPM[2,3] = 1 - TPM[2,2];
  // TPM[3,1] = 0;
  // TPM[3,2] = 0;
  // TPM[3,3] = exp(0);


  D_pow[1] = generate_p_t(n_state, rate, 1);
  TPM[1] = D_pow[1];
  DM_pow[1] = M * D_pow[1];//TPM;
  // M should be multiplied in rate
  for (i in 2:max(time_obs)){
    TPM[i] = generate_p_t(n_state, rate, i);
    D_pow[i] = TPM[i] * D_pow[i-1];
    DM_pow[i] = M * D_pow[i];
  }
}


model {
  for(i in 1:N){
    if (time_obs[i] ==1){
      target += -(D_pow[time_obs[i]]*initial - state_obs[i])'*(D_pow[time_obs[i]]*initial - state_obs[i]); //how to prevent DM_pow[0]?
    }
    else{
      target += -(DM_pow[time_obs[i]-1]  * initial - state_obs[i])'*(DM_pow[time_obs[i]-1] * initial - state_obs[i]);
    }
  }
}
