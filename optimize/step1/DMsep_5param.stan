
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
  real<lower=0> rate[4,3];
  real<lower=0, upper=1> p;
}

transformed parameters {
  matrix[n_state, n_state] DM_pow[31];
  matrix[n_state, n_state] D[4];
  matrix[n_state, n_state] M;
  
  
  for(i in 1:n_state){
    for(j in 1:n_state){
      if ((i==j && i< 3)){
        M[i, i] = 1;
      }
      else if (i==1 && j ==3){
        M[i, j] = p;
      }
      else if (i==2 && j ==3){
        M[i, j] = 1-p;
      }
      else M[i, j] = 0;
    }
  }
  
  for(i in 1:4){
    D[i][1,1] = exp(-(rate[i,1]+ rate[i,2]));
    D[i][2,1] = rate[i,1] * exp(-rate[i,3]) * (1-exp(-(rate[i,1]+ rate[i,2] - rate[i,3]))) / (rate[i,1]+ rate[i,2] - rate[i,3]);
    D[i][3,1] = 1 - D[i][1,1] - D[i][1,2];
    D[i][1,2] = 0;
    D[i][2,2] = exp(-rate[i,3]);
    D[i][3,2] = 1 - D[i][2,2];
    D[i][1,3] = 0;
    D[i][2,3] = 0;
    D[i][3,3] = exp(0);
  }
  
  DM_pow[1] = D[1];

  for (i in 2:max(time_obs)){
    if (i <= 8){
      DM_pow[i] = D[1] * M * DM_pow[i-1];
    }
    else if (i <=20){
      DM_pow[i] = D[2] * M * DM_pow[i-1];
    }
    else if (i <=26){
      DM_pow[i] = D[3] * M * DM_pow[i-1];
    }
    else{
      DM_pow[i] = D[4] * M * DM_pow[i-1];
    }
  }
}


model {
  for(i in 1:N){
    if (time_obs[i] ==1){
      target += -(D[1]*initial - state_obs[i])'*(D[1]*initial - state_obs[i]); //how to prevent DM_pow[0]?
    }
    else{
      target += -(DM_pow[time_obs[i]]  * initial - state_obs[i])'*(DM_pow[time_obs[i]] * initial - state_obs[i]);
    }
  }
}
