//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // numbmer of observations
  int<lower=1>n_states; // number of states
  vector[5] state_obs[N];
  int time_obs[N];
  int max_allowed_state;
  int repair_state;
  int<lower=1, upper=n_states> initial_state;
}

transformed data {
  vector[n_states] initial;
  matrix[n_states, n_states] M;

  for(i in 1:n_states){
    initial[i] = 0;
    if(i == initial_state){
      initial[i] = 1;
    }
  }

  for(i in 1:n_states){
    for(j in 1:n_states){
      M[i, j] = 0;
    }
  }
  for(i in 1:max_allowed_state){
    M[i, i] = 1;
  }
  for(i in max_allowed_state+1:n_states){
    M[i, repair_state] = 1;
  }
}


parameters {
  simplex[n_states] probs[n_states];
}

transformed parameters {
  matrix[n_states, n_states] D;
  matrix[n_states, n_states] D_pow[N];
  for(i in 1:n_states){
   for(j in 1:n_states){
     D[i,j] = probs[i, j];
   }
  }
  for(i in 1:N){
    D_pow[i] = D * M;
    for(j in 2:time_obs[i]){
      D_pow[i] = D_pow[i] * D_pow[i];
    }
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  for(i in 1:N){
    target += -(initial'  * D_pow[i] * D - state_obs[i]');
  }
}

