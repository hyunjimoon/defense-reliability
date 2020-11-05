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
  int<lower=0> N;
  vector[5] state_obs[N];
  int time_obs[N];
}

transformed data {
  vector[5] initial;
  initial[1] = 1;
  for(i in 2:5){
    initial[i] = 0;
  }
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  simplex[5] state_1;
  simplex[4] state_2;
  simplex[3] state_3;
  simplex[2] state_4;
}

transformed parameters {
  matrix[5, 5] D;
  matrix[5,5] D_pow[N];
  D[2, 1] = 0;
  D[3, 1] = 0;
  D[3, 2] = 0;
  D[4, 1] = 0;
  D[4, 2] = 0;
  D[4, 3] = 0;
  for(i in 1:5){
    D[1, i] = state_1[i];
  }
  for(i in 2:5){
    D[2, i] = state_2[i-1];
  }
  for(i in 3:5){
    D[3, i] = state_3[i-2];
  }
  for(i in 4:5){
    D[4, i] = state_4[i-3];
  }
  for(i in 1:4){
    D[5, i] = 0;
  }
  D[5, 5] = 1;
  
  
  for(i in 1:N){
    D_pow[i] = D;
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
    target += -(D_pow[i] * initial - state_obs[i]);
  }
}

