data {
  int<lower=0> N; // numbmer of obs
  int<lower=0> T; // number of time bins = max age
  int<lower=1>S; // number of states
  int<lower=1> P; // number of periods for inhomogenous dtr.
  int<lower=0, upper =S> states[N]; //observed state
  int obs2age[N]; // map obs to age
  int<lower=0, upper=S> initial_state;
}

transformed data {
  vector[S] initial;
  vector<lower=0>[S] alpha;
  vector<lower=0>[S] beta;
  initial = rep_vector(0, S);
  initial[initial_state] = 1;
  alpha = rep_vector(3,3);
  beta = rep_vector(3,3);
}

parameters {
  real<lower=0> rate[P, S];
  real<lower=0, upper=1> p21;
  real<lower=0, upper=1> p31;
  //real<lower=0, upper=1> r;
  //simplex[K] theta[K];        // transit probs
}

transformed parameters {
  simplex[S] Dtr[P,S]; // Deterioration
  simplex[S] Mnt[S]; // Maintenance
  simplex[S] DM_pow[T, S];
  matrix[S, S] tmp_p;
  // Deterioration by period
  // is left stoch. matrix, transpose of eq.13 from the paper
  for(p in 1:P){
    tmp_p[1,1] = exp(-rate[p,1]- rate[p,2]);
    tmp_p[2,1] = rate[p,1] * exp(-rate[p,3]) * (1-exp(-(rate[p,1]+ rate[p,2] - rate[p,3]))) / (rate[p,1]+ rate[p,2] - rate[p,3]);
    tmp_p[3,1] = exp(-rate[p,3]);
    Dtr[p,1] = to_vector([tmp_p[1,1], tmp_p[2,1],1 - tmp_p[1,1] - tmp_p[2,1]]);
    Dtr[p,2] = to_vector([0, tmp_p[3,1], 1 - tmp_p[3,1]]);
    Dtr[p,3] = to_vector([0, 0,1]);
  }
  // Inhomogenous Dtr
  DM_pow[1] = Dtr[1];
  for (t in 2:T){
    if (t <= 8){DM_pow[t] = Dtr[1] * Mnt * DM_pow[t-1];}
    else if (t <=20){DM_pow[t] = Dtr[2] * Mnt * DM_pow[t-1];}
    else if (i <=26){DM_pow[t] = Dtr[3] * Mnt * DM_pow[t-1];}
    else{DM_pow[t] = Dtr[4] * M * DM_pow[t-1];}
  }
  // Maintenance
  // is left stoch. matrix, transpose of eq.14 from the paper
  Mnt[1] = to_vector([1, 0, 0]);
  Mnt[2] = to_vector([p21,1- p21, 0]);
  Mnt[3] = to_vector([p31, 1- p31, 0]);
}
model {
  for (s in 1:S){
    Mnt[s] ~  dirichlet(alpha);
    for (p in 1:P){
      Dtr[p, s] ~ dirichlet(beta);
    }
  }
  for (t in 1:T){
    states[t] ~ categorical(DM_pow[t, states[t-1]]);
  }
}
