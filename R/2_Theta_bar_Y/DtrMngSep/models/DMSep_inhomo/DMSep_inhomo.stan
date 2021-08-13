data {
  int<lower=0> N; // numbmer of obs
  int<lower=0> T; // number of time bins = max age
  int<lower=1>S; // number of states
  int<lower=1> P; // number of periods for inhomogenous dtr.
  int<lower=0, upper =S> states[N]; //observed state
  int obs2time[N]; // map obs to age
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
  matrix[S, S] Dtr[P]; // Deterioration
  matrix[S, S] Dtr_pow[T]; // Deterioration
  matrix[S, S] Mnt; // Maintenance
  simplex[S] latent_states[T];
  real<lower = 0> tmp_p11;
  real<lower = 0> tmp_p21;
  real<lower = 0> tmp_p31;
    // Maintenance
    // is left stoch. matrix, transpose of eq.14 from the paper
    Mnt = [[1, p21, p31],
           [0, 1-p21, 1-p31],
           [0, 0, 0]];
  // Deterioration by period
  // is left stoch. matrix, transpose of eq.13 from the paper
    for(p in 1:P){
      tmp_p11 = exp(-rate[p,1]- rate[p,2]);
      tmp_p21 = rate[p,1] * exp(-rate[p,3]) * (1-exp(-(rate[p,1]+ rate[p,2] - rate[p,3]))) / (rate[p,1]+ rate[p,2] - rate[p,3]);
      tmp_p31 = exp(-rate[p,3]);
      Dtr[p] = [[tmp_p11, 0, 0],
                [tmp_p21, tmp_p31, 0],
                [1 - tmp_p11 - tmp_p21, 1 - tmp_p31, 1]];
  }
    // (In) or homogenous Dtr
    latent_states[1] = Dtr[1] * initial;
    for (t in 2:T){
      if (t <= 8){latent_states[t] = Dtr[1] * Mnt * latent_states[t-1];}
      else if (t <=20){latent_states[t] = Dtr[2] * Mnt * latent_states[t-1];}
      else if (t <=26){latent_states[t] = Dtr[3] * Mnt * latent_states[t-1];}
      else{latent_states[t] = Dtr[4] * Mnt * latent_states[t-1];}
    }
}
model {
  for (n in 1:N){
    states[n] ~ categorical(latent_states[obs2time[n]]);
  }
}
generated quantities{
int<lower=0, upper =S> gen_states[N];
  for (n in 1:N){
    gen_states[n] = categorical_rng(latent_states[obs2time[n]]);
  }
}
