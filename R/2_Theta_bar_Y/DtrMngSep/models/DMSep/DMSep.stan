data {
  int<lower=0> N; // numbmer of obs
  int<lower=0> T; // number of time bins = max age
  int<lower=1>S; // number of states
  int<lower=1> P; // number of periods for inhomogenous dtr.
  int<lower=0, upper =S> states[N]; //observed state
  int obs2time[N]; // map obs to time
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
  real<lower=0> rate[S];
  real<lower=0, upper=1> p21;
  real<lower=0, upper=1> p31;
}

transformed parameters {
  matrix[S, S] Dtr; // Deterioration
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
  tmp_p11 = exp(-rate[1]- rate[2]);
  tmp_p21 = rate[1] * exp(-rate[3]) * (1-exp(-(rate[1]+ rate[2] - rate[3]))) / (rate[1]+ rate[2] - rate[3]);
  tmp_p31 = exp(-rate[3]);
  Dtr = [[tmp_p11, 0, 0],
            [tmp_p21, tmp_p31, 0],
            [1 - tmp_p11 - tmp_p21, 1 - tmp_p31, 1]];
  // (In) or homogenous Dtr
  latent_states[1] = Dtr * initial;
  for (t in 2:T){
    latent_states[t] =  (Dtr * Mnt) *latent_states[t-1]; //matrix_power((Dtr * Mnt), (t-1)) * initial;
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
