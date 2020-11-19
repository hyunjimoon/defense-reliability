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
  matrix[n_state, n_state] D_rate;
}

transformed data {
  int cm_cost = 1;
  int pm_cost = 2;
  //int change_t = 21; //[2];
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
   vector <lower = 1, upper = max(time_obs)>[3] I_t;
}

transformed parameters {
  matrix[n_state, n_state] D_init = diag_matrix(rep_vector(1, n_state));
  matrix[n_state, n_state] D_pow[max(time_obs)];
  matrix[n_state, n_state] DM_pow[max(time_obs)];
  vector[n_state] t_p[max(time_obs)];

  for (i in 1:max(time_obs)){
    D_pow[i] = scale_matrix_exp_multiply(i, D_rate, D_init); // state prob matrix after time i given chain started at time 0
    if (i ==1){
      t_p[i] = (D_pow[i])' * initial; // need to extract entire row not column, so we transpose
      DM_pow[i] = D_init * matrix_exp(D_rate);
    }else{
      if (i == I_t[1]|| i == I_t[2]|| i == I_t[3]){ // !!!! ALWAYS FALSE(comparison of integer to a real value)
        DM_pow[i] = DM_pow[i-1] * matrix_exp(D_rate) * M;
      }else{
        DM_pow[i] = DM_pow[i-1] * matrix_exp(D_rate);
      }
    t_p[i] = (DM_pow[i])' * initial;
    }
  }
}

model {
  for(i in 1:N){
    target += pm_cost * (t_p[time_obs[i]])[pm_state] + cm_cost * (t_p[time_obs[i]])[cm_state];
  }
}

generated quantities{
}
