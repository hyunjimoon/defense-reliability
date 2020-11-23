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
  int interval_cnt = 10;
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
      M[i, j] = D_rate[i, j];
    }
  }
  for(i in (max_allowed_state+1):n_state){
    M[i, min(n_state, i+1)] = 0;
    M[i, repair_state] = fabs(M[i, i]);
  }
}

parameters {
   vector <lower = 1, upper = max(time_obs)>[interval_cnt] I_t; // interval_cnt is in transformed data
}


transformed parameters {
  matrix[n_state, n_state] D_init = diag_matrix(rep_vector(1, n_state));
  matrix[n_state, n_state] t_p[max(time_obs)]; // transition probability at time t
  vector[n_state] state_t[max(time_obs)]; // state at time t
  for (i in 1:max(time_obs)){
    if (i ==1){
      t_p[i] = matrix_exp(D_rate); // matrix * bf_state -> state
    }else{
      for(j in 1:interval_cnt){
        if(i <= I_t[j] && I_t[j] < (i+1)){
          t_p[i] =  matrix_exp(M) * t_p[i-1];
          break;
        }
        else if(j == interval_cnt){
          t_p[i] = matrix_exp(D_rate) * t_p[i-1];
        }
      }
    }
    state_t[i] = t_p[i]' * initial;  // Need to transpose t_p before multiplication! initial is a column vector.
  }
}

model {
  for(i in 1:N){
    target += -(pm_cost * (state_t[time_obs[i]])[pm_state] + cm_cost * (state_t[time_obs[i]])[cm_state]);
  }
}

generated quantities{
}

