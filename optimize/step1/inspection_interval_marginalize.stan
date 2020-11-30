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
  int cm_cost = 2;
  int pm_cost = 1;
  int inspection_cost = 7;
  int n_era = 2;
  int max_interval = 10;
  int b_era[n_era];
  vector[n_state] initial = rep_vector(0, n_state);
  matrix[n_state, n_state] M;
  initial[initial_state] = 1;

  b_era[1] = 20;
  b_era[2] = max(time_obs);

  M = rep_matrix(0, n_state, n_state);
  for(i in 1:max_allowed_state){
    M[i, i] = 1;
  }
  for(i in (max_allowed_state+1):n_state){
    M[i, repair_state] = 1;
  }
}

parameters {
}

transformed parameters {
}

model {
}

generated quantities{
  matrix[n_state, n_state] D_init = diag_matrix(rep_vector(1, n_state));
  matrix[n_state, n_state] t_p[max(time_obs)]; // transition probability at time t
  matrix[n_state, n_state] t_p_raw[max(time_obs)]; // transition probability without repairs at all
  vector[n_state] state_t[max(time_obs)]; // state at time t
  matrix[max_interval, max_interval] cost = rep_matrix(0, max_interval, max_interval);

  for(i1 in 1:max_interval){ // interval 1
    for(i2 in 1:max_interval){ // interval 2
      for (t in 1:max(time_obs)){
        if (t == 1){
          t_p[t] = matrix_exp(D_rate); // matrix * bf_state -> state
          t_p_raw[t] = matrix_exp(D_rate);
        }
        else{
          if(t <= b_era[1]){ // interval 1
            if(t % i1 == 0){
              t_p[t] = M * matrix_exp(D_rate) * t_p[t-1];
              cost[i1, i2] = cost[i1, i2] + inspection_cost;
            }
            else{
              t_p[t] = matrix_exp(D_rate) * t_p[t-1];
            }
          }
          else if(t > b_era[1] && t <= b_era[2]){ // interval 2
            if((t - b_era[1]) % i2 == 0){
              t_p[t] = M * matrix_exp(D_rate) * t_p[t-1];
              cost[i1, i2] = cost[i1, i2] + inspection_cost;
            }
            else{
              t_p[t] = matrix_exp(D_rate) * t_p[t-1];
            }
          }
          t_p_raw[t] = matrix_exp(D_rate) * t_p_raw[t-1];
        }
        state_t[t] = (t_p[t]' * initial) - (t_p_raw[t]' * initial);
        
      }
      for(i in 1:N){
        cost[i1, i2] = cost[i1, i2] + sum(state_t[time_obs[N]][pm_state:cm_state]) * pm_cost + sum(state_t[time_obs[N]][cm_state:]) * cm_cost;
      }
    }
  }
}
