data {
  int<lower=0> N; // numbmer of observations
  int<lower=1>n_state; // number of states
  vector[n_state] state_obs[N];
  int time_obs[N];
  int pm_repair;
  int cm_repair;
  int pm_init;
  int cm_init;
  int<lower=0, upper=n_state> initial_state;
  matrix[n_state, n_state] D_rate;
}

transformed data {
  int cm_cost = 40;
  int pm_cost = 10;
  int inspection_cost = 1;
  int n_era = 2;
  int max_interval = 10;
  int b_era[n_era];  // upper bound for each interval
  vector[n_state] initial = rep_vector(0, n_state);
  matrix[n_state, n_state] M;
  initial[initial_state] = 1;

  b_era[1] = 20;
  b_era[2] = max(time_obs);

  M = rep_matrix(0, n_state, n_state);
  for(i in 1:(pm_init - 1)){
    M[i, i] = 1;
  }
  for(i in pm_init:(cm_init -1)){
    M[i, pm_repair] = 1;
  }
  for(i in cm_init:n_state){
    M[i, cm_repair] = 1;
  }
  for(r in 1:n_state){  // Verify Maintenance matrix for errors
    if(sum(M[r, 1:n_state]) != 1){
      reject("A row of maintenance failed rowSum() == 1 check. row number=", r);
    }
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
  // matrix[n_state, n_state] M_test;
  // for(r in 1:n_state){
  //   for(c in 1:n_state){
  //     M_test[r, c] = M[r, c];
  //   }
  // }
  for(i1 in 1:max_interval){ // interval 1
    for(i2 in 1:max_interval){// interval 2
      if (b_era[1]/i1+(b_era[2]-b_era[1])/i2 < 10) {
        cost[i1,i2] = positive_infinity();
        continue;
      }
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
        cost[i1, i2] = cost[i1, i2] + sum(state_t[time_obs[N]][pm_init:(cm_init - 1)]) * pm_cost + sum(state_t[time_obs[N]][cm_init:]) * cm_cost;
      }
    }
  }
}
