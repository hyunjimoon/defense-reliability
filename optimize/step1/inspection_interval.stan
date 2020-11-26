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
  int cm_cost = 100;
  int pm_cost = 1;
  int interval_cnt = 10;
  int n_era = 2;
  row_vector[n_era] b_era = [20, max(time_obs)];
  //int change_t = 21; //[2];
  vector[n_state] initial = rep_vector(0, n_state);
  matrix[n_state, n_state] M;
  int cnt = 0;
  initial[initial_state] = 1;
  // for(i in 1:n_state){
  //   initial[i] = 0;
  //   if(i == initial_state){
  //     initial[i] = 1;
  //   }
  // }

  M = rep_matrix(0, n_state, n_state);
  for(i in 1:max_allowed_state){
    M[i, i] = 1;
  }
  for(i in (max_allowed_state+1):n_state){
    M[i, repair_state] = 1;
  }
}


parameters {
   //vector <lower = 1, upper = max(time_obs)> [interval_cnt] I_t; //  <lower = 1, upper = max(time_obs)> interval_cnt is in transformed data
   vector <lower =1, upper = 10>[max(time_obs)] interval;
}

transformed parameters {
  matrix[n_state, n_state] D_init = diag_matrix(rep_vector(1, n_state));
  matrix[n_state, n_state] t_p[max(time_obs)]; // transition probability at time t
  vector[n_state] state_t[max(time_obs)]; // state at time t

  for (i in 1:max(time_obs)){
    if (i ==1){
      t_p[i] = matrix_exp(D_rate); // matrix * bf_state -> state
    }
    else{
      for(j in 1:n_era){
        if(i <= b_era[j]){
          if(fmod(i, round(interval[j])) == 0){  // Not recommended, not good :(
            //t_p[i] = matrix_exp(D_rate) * t_p[i-1];
            t_p[i] = M *  matrix_exp(D_rate) * t_p[i-1];
          }else{
            t_p[i] = matrix_exp(D_rate) * t_p[i-1];
          }
          break;
        }
        else{
          continue;
        }
      }
    }
    state_t[i] = t_p[i]' * initial;
  }
}
model {
  for(i in 1:N)
    #target += - (rep_row_vector(pm_cost, cm_state-pm_state+1) * (state_t[time_obs[i]])[pm_state:cm_state] + rep_row_vector(cm_cost, n_state-cm_state) * (state_t[time_obs[i]])[cm_state+1:]);
    target += - (sum(pm_cost * (state_t[time_obs[i]])[pm_state:cm_state]) + sum(cm_cost * (state_t[time_obs[i]])[cm_state+1:]));
  }

generated quantities{
}
