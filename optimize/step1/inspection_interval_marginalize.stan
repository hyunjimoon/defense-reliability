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
  int max_interval = 8;
  int b_era[n_era];  // upper bound for each interval
  vector[n_state] state_0 = rep_vector(0, n_state);
  matrix[n_state, n_state] M;
  state_0[initial_state] = 1;

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
  matrix[n_state, n_state] t_p[max(time_obs)]; // transition probability at time t
  vector[n_state] state_t[max(time_obs)]; // state at time t
  vector[n_state] state_t_NM;
  matrix[max_interval, max_interval] cost[max_interval, max_interval]; //check init to 0
  real M_cost;
  real NM_cost;
  //matrix[max_interval,max_interval,max_interval,max_interval]cost; : rep_matrix(0, max_interval, max_interval),max_interval,max_interval);
  for (a in 1:max_interval){
    for (b in 1:max_interval){
     cost[a,b] = rep_matrix(0, max_interval, max_interval);
    }
  }
  for(i1 in 1:max_interval){ // interval 1
    for(i2 in 1:max_interval){// interval 2
      for(i3 in 1:max_interval){ // interval 1
        for(i4 in 1:max_interval){// interval 2
          if (i1 + i2 + i3 + i4 > 20) { // TODO could be improved
            cost[i1,i2,i3,i4] = positive_infinity();
            continue;
          }
          M_cost = 0;
          NM_cost = 0;
          for (t in 1:max(time_obs)){
            if (t == 1){
              t_p[t] = matrix_exp(D_rate); // matrix * bf_state -> state
              state_t[t] = t_p[t]' * state_0; // matrix * bf_state -> state #t_p[t];
            }
            
            else{
              t_p[t] = matrix_exp(D_rate) * t_p[t-1];
              state_t[t] = t_p[t]' * state_t[t-1];
              state_t_NM = t_p[t]' * state_t[t-1];
              if ((t%9==0)||(t % i1 == 0 && t <= 7) || ((t-7) % i2 == 0 && t > 7 && t <= 20) || ((t-20) % i3 == 0 &&  t > 20 && t <= 27) || ((t-27) % i4 == 0 && t >27)){
                state_t[t] = M' * state_t[t];
                M_cost = M_cost + inspection_cost + pm_cost * sum(state_t[t][pm_init:(cm_init - 1)]) + cm_cost * sum(state_t[t][cm_init:]);
                NM_cost = NM_cost + inspection_cost + pm_cost * sum(state_t_NM[pm_init:(cm_init - 1)]) + cm_cost * sum(state_t_NM[cm_init:]);
                cost[i1, i2, i3, i4] = cost[i1, i2, i3, i4] + NM_cost - M_cost; // NM - M vs M - NM?
              }
            }
          }
        }
      }
    }
  }
}
