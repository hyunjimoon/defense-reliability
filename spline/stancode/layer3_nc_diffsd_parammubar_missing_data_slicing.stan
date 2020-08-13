data {
    int <lower = 1> K; // 11
    //int <lower = 1> N; // number of total values(653)

    int N_obs; // number of observed values
    int N_mis; // number of missing values
    int ii_obs[N_obs]; // index of observed values
    int ii_mis[N_mis]; // index of missing values
    vector [N_obs] Y_obs; // Y observed values

    int <lower = 1> T; // time length of data(31)
    int <lower = 1> S; // number of ships (99)
    int <lower = 1> E; // number of unique engines (5)
    int <lower = 1> age[N_obs + N_mis];
    int <lower = 1> ship[N_obs + N_mis]; // ship type mapping
    int <lower = 1> engine[S]; //ship_ID to engine map
    matrix[T,K] B; // spline values, 2d array

    
    // below data are used for posterior predictive sampling
    int N_hat; // number of posterior prediction datapoints
    int age_hat[N_hat];
    int ship_hat[N_hat];
}

transformed data{
    int<lower = 0> N = N_obs + N_mis;
}

parameters {
    real mu_a_bar;
    vector[K] mu_w_bar;
    real<lower=0> s_a; // model
    real<lower=0> s_w; // model
    vector[E] s_a_tilde; // model
    vector[E] s_w_tilde; // model
    real s_a_bar;
    real s_w_bar;
    vector[E] s_a_bar_tilde; // model
    vector[E] s_w_bar_tilde; // model
    real<lower=0> s_Y; // model
    
    vector[N_mis] Y_mis; // parameter to replace missing data(ii_mis)
    
}

transformed parameters { 
    vector[E] a_bar;
    vector[K] w_bar[E];
    vector[S] a;
    vector[K] w[S];
    vector[N] mu;

    vector[N] Y;  // combined data array
    Y[ii_obs] = Y_obs;
    Y[ii_mis] = Y_mis;

    for (e in 1:E){
        //a_bar[e] ~ normal(mu_a_bar, s_a_bar);
        a_bar[e] = mu_a_bar + s_a_bar * s_a_bar_tilde[e];
        //w_bar[e] ~ normal(mu_w_bar,s_w_bar);
        w_bar[e] = mu_w_bar + s_w_bar * s_w_bar_tilde[e];
    }
    for (s in 1:S){
        //a[s] ~ normal(a_bar[S2F[s]], s_a)
        a[s] = a_bar[engine[s]] + s_a * s_a_tilde[engine[s]];
        w[s] = w_bar[engine[s]] + s_w * s_w_tilde[engine[s]];
    }
    
  
    for (n in 1: N){
        mu[n] = a[ship[n]] + B[age[n]] * w[ship[n]];
    }
}

model {
    Y_mis ~ normal(0, 1);
    mu_a_bar ~ normal(0, 1);
    mu_w_bar ~ normal(0, 1);
    //s_a_bar ~ exponential(1);
    //s_w_bar ~ exponential(1);
    s_a_bar ~ gamma(4, 0.25);
    s_w_bar ~ gamma(4, 0.25);
    //s_a_bar ~ gamma(1, 1);
    //s_w_bar ~ gamma(1, 1);
    s_a_bar_tilde ~ normal(0,1);
    s_w_bar_tilde ~ normal(0,1);

    //s_a ~ gamma(10,10);
    //s_a ~ gamma(1, 1);
    s_a ~ gamma(4,0.25);
    //s_w ~ gamma(10,10);
    //s_w ~ gamma(1, 1);
    s_w ~ gamma(4, 0.25);
    s_a_tilde ~ normal(0,1);
    s_w_tilde ~ normal(0,1);

    s_Y ~ exponential(1);

    Y ~ normal(mu, s_Y);
}
generated quantities{
    /*vector[N] log_likelihood;
    for (i in 1:N) {
        log_likelihood[i] = normal_lpdf(Y[i]|mu[i], s_Y);
    }*/
    
    //real y_rep[N] = normal_rng(mu, s_Y);  // randomly sampled posterior predictive distribution

    real y_hat[N_hat]; // posterior predictive sampling
    for (i in 1:N_hat) {
        y_hat[i] = normal_rng(a[ship_hat[i]] + B[age_hat[i]] * w[ship_hat[i]], s_Y);
    }
}