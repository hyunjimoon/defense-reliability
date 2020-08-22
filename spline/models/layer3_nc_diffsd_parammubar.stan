data {
    int <lower = 1> K; // 11
    int <lower = 1> N; // number of total values(653)
    int <lower = 1> T; // time length of data(31)
    int <lower = 1> S; // number of ships (99)
    int <lower = 1> E; // number of unique engines (5)
    int <lower = 1> age[N];
    int <lower = 1> ship[N]; // ship type mapping
    int <lower = 1> engine[S]; //ship_ID to engine map
    matrix[T,K] B; // spline values, 2d array
    vector [N] Y; // actual values
    
    // below data are used for posterior predictive sampling
    int N_hat; // number of posterior prediction datapoints
    int age_hat[N_hat];
    int ship_hat[N_hat];
}

parameters {
    real mu_a_bar;
    vector[K] mu_w_bar;
    real<lower=0> s_a[S]; // model
    real<lower=0> s_w[S]; // model
    vector[S] s_a_tilde; // model
    
    vector[S] s_w_tilde; // model
    
    real s_a_bar;
    real s_w_bar;
    vector[E] s_a_bar_tilde; // model
    
    vector[E] s_w_bar_tilde; // model
    
    real<lower=0> s_Y; // model
}

transformed parameters { 
    vector[E] a_bar;
    vector[K] w_bar[E];
    vector[S] a;
    vector[K] w[S];
    vector[N] mu;
    
    for (e in 1:E){
        //a_bar[e] ~ normal(mu_a_bar, s_a_bar);
        a_bar[e] = mu_a_bar + s_a_bar * s_a_bar_tilde[e];
        //w_bar[e] ~ normal(mu_w_bar,s_w_bar);
        w_bar[e] = mu_w_bar + s_w_bar * s_w_bar_tilde[e];
    }
    for (s in 1:S){
        //a[s] ~ normal(a_bar[S2F[s]], s_a)
        a[s] = a_bar[engine[s]] + s_a[s] * s_a_tilde[s];
        w[s] = w_bar[engine[s]] + s_w[s] * s_w_tilde[s];
    }
    
  
    for (n in 1: N){
        mu[n] = a[ship[n]] + B[age[n]] * w[ship[n]];
    }
}

model {
    mu_a_bar ~ normal(0, 1);
    mu_w_bar ~ normal(0, 1);
    s_a_bar ~ gamma(5, 1);
    s_w_bar ~ gamma(5, 1);
    s_a_bar_tilde ~ normal(0,1);
    
    s_w_bar_tilde ~ normal(0,1);


    //s_a ~ gamma(4,1);
    s_a ~ exponential(0.5);
    //s_w ~ gamma(4, 1);
    s_w ~ exponential(0.5);
    s_a_tilde ~ normal(0,1);
    
    s_w_tilde ~ normal(0,1);

    s_Y ~ exponential(0.5);

    Y ~ normal(mu, s_Y);
}
generated quantities{
    vector[N] log_lik;
    for (i in 1:N) {
        log_lik[i] = normal_lpdf(Y[i]|mu[i], s_Y);
    }
    /*vector[S] ship_test;
    for(s in 1:S){
        ship_test[s] = s_a * s_a_tilde[engine[s]];
    }*/
    //real y_rep[N] = normal_rng(mu, s_Y);  # randomly sampled posterior predictive distribution

    real y_new_pred[N_hat]; // posterior predictive sampling
    for (i in 1:N_hat) {
        y_new_pred[i] = normal_rng(a[ship_hat[i]] + B[age_hat[i]] * w[ship_hat[i]], s_Y);
    }
}