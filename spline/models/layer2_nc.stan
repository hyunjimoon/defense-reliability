data {
    int <lower = 1> K; // number of basis functions
    int <lower = 1> N; // number of total values
    int <lower = 1> T; // time length of data 
    int <lower = 1> S; // number of obj  e.g ship
    int <lower = 1> E; // number of unique cat  e.g engine
    int <lower = 1> age[N];
    int <lower = 1> obj[N]; // obj type mapping
    int <lower = 1> cat[S]; //obj_ID to category map
    matrix[T,K] B; // spline values, 2d array
    vector [N] Y; // actual values
}

parameters {
    real mu_a_bar; // model
    vector[K] mu_w_bar; // model
    real<lower=0> s_a; // model
    real<lower=0> s_w; // model
    vector[S] s_a_tilde; // model
    
    vector[S] s_w_tilde; // model
    
    real<lower = 0> s_a_bar; // model
    real<lower = 0> s_w_bar; // model
    
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
        a[s] = a_bar[cat[s]] + s_a * s_a_tilde[s];
        w[s] = w_bar[cat[s]] + s_w * s_w_tilde[s];
    }
    
    
  
    for (n in 1: N){
        mu[n] = a[obj[n]] + B[age[n]] * w[obj[n]];
    }
}

model {

    mu_a_bar ~ normal(0, 1);
    mu_w_bar ~ normal(0, 1);

    s_a_bar ~ normal(0, 2);
    s_w_bar ~ gamma(2, 1.7);

    s_a_bar_tilde ~ normal(0, 1);
    s_w_bar_tilde ~ normal(0, 1);

    //s_a ~ exponential(0.5);
    s_a ~ gamma(5, 2);
    //s_w ~ exponential(0.5);
    s_w ~ gamma(5, 2);
    
    s_a_tilde ~ normal(0, 1); // normal(0,1)
    s_w_tilde ~ normal(0, 1);    

    s_Y ~ exponential(1); // exponential(4)
    
    Y ~ normal(mu, s_Y);
}

generated quantities{
    vector[N] log_lik;
    for (i in 1:N) {
        log_lik[i] = normal_lpdf(Y[i]|mu[i], s_Y);
    }
    real y_new_pred[N]; // posterior predictive sampling
    for (i in 1:N) {
        y_new_pred[i] = normal_rng(a[obj_hat[i]] + B[age_hat[i]] * w[obj_hat[i]], s_Y);
    }
}