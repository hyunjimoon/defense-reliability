//order E > S -- Y
data {
    int <lower = 1> K; // number of basis functions
    int <lower = 1> N; // number of total values
    int <lower = 1> T; // time length of data
    int <lower = 1> S; // number of x_atom  e.g ship
    int <lower = 1> E; // number of unique cat  e.g engine
    int <lower = 1> age[N];
    int <lower = 1> x[N]; // x type mapping (id to seriaal #)
    int <lower = 1> cat[S]; //x_ID to category map
    matrix[T,K] B; // spline values, 2d array
    vector [N] Y; // actual values
}

parameters {
    //first layer - meta cat
    real<lower = 0> s_a_bar;
    real<lower = 0> s_w_bar;
    vector[E] eta_a_bar;
    vector[E] eta_w_bar;
    real mu_a_bar;
    vector[K] mu_w_bar;
    //second layer- per cat
    real<lower=0> s_a;
    real<lower=0> s_w;
    vector[S] eta_a;
    vector[S] eta_w;

    real<lower=0> s_Y;
}

transformed parameters {
    vector[E] a_bar;
    vector[K] w_bar[E];
    vector[S] a;
    vector[K] w[S];
    vector[N] mu;

    for (e in 1:E){
        //a_bar[e] ~ normal(mu_a_bar, s_a_bar);
        a_bar[e] = mu_a_bar + s_a_bar * eta_a_bar[e];
        //w_bar[e] ~ normal(mu_w_bar,s_w_bar);
        w_bar[e] = mu_w_bar + s_w_bar * eta_w_bar[e];
    }
    for (s in 1:S){
        //a[s] ~ normal(a_bar[S2F[s]], s_a)
        a[s] = a_bar[cat[s]] + s_a * eta_a[s];
        w[s] = w_bar[cat[s]] + s_w * eta_w[s];
    }
    for (n in 1: N){
        mu[n] = a[x[n]] + B[age[n]] * w[x[n]];
    }
}

model {
    mu_a_bar ~ normal(0, 1);
    mu_w_bar ~ normal(0, 1);
    s_a_bar ~ normal(0, 2);
    s_w_bar ~ gamma(2, 1.7);

    eta_a_bar ~ normal(0, 1);
    eta_w_bar ~ normal(0, 1);

    //s_w, s_a ~ exponential(0.5);
    s_a ~ gamma(5, 2);
    s_w ~ gamma(5, 2);
    eta_a ~ normal(0, 1);
    eta_w ~ normal(0, 1);
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
            y_new_pred[i] = normal_rng(mu, s_Y);
    }
}
