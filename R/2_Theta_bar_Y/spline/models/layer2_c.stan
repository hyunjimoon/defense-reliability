
data {
    int<lower=1> K; // number of knots
    int<lower=1> N; // number of datapoints
    int<lower=1> T; // maximum age
    int<lower=1> S; // number of objs
    int<lower=1> E; // number of cats
    int<lower=1> age[N];
    int<lower=1> obj[N];
    int<lower=1> cat[S];
    matrix[T,K] B;
    vector[N] Y;
}
parameters {
    // first layer
    real<lower=0> mu_alpha_bar;
    real<lower=0> sigma_alpha_bar;
    vector[K] mu_w_bar;
    real<lower=0> sigma_w_bar;
    // second layer
    real<lower=0> alpha_bar[E];
    real<lower=0> sigma_alpha;
    vector[K] w_bar[E];
    real<lower=0> sigma_w;
    // third layer
    real<lower=0> alpha[S];
    vector[K] w[S];
    real<lower=0> sigma_y;
}
transformed parameters {
    vector[N] mu;
    for (n in 1:N) {
        mu[n] = alpha[obj[n]] + B[age[n]] * w[obj[n]];
    }
}
model {
    mu_alpha_bar ~ normal(0.22, 40);
    mu_w_bar ~ normal(0.22, 1);
    sigma_alpha_bar ~ exponential(1);
    sigma_w_bar ~ exponential(1);
    for (e in 1:E) {
        alpha_bar[e] ~ normal(mu_alpha_bar, sigma_alpha_bar);
        w_bar[e] ~ normal(mu_w_bar, sigma_w_bar);
    }
    sigma_alpha ~ gamma(10,10);
    sigma_w ~ gamma(10,10);
    for (s in 1:S) {
        alpha[s] ~ normal(alpha_bar[cat[s]], sigma_alpha);
        w[s] ~ normal(w_bar[cat[s]], sigma_w);
    }
    sigma_y ~ exponential(1);

    Y ~ normal(mu, sigma_y);


}
generated quantities{
    vector[N] log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_lpdf(Y[n]|mu[n], sigma_y);
    }
}
