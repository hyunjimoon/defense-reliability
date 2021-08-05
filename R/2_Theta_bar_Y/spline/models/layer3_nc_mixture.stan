// order: C > E > S -- Y
data {
    int <lower = 1> K; // number of basis functions
    int <lower = 1> N; // number of total values
    int <lower = 1> T; // time length of data
    int <lower = 1> S; // number of x_atom  e.g ship
    int<lower=1> C; // number of categories
    int <lower = 1> E; // max number of types in each cat e.g 5 engine
    int <lower = 1> age[N];
    int <lower = 1> x[N]; // x_id to feature mapping
    int <lower = 1> cat[C, S]; // to category map
    matrix[T,K] B; // spline values, 2d array
    vector[N] Y; // obs
}

parameters {
  // first layer - meta cat
  real<lower = 0> s_a_bar[C];
  real<lower = 0> s_w_bar[C];
  vector[E] eta_a_bar[C];
  vector[E] eta_w_bar[C];
  real mu_a_bar;
  vector[K] mu_w_bar;
  // second layer- per cat
  real<lower=0> s_a;
  real<lower=0> s_w;
  vector[S] eta_a;
  vector[S] eta_w;
  // btw cat
  simplex[C] theta;
  real<lower=0> s_Y;
}

transformed parameters {
  vector[E] a_bar[C];
  vector[K] w_bar[C, E];
  vector[S] a[C];
  vector[K] w[C, S];
  matrix [C, N] mu;
  vector[N] yhat;
  for (c in 1:C){
    for (e in 1:E){
      a_bar[c,e] = mu_a_bar + s_a_bar[c] * eta_a_bar[c,e];// e = 1+ 1 * e
      w_bar[c, e] = mu_w_bar + s_w_bar[c] * eta_w_bar[c, e]; //k = k + 1*1, K broadcast
    }
    for (s in 1:S){
      a[c, s] = a_bar[c, cat[c, s]] + s_a * eta_a[s];
      w[c, s] = w_bar[c, cat[c, s]] + s_w * eta_w[s];
    }
    for (n in 1:N){
      mu[c, n] = a[c, x[n]] + B[age[n]] * w[c, x[n]]; // 1*k k*1
    }
  }
  yhat = mu' * theta; // c*n
}

model {
    // mean effect
    mu_a_bar ~ normal(0, 1);
    mu_w_bar ~ normal(0, 1);
    // inv-var effect
    s_a_bar ~ normal(0, 2);// tuning needed
    s_w_bar ~ gamma(2, 1);
    s_a ~ gamma(5, 2); // tuning needed //s_w,a ~ exp(0.5);
    s_w ~ gamma(5, 2);
    // std error of effect
    for (s in 1:S) {
      eta_a[s] ~  normal(0, 1);
      eta_w[s] ~ normal(0, 1);
    }
    for (c in 1:C) {
     for (e in 1:E){
      eta_a_bar[c, e] ~ normal(0, 1);
      eta_w_bar[c, e] ~ normal(0, 1);
      }
    }
    s_Y ~ exponential(1); // tuning needed . exponential(4)
    Y ~ normal(yhat, s_Y);
}

generated quantities{
    vector[N] log_lik;
    for (i in 1:N) {
        log_lik[i] = normal_lpdf(Y[i]|yhat[i], s_Y);
    }
}
