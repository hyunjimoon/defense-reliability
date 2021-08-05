data {
    int <lower = 1> K;
    int <lower = 1> N;
    int <lower = 1> T; // time length of data
    int <lower = 1> S; // number of ships (99)
    int <lower = 1> E; // number of unique engines (5)
    int <lower = 1> Age[N];
    int <lower = 1> Ship[N];
    int <lower = 1> S2F[S]; //ship_ID to feature map
    matrix[T,K] B;
    vector [N] Y;
}


generated quantities{
    vector[K] w_bar[E];  // defined in transformed parameters
    vector[S] a;  // defined in transformed parameters
    vector[K] w[S];  // defined in transformed parameters
    vector [N] mu;  // defined in transformed parameters
    real s_a_tilde[E];
    real s_w_tilde[E];
    real s_a_bar_tilde[E];
    real s_w_bar_tilde[E];
    real a_bar[E];
    vector[K] mu_w_bar;
    real Y_sim[N];
    real<lower=0> s_a = gamma_rng(10, 10);
    real<lower=0> s_w = gamma_rng(10, 10);
    
    
    real<lower=0> s_Y = exponential_rng(1);
    real mu_a_bar = normal_rng(0, 1);
       real s_a_bar = exponential_rng(1);
    real s_w_bar = exponential_rng(1);
    for(i in 1:K){
        mu_w_bar[i] = normal_rng(0, 1);
    }

    
    for(i in 1:E){
        s_w_bar_tilde[i] = normal_rng(0, 1);
        s_a_bar_tilde[i] = normal_rng(0, 1);
        s_a_tilde[i] = normal_rng(0, 1);
        s_w_tilde[i] = normal_rng(0, 1);
    }


    for (e in 1:E){
        a_bar[e] = mu_a_bar + s_a_bar * s_a_bar_tilde[e];
        w_bar[e] = mu_w_bar + s_w_bar * s_w_bar_tilde[e];
    }
    for (s in 1:S){
        a[s] = a_bar[S2F[s]] + s_a * s_a_tilde[S2F[s]];
        w[s] = w_bar[S2F[s]] + s_w * s_w_tilde[S2F[s]];
    }
  
    for (n in 1: N){
        mu[n] = a[Ship[n]] + B[Age[n]] * w[Ship[n]];
    }
    Y_sim = normal_rng(mu, s_Y);
}