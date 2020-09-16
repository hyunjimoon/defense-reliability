data {
    int <lower = 1> K; // number of basis functions
    int <lower = 1> N; // 
    int <lower = 1> T; // time length of data(31)
    matrix[T,K] B; // spline values, 2d array
    int <lower = 1> K_i;
    vector [N] y;
    int N_obs;
    vector[N_obs] y_obs;
}

parameters {
    real mu_a_bar;
    real mu_w_bar[K];
    
    real s_a;
    real s_w;
    real s_a_tilde;
    real s_w_tilde;
    
    real <lower = 0> s_a_bar;
    real <lower = 0> s_w_bar;
    
    real s_a_bar_tilde;
    real s_w_bar_tilde;
    
    real <lower=0> s_Y;
    
    
}

transformed parameters {
    real mu;
    vector[K] w;
    for(i in 1:K){
        w[i] = mu_w_bar[i] + s_w_bar * s_w_bar_tilde + s_w * s_w_tilde;
    }
    if(is_nan(B[K_i] * w)){
        print("dot product of B * w is nan: ", B[K_i] * w);
        print(mu_w_bar);
        print(s_w_bar);
        print(s_w_bar_tilde);
        print(s_w);
        print(s_w_tilde);
    }
    
    mu = mu_a_bar + s_a_bar * s_a_bar_tilde + s_a * s_a_tilde + B[K_i] * w;
    if(is_nan(mu)){
    print("mu is also nan: ", B[K_i] * w)
    }
}

model {
    mu_a_bar ~ normal(0, 1);
    for(i in 1:K){
        mu_w_bar[i] ~ normal(0, 1);
    }

    s_a ~ gamma(5, 2);
    s_w ~ gamma(5, 2);
    
    s_a_bar ~ normal(0, 2);
    s_w_bar ~ gamma(2, 1.7);

    s_a_bar_tilde ~ std_normal();
    s_w_bar_tilde ~ std_normal();
    
    s_a_tilde ~ std_normal(); // normal(0,1)
    s_w_tilde ~ std_normal();

    s_Y ~ exponential(1); // exponential(4)
    
    y ~ normal(mu, s_Y);
}

generated quantities{
    vector[N] post_y;
    for(i in 1:N){
        post_y[i] = normal_rng(mu, s_Y);
    }
}
