data {
    int <lower = 1> K; // number of basis functions
    int <lower = 1> N; // number of total values(653)
    int <lower = 1> T; // time length of data(31)
    matrix[T,K] B; // spline values, 2d array
    int <lower = 1> K_i;
    int <lower = 1> N_obs; // number of total values(653)
    vector [N_obs] y_data; // actual values
}

generated quantities{
    real mu_a_bar = normal_rng(0, 1);
    real mu_w_bar[K];
    
    real s_a = gamma_rng(5, 2);
    real s_w = gamma_rng(5, 2);
    real s_a_tilde = normal_rng(0, 1);
    real s_w_tilde = normal_rng(0, 1);
    
    real <lower = 0> s_a_bar = abs(normal_rng(0, 2));
    real <lower = 0> s_w_bar = gamma_rng(2, 1.7);
    
    real s_a_bar_tilde = normal_rng(0, 1);
    real s_w_bar_tilde = normal_rng(0, 1);
    
    real <lower=0> s_Y = exponential_rng(1);
    vector[K] w;
    vector[N] sim_y;
    real mu;
    for(i in 1:K) {
        mu_w_bar[i] = normal_rng(0, 1);
    }
    
    for(i in 1:K){
        w[i] = mu_w_bar[i] + s_w_bar * s_w_bar_tilde + s_w * s_w_tilde;
    }
    
    mu = mu_a_bar + s_a_bar * s_a_bar_tilde + s_a * s_a_tilde + B[K_i] * w;
    
    
    for(i in 1:N){
        sim_y[i] = normal_rng(mu, s_Y);
    }
    
}