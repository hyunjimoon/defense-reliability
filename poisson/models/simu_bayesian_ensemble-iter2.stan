functions{
    real failure_form(real shape, real scale, real age){
    
        return fmin((shape/scale) * pow(age/scale, shape-1) * exp(-pow(age/scale, shape)), 2);
    }
    real truncated_normal_rng(real mu, real sigma, real lb){
        real p = normal_cdf(lb, mu, sigma);  // cdf for bounds
        real u = uniform_rng(p, 1);
        return (sigma * inv_Phi(u)) + mu;  // inverse cdf for value
    }
}
data {
    int N;
    real complexity;
    real<lower=0,upper=1> age;
    real relative_displacement;
    int engine_count;
    
}

generated quantities{
    real<lower=0> alpha = truncated_normal_rng(0.05, 0.3, 0.0);
    real<lower=0> beta = truncated_normal_rng(1.2, 0.05, 0.0);
    real<lower=0> gamma = truncated_normal_rng(2.5, 0.5, 0.0);
    real<lower=0> delta = truncated_normal_rng(1.2, 0.1, 0.0);

    real eta = normal_rng(0.5, 0.7);
    real theta = normal_rng(0.2, 0.3);

    real wear = complexity * engine_count * eta;// - theta * log(relative_displacement);
    real lambda = failure_form(alpha, beta, age) + wear * failure_form(gamma, delta, age);
    int y[N];
    for(n in 1:N){
        y[n] = poisson_rng(exp(lambda));
        if(n % 100 == 0){
            print(complexity * engine_count * eta, " ",theta * log(relative_displacement), " ", relative_displacement);
        }
    }
}