functions{
    real failure_form(real shape, real age){
        return (pow(exp(shape), age) - 1)/(exp(shape) - 1);
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
    real<lower=1> phi = truncated_normal_rng(10, 5, 1.0);
    real<lower=1> rho = truncated_normal_rng(5, 3, 1.0);
    real alpha = normal_rng(1.5, 1);
    real beta = normal_rng(0, 0.5);
    real gamma = normal_rng(0, 1);
    real delta = normal_rng(0, 0.5);
    real eta = normal_rng(0, 1);

    real early = complexity * alpha + beta * log(relative_displacement);
    real wear = complexity * engine_count * gamma + delta * log(relative_displacement);
    real lambda = early * failure_form(phi, -age+1) + wear * failure_form(rho, age) + eta;
    int y[N];
    for(n in 1:N){
        y[n] = poisson_rng(exp(lambda));
    }
}