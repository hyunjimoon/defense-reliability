functions{
    real failure_form(real shape, real scale, real age){
    
        return fmin((shape/scale) * pow(age/scale, shape-1) * exp(-pow(age/scale, shape)), 2);
    } 
}

data {
    int N;
    int y[N];
    real complexity;
    real<lower=0,upper=1> age;
    real relative_displacement;
    int engine_count;
    
}

parameters {
    real<lower=0> alpha;
    real<lower=0> beta;
    real<lower=0> gamma;
    real<lower=0> delta;

    real epsilon;
    real zeta;
    real eta;
    real theta;
}

transformed parameters {
    real early;
    real wear;
    real lambda;
    early = complexity * epsilon + zeta * log(relative_displacement);
    wear = complexity * engine_count * eta + theta * log(relative_displacement);
    lambda = early * failure_form(alpha, beta, age) + wear * failure_form(gamma, delta, age);
}

model {
    alpha ~ lognormal(0, 0.5);
    beta ~ normal(1.5, 0.3);
    gamma ~ normal(4, 1.5);
    delta ~ normal(1.5, 0.3);
    
    epsilon ~ normal(1.5, 1);
    zeta ~ normal(0.5, 0.5);
    eta ~ normal(1.5, 1);
    theta ~ normal(0.5, 0.5);
    
    y ~ poisson(exp(lambda));
}

generated quantities{
    int post_y[N];
    for(i in 1:N){
        post_y[i] = poisson_rng(exp(lambda));
    }
}