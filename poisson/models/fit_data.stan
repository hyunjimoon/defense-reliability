functions{
    real failure_form(real phi, real age){
        return (pow(exp(phi), age) - 1)/(exp(phi) - 1);
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
    real<lower=1> phi;
    real alpha;
    real beta;
    real gamma;
    real delta;
    real eta;
}

transformed parameters {
    real early;
    real wear;
    real lambda;
    early = complexity * alpha + beta * log(relative_displacement);
    wear = complexity * engine_count * gamma + delta * log(relative_displacement);
    lambda = early * failure_form(phi, -age + 1) + wear * failure_form(phi, age) + eta;
}

model {
    phi ~ normal(10, 5);
    alpha ~ normal(1.5, 1);
    beta ~ normal(0, 0.5);
    gamma ~ normal(0, 1);
    delta ~ normal(0, 0.5);
    eta ~ normal(0, 1);
    
    y ~ poisson(exp(lambda));
}

generated quantities{
    int post_y[N];
    for(i in 1:N){
        post_y[i] = poisson_rng(exp(lambda));
    }
}