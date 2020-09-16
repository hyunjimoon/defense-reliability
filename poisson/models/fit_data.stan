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
    real<lower=0> alpha;
    real<lower=0> beta;
    real<lower=0> gamma;
    real<lower=0> delta;
    real eta;
}

transformed parameters {
    real early;
    real wear;
    real<lower=0> lambda;
    early = complexity * alpha + beta * log(relative_displacement);
    wear = complexity * engine_count * gamma + delta * log(relative_displacement);
    lambda = early * (pow(phi, age) - 1)/(phi - 1) + wear * (pow(phi, -age + 1) - 1)/(phi - 1) + eta;a
}

model {
    phi ~ normal(3, 1.5);
    alpha ~ normal(3, 5);
    beta ~ normal(3, 5);
    gamma ~ normal(3, 5);
    delta ~ normal(3, 5);
    eta ~ normal(5, 1);
    
    y ~ poisson(lambda);
}

generated quantities{
    int post_y[N];
    for(i in 1:N){
        post_y[i] = poisson_rng(lambda);
    }
}