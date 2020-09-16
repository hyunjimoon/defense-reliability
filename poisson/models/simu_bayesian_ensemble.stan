data {
    int N;
    real complexity;
    real<lower=0,upper=1> age;
    real relative_displacement;
    int engine_count;
    
}

generated quantities{
    real<lower=1> phi = fabs(normal_rng(3, 1.5));
    real<lower=0> alpha = fabs(normal_rng(3, 5));
    real<lower=0> beta = fabs(normal_rng(3, 5));
    real<lower=0> gamma = fabs(normal_rng(3, 5));
    real<lower=0> delta = fabs(normal_rng(3, 5));
    real eta = normal_rng(5, 1);

    real early = complexity * alpha + beta * log(relative_displacement);
    real wear = complexity * engine_count * gamma + delta * log(relative_displacement);
    real<lower=0> lambda = early * (pow(phi, age) - 1)/(phi - 1) + wear * (pow(phi, -age + 1) - 1)/(phi - 1) + eta;
    int y[N];
    for(n in 1:N){
        y[n] = poisson_rng(lambda);
    }
}