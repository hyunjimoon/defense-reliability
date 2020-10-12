functions{
    real failure_form(real shape, real scale, real age){
        if(shape == 0 && scale == 0){
            return 0;
        }
        return fmin((shape/scale) * pow(age/scale, shape-1) * exp(-pow(age/scale, shape)), 2);
    } 
}

data {
    int N;
    int engine_types; // number of unique engines in engine_type
    int y[N];
    real complexity[N];
    real<lower=0,upper=1> age[N];
    int engine_type[N];
    real relative_displacement[N];
    int engine_count[N];
    int ship_number[N];
    int ship_number_max;

    int N_pred;
    real<lower=0, upper=1> age_pred[N_pred];
    int engine_type_pred[N_pred];
    real complexity_pred[N_pred];
    real relative_displacement_pred[N_pred];
    int engine_count_pred[N_pred];
    int ship_number_pred[N_pred];
}

parameters {
    real<lower=0> alpha[engine_types];
    real<lower=0> beta[engine_types];
    real<lower=0> gamma[engine_types];
    real<lower=0> delta[engine_types];

    real epsilon[engine_types];
    real zeta[ship_number_max];
    real eta[engine_types];
    real theta[ship_number_max];
}

transformed parameters {
    real early[N];
    real wear[N];
    //real lambda[N];
    real lambda[N];
    for(i in 1:N){
        early[i] = complexity[i] * epsilon[engine_type[i]] + zeta[ship_number[i]] * log(relative_displacement[i]);
        wear[i] = complexity[i] * engine_count[i] * eta[engine_type[i]] + theta[ship_number[i]] * log(relative_displacement[i]);

        lambda[i] = early[i] * failure_form(alpha[engine_type[i]], beta[engine_type[i]], age[i]) + wear[i] * failure_form(gamma[engine_type[i]], delta[engine_type[i]], age[i]);

        lambda[i] = exp(lambda[i]);
        /*if(i % 100 == 0){
            print(lambda[i]);
        }*/
    }
}

model {
    for(i in 1:engine_types){
        alpha[i] ~ normal(0.5, 0.3);
        beta[i] ~ normal(1.2, 0.05);
        gamma[i] ~ normal(4, 0.3);
        delta[i] ~ normal(1.2, 0.05);
        
        epsilon[i] ~ normal(1.5, 1);
        eta[i] ~ normal(1.5, 1);
    }
    for(i in 1:ship_number_max){
        zeta[i] ~ normal(0.5, 0.5);
        theta[i] ~ normal(0.5, 0.5);
    }
    
    
    for(i in 1:N){
        y[i] ~ poisson(lambda[i]);
    }
    
}

generated quantities{
    /*int y_post_pred[N_pred];
    for(i in 1:N_pred){
        real early_pred = complexity_pred[i] * epsilon[engine_type_pred[i]] + zeta[ship_number_pred[i]] * log(relative_displacement_pred[i]);

        real wear_pred = complexity_pred[i] * engine_count_pred[i] * eta[engine_type[i]] + theta[ship_number_pred[i]] * log(relative_displacement_pred[i]);

        real lambda_pred = early_pred * failure_form(alpha[engine_type[i]], beta[engine_type[i]], age_pred[i]) + wear_pred * failure_form(gamma[engine_type[i]], delta[engine_type[i]], age_pred[i]);

        y_post_pred[i] = poisson_rng(exp(lambda_pred));
    }*/
}