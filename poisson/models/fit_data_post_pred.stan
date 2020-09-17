functions{
    real failure_form(real phi, real age){
        return (pow(exp(phi), age) - 1)/(exp(phi) - 1);
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
    int ship_number_pred[N];
    
}

parameters {
    real<lower=1> phi;
    real alpha[engine_types];
    real beta[ship_number_max];
    real gamma[engine_types];
    real delta[ship_number_max];
    real eta[ship_number_max];
}

transformed parameters {
    real early[N];
    real wear[N];
    //real lambda[N];
    real lambda[N];
    for(i in 1:N){
        early[i] = complexity[i] * alpha[engine_type[i]] + beta[ship_number[i]] * log(relative_displacement[i]);
        wear[i] = complexity[i] * engine_count[i] * gamma[engine_type[i]] + delta[ship_number[i]] * log(relative_displacement[i]);
        lambda[i] = early[i] * failure_form(phi, -age[i]+1) + wear[i] * failure_form(phi, age[i]) + eta[ship_number[i]];
        lambda[i] = exp(lambda[i]);
    }
}

model {
    phi ~ normal(10, 5);
    for(i in 1:engine_types){
        alpha[i] ~ normal(1.5, 1);
        gamma[i] ~ normal(0, 1);
    }
    for(i in 1:ship_number_max){
        beta[i] ~ normal(0, 0.5);
        delta[i] ~ normal(0, 0.5);
        eta[i] ~ normal(0, 1);
    }
    
    
    for(i in 1:N){
        y[i] ~ poisson(lambda[i]);
    }
    
}

generated quantities{
    int y_post_pred[N_pred];
    for(i in 1:N_pred){
        real early_pred = complexity_pred[i] * alpha[engine_type_pred[i]] + beta[ship_number_pred[i]] * log(relative_displacement_pred[i]);
        real wear_pred = complexity_pred[i] * engine_count_pred[i] * gamma[engine_type_pred[i]] + delta[ship_number_pred[i]] * log(relative_displacement_pred[i]);
        real lambda_pred = early_pred * failure_form(phi, -age_pred[i]+1) + wear_pred * failure_form(phi, age_pred[i]) + eta[ship_number_pred[i]];
        y_post_pred[i] = poisson_rng(exp(lambda_pred));
    }
}