functions{
    real failure_form(real shape, real scale, real age){
    
        return fmin((shape/scale) * pow(age/scale, shape-1) * exp(-pow(age/scale, shape)), 2);
        //print(shape, " ", scale, " ", age);
        //return shape + scale + age;
    } 
}

data {
    int N; // number of observed data values
    int engine_types; // number of unique engines in engine_type
    int y[N]; // observed data
    real complexity[N]; // feature1
    real<lower=0,upper=1> age[N]; // age at failure
    int engine_type[N]; // unique engine type
    real relative_displacement[N]; // feature2
    int engine_count[N]; // feature3
    int ship_number[N]; // ship type
    int ship_number_max; // number of unique ship types

    int N_pred;
    real<lower=0, upper=1> age_pred[N_pred];
    int engine_type_pred[N_pred];
    real complexity_pred[N_pred];
    real relative_displacement_pred[N_pred];
    int engine_count_pred[N_pred];
    int ship_number_pred[N_pred];
}

parameters {
    real<lower=0> alpha[ship_number_max];
    real<lower=0> beta[ship_number_max];
    real<lower=0> gamma[engine_types];
    real<lower=0> delta[engine_types];

    real eta[engine_types];
    real theta[ship_number_max];
    
    //real early_p1[engine_types];
    //real early_p2[engine_types];
}

transformed parameters {
    real wear[N];
    //real early[N];
    real lambda[N];
    
    for(i in 1:N){
        //early[i] = complexity[i] * early_p1[engine_type[i]] + early_p2[engine_type[i]] * log(relative_displacement[i]);
        wear[i] = complexity[i] * engine_count[i] * eta[engine_type[i]] + theta[ship_number[i]] * log(relative_displacement[i]);

        lambda[i] = failure_form(alpha[ship_number[i]], beta[ship_number[i]], age[i]) * wear[i] + failure_form(gamma[engine_type[i]], delta[engine_type[i]], age[i]) * wear[i];

        lambda[i] = exp(lambda[i]) + 0.00001;
        if(lambda[i] <= 0){
            print("ERROR ERROR", lambda);
        }
    }
}

model {
    for(i in 1:engine_types){
        gamma[i] ~ normal(4, 0.3); //weibull shape
        delta[i] ~ normal(1.2, 0.05); // weibull scale
        
        eta[i] ~ normal(0, 1);
        
        //early_p1[i] ~ normal(1.5, 1);
        //early_p2[i] ~ normal(0, 0.5);
    }
    for(i in 1:ship_number_max){
        alpha[i] ~ normal(0.5, 0.3); //weibull shape
        beta[i] ~ normal(1.2, 0.05); //weibull scale
        
        theta[i] ~ normal(0, 0.5);
    }
    
    
    for(i in 1:N){
        y[i] ~ poisson(lambda[i]);
    }
    
}

generated quantities{
    int y_post_pred[N_pred];
    for(i in 1:N_pred){

        //real early_pred = complexity_pred[i] * early_p1[engine_type_pred[i]] + early_p2[engine_type_pred[i]] * log(relative_displacement_pred[i]);

        real wear_pred = complexity_pred[i] * engine_count_pred[i] * eta[engine_type_pred[i]] + theta[ship_number_pred[i]] * log(relative_displacement_pred[i]);

        real lambda_pred = failure_form(alpha[ship_number_pred[i]], beta[ship_number_pred[i]], age_pred[i]) * wear_pred + wear_pred * failure_form(gamma[engine_type_pred[i]], delta[engine_type_pred[i]], age_pred[i]);

        y_post_pred[i] = poisson_rng(exp(lambda_pred));
    }
}