library(SBC)
library(mice)
library(cmdstanr)
library(dplyr)
source(file.path(getwd(), "R/1_Y_bar_y/impute/engine_fail_imputation.R"))
source("R/utils/functions.r")
# Three simulators for prior, data, posterior
# 1. prior simulator
true <- as_tibble(read.csv("./data/DMSep_validation_trueparam.csv"))
# 2. data simulator
# DGP with truth-point benchmarking
# Metadata
N = 3069
T = 31
S = 3
P = 4
p21 <-true$x[true$X== "p21"]  #0.1 #
p31 <- true$x[true$X== "p31"] #0.2 #
rate <- rep(NA, S);
for (j in 1:S) rate[j] <- true$x[true$X== sprintf("rate[%s]", j)]
true_pars <-list(p21, p31, rate)

generator <- function(){
  function(iter, true_pars){
    # Metadata
    N = 3069
    T = 31
    S = 3
    P = 4
    obs2time = complete(generateMice(), 1)$age_ind
    initial_state = 1
    init_state_vec = rep(0, S);
    init_state_vec[initial_state] = 1
    latent_states = array(NA, dim=c(T,S))
    # Hyperparameter -> receive as input
    p21 = true_pars[[1]]
    p31 = true_pars[[2]]
    rate = true_pars[[3]]
    # Maintenance
    Mnt = array(c(1, 0, 0, p21, 1- p21, 0, p31, 1-p31,0), dim=c(S,S))
    # Deterioration
    tmp_p <- array(rep(0, S * S), dim=c(S,S))
    tmp_p[1,1] = exp(-rate[1]- rate[2]);
    tmp_p[2,1] = rate[1] * exp(-rate[3]) * (1-exp(-(rate[1]+ rate[2] - rate[3]))) / (rate[1]+ rate[2] - rate[3]);
    tmp_p[3,1] = exp(-rate[3]);
    Dtr <- array(c(tmp_p[1,1], tmp_p[2,1], 1- tmp_p[1,1], 0, tmp_p[3,1], 1 - tmp_p[3,1], 0,0,1), dim=c(S, S))
    latent_states[1,] = Dtr %*% init_state_vec;
    for (t in 2:T){
      latent_states[t,] =  (Dtr %*% Mnt) %*% latent_states[t-1,];
    }
  list(
    generated = lapply(seq(1:length(obs2time)), function(n){sample(c(1,2,3), 1, prob = latent_states[obs2time[n],])}),
    parameters = list(
      p21 = p21,
      p31 = p31,
      rate = rate
    )
  )
  }
}

# 3. posterior simulator: prior and data simulator embedded in stan with supplied data
# 3-1. prior and data simulator
modelName = "DMSep"
Dir <- set_get_Dir(modelName,  "R/2_Theta_bar_Y/DtrMngSep")
DMSep = cmdstanr::cmdstan_model(Dir$file)

# 3-2. simulated y data template with init as real y
obs2time <- complete(generateMice(), 1)$age_ind # data template only
states <- rep(3, N) # data template
initial_state <- 1
stan_data <- list(N= N,T = max(obs2time), S = 3, P = 4, states=states, obs2time=obs2time, initial_state=initial_state)
SBC_N = 1
SBC_M = 500 # 40
workflow <- SBCWorkflow$new(DMSep, generator())
workflow$simulate(SBC_N, true_pars)
stan_data$states <- unlist(posterior::draws_of(posterior::subset_draws(workflow$simulated_y, variable="y")$y))
write.csv(stan_data$states, file = "sim_y.csv")
summary(stan_data$states)
sampling_res<-DMSep$sample(data = stan_data,iter_warmup =1000,  iter_sampling = 1000, parallel_chains = 4)
