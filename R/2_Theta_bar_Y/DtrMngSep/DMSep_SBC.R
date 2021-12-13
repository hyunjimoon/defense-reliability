library(SBC)
library(mice)
library(cmdstanr)
library(dplyr)
source(file.path(getwd(), "R/1_Y_bar_y/impute/engine_fail_imputation.R"))
source("R/utils/functions.r")
# Three simulators for prior, data, posterior
# 1. prior simulator
true <- as_tibble(read.csv("./data/DMSep_validation_trueparam.csv"))
true_ratepq<-as_tibble(read.csv("./data/trueparam_lambdapq.csv")[,-c(1,4)])
true_ratepq
# 2. data simulator
# DGP with truth-point benchmarking
# Metadata
N = 3069
T = 31
S = 3
P = 4
p21 <-true_ratepq$mean[4] #0.1 #
p31 <- true_ratepq$mean[5] #0.2 #
rate <- true_ratepq$mean[1:3]
true_pars <-list(p21, p31, rate)


generator_single <- function(N,T,S,P){

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
  states <- sapply(seq(1:length(obs2time)), function(n){sample(c(1,2,3), 1, prob = latent_states[obs2time[n],])})
  stan_data <- list(N= N,T = 31, S = 3, P = 4, states=states, obs2time=obs2time, initial_state=initial_state)
  list(
    generated = stan_data,
    parameters = list(
      p21 = p21,
      p31 = p31,
      rate = rate
    )
  )
}

# Create SBC_Datasets from generator

set.seed(54882235)
n_datasets <- 100  # Number of SBC iterations to run

generator <- SBC_generator_function(generator_single, N=N, T=T, S=S, P=P)
dataset <- generate_datasets(generator, n_datasets)

# Defining backend

use_cmdstanr<-TRUE

model_DMsep <-cmdstanr::cmdstan_model(file.path(getwd(), "R/2_Theta_bar_Y/DtrMngSep/models/DMSep/DMSep.stan"))

if(use_cmdstanr) {
  backend <- SBC_backend_cmdstan_sample(
    model_DMsep, iter_warmup = 1000, iter_sampling = 1000, chains = 4)
} else {
  poisson_backend <- SBC_backend_rstan_sample(
    rstan_model, iter = 2000, warmup = 1000, chains = 2)  
}

# Computing Ranks

cache_dir <- "R/2_Theta_bar_Y/DtrMngSep/basic_usage_SBC_cache"
if(!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}
#results <- compute_results(dataset, backend, 
#                           cache_mode = "results", 
#                           thin_ranks = 30,
#                           cache_location = file.path(cache_dir, "results"))

results <- compute_results(dataset, backend, 
                           thin_ranks = 100)

# Viewing Results

results$stats

# Plots
png(file="R/2_Theta_bar_Y/DtrMngSep/figure/SBC_thin_100.png")
plot_rank_hist(results)
dev.off()




