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
  obs2time <- complete(generateMice(), 1)$age_ind # data template only
  states <- rep(3, N) # data template
  initial_state <- 1
  stan_data <- list(N= N,T = max(obs2time), S = 3, P = 4, states=states, obs2time=obs2time, initial_state=initial_state)
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
results <- compute_results(dataset, backend, 
                           cache_mode = "results", 
                           cache_location = file.path(cache_dir, "results"))

# Viewing Results

results$stats

# Plots
png(file="R/2_Theta_bar_Y/DtrMngSep/SBC_results.png",
    width=600, height=350)
plot_rank_hist(results)
dev.off()




