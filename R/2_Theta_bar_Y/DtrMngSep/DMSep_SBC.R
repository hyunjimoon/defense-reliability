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
p21 <- 0.1 #true$x[true$X== "p21"]
p31 <- 0.2 #true$x[true$X== "p31"]
rate <- matrix(nrow = P, ncol = S)
for (i in 1:P) for (j in 1:S) rate[i,j] <- true$x[true$X== sprintf("rate[%s,%s]",  i, j)]
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
    # Hyperparameter -> receive as input
    p21 = true_pars[[1]]
    p31 = true_pars[[2]]
    rate = true_pars[[3]]
    # Maintenance
    Mnt = array(c(1, 0, 0, p21, 1- p21, 0, p31, 1-p31,0), dim=c(S,S))
    # Deterioration
    DM_pow <- array(rep(1, S * S * T), dim=c(S, S, T))
    tmp_p <- array(rep(0, S * S), dim=c(S,S))
    Dtr <- array(rep(0, S * S * P), dim=c(S, S, P))
    for(p in 1:P){
      tmp_p[1,1] = exp(-rate[p,1]- rate[p,2]);
      tmp_p[2,1] = rate[p,1] * exp(-rate[p,3]) * (1-exp(-(rate[p,1]+ rate[p,2] - rate[p,3]))) / (rate[p,1]+ rate[p,2] - rate[p,3]);
      tmp_p[3,1] = exp(-rate[p,3]);
      Dtr[,,p] = array(c(tmp_p[1,1], tmp_p[2,1], 1- tmp_p[1,1], 0, tmp_p[3,1], 1 - tmp_p[3,1], 0,0,1))
    }
    DM_pow[,,1] = Dtr[,,1];
    for (t in 2:T){
      if (t <= 8){ DM_pow[,,t] = Dtr[,,1] * Mnt * DM_pow[,,t-1];} # colsum is 1 for right stoch.matrix
      else if (t <=20){ DM_pow[,,t] = Dtr[,,2] * Mnt * DM_pow[,,t-1];}
      else if (t <=26){DM_pow[,,t] = Dtr[,,3] * Mnt * DM_pow[,,t-1];}
      else{DM_pow[,,t] = Dtr[,,4] * Mnt * DM_pow[,,t-1];}
    }
  list(
    generated = lapply(seq(1:length(obs2time)), function(i){sample(c(1,2,3), 1, prob = c(DM_pow[,,obs2time[i]] %*% init_state_vec))}),
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
#write.csv(unlist(as_draws_df(workflow$simulated_y))$y, file = "sim_y.csv")
summary(stan_data$states)
sampling_res<-DMSep$sample(data = stan_data,iter_warmup =1000,  iter_sampling = 1000)
res_df<-as.data.frame(sampling_res)
sample_mean<-apply(res_df,2,mean)

# using SBC package TODO
workflow$fit_model(sample_iterations = SBC_M, warmup_iterations = 500, stan_data)
prior <- workflow$prior_samples
post <- workflow$posterior_samples
workflow$calculate_rank()
plot_ecdf(workflow, var = "theta")
plot_ecdf_diff(workflow, var="theta")
