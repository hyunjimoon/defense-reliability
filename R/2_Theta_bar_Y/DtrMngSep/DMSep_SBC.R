devtools::install_github("hyunjimoon/SBC",ref="api-variant")
library(SBC)
library(mice)
library(cmdstanr)
library(dplyr)
library("ConnMatTools")
source(file.path(getwd(), "R/1_Y_bar_y/impute/engine_fail_imputation.R"))
source("R/utils/functions.r")
source("R/2_Theta_bar_Y/DtrMngSep/DMsep_5param_test.R")
modelName = "DMSep"
Dir <- set_get_Dir(modelName,  "~/Dropbox/21S_paper/defense-reliability/R/2_Theta_bar_Y/DtrMngSep")
DMSep = cmdstanr::cmdstan_model(Dir$file)

mice_imp <- generateMice()
imputed_data<- complete(mice_imp, 1)
imputed_data_gp <- read.csv("data/y_pred_5var.csv")[,-1]
imputed_data$y_data <- unlist(c(imputed_data_gp))
n_state = 3
initial_state = 1
generate_state_matrix <- function(data, n){
  #state <- cut(data, breaks=c(0, 80, 160, max(data)), labels=1:n, include.lowest = TRUE)
  state <- cut(data, breaks=quantile(data,c(0,1/3,2/3,1)), labels=1:n, include.lowest = TRUE)
  state<-as.numeric(state)
  matrix(state,nrow=31)
}
state_matrix <- generate_state_matrix(imputed_data$y_data, n_state)
states <- as.vector(t(state_matrix))

iter=2000
MSE_df<-data.frame(index=rep(0,iter),p=rep(0,iter),q=rep(0,iter),train_MSE=rep(0,iter),test_MSE=rep(0,iter))
rate_array<-data.frame(period=c(),index=c(),rate=c())
D_array <- array(0, dim=c(iter,4,3,3))
ship_ind_df<-matrix(0,nrow=iter,ncol=5)
# Three simulators for prior, data, posterior
# 1. prior simulator
# 1.1. true_point benchmarking
true_pointest <- as_tibble(read.csv("./data/DMSep_validation_trueparam.csv"))
# 1.2. true_distribution
summ <- summarise_draws(gpImputeDMSepFit(1:99), "mean", "sd")
#true_sample <- sampling_res summ <- summarise_draws(true_sample,  "mean", "sd")
true_p21_hp <- summ[4,c(2,3)]
true_p31_hp <- summ[5,c(2,3)]
true_rate_hp <- summ[c(1,2,3),c(2,3)]
##############
# truth-point benchmarking
##############
true_hp <- list(true_p21_hp, true_p31_hp, true_rate_hp)
p21 <-true$x[true$X== "p21"]  #0.1 #
p31 <- true$x[true$X== "p31"] #0.2 #
rate <- rep(NA, S);
for (j in 1:S) rate[j] <- true$x[true$X== sprintf("rate[%s]", j)]
true_pars <-list(p21, p31, rate)
point_generator <- function(){
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
    tmp_p[1,1] = exp(-rate[1]- rate[2])
    tmp_p[2,1] = rate[1] * exp(-rate[3]) * (1-exp(-(rate[1]+ rate[2] - rate[3]))) / (rate[1]+ rate[2] - rate[3])
    tmp_p[3,1] = exp(-rate[3])
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
##############
# truth-dist. benchmarking
##############
true_hp <- lapply(seq(1:length(true_hp)), function(x){as.numeric(unlist(true_hp[[x]]))})
# 2. data simulator
# DGP with truth-point benchmarking
# Metadata
N = 3069
T = 31
S = 3
P = 4
dist2samp <- function(true_hp){
  estBetaParams <- function(mu, var) {
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
  }
  p21 = rnorm(1, true_hp[[1]][1], true_hp[[1]][2])  #rnorm(1, true_hp[[1]][1], true_hp[[1]][2]) #rbeta(1,p21_hp$alpha, p21_hp$beta)
  p31 = rnorm(1, true_hp[[2]][1], true_hp[[2]][2])  #rbeta(1, p31_hp$alpha, p31_hp$beta)
  for (i in 1:3){
    #ratet_hp <- gammaParamsConvert(mean=true_hp[[3]][1],sd=true_hp[[3]][2])
    rate[i] = rnorm(1,true_hp[[3]][1], true_hp[[3]][2]/3)  #rgamma(1, ratet_hp$shape, ratet_hp$scale)
  }
  return(list(p21, p31, rate))
}
#lapply(1:5, function(x) {dist2samp(true_hp)
true_samples <- list()
for (i in 1:SBC_N){
  true_samples <- c(tmp, list(dist2samp(true_hp)))
}
custom_SBC_generator()
generator_dist <- function(true_hp){
    print(true_hp)
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
    for (i in SBC_N){
      p21 = true_hp[[1]] #rbeta(1,p21_hp$alpha, p21_hp$beta)
      p31 = rnorm(1, true_hp[[2]][1], true_hp[[2]][2] )  #rbeta(1, p31_hp$alpha, p31_hp$beta)
      for (i in 1:3){
        #ratet_hp <- gammaParamsConvert(mean=true_hp[[3]][1],sd=true_hp[[3]][2])
        rate[i] = true_hp[[3]][i]  #rnorm(1,true_hp[[3]][1], true_hp[[3]][2]/3)  #rgamma(1, ratet_hp$shape, ratet_hp$scale)
      }
      # Maintenance
      Mnt = array(c(1, 0, 0, p21, 1- p21, 0, p31, 1-p31,0), dim=c(S,S))
      # Deterioration
      tmp_p <- array(rep(0, S * S), dim=c(S,S))
      tmp_p[1,1] = exp(-rate[1]- rate[2])
      tmp_p[2,1] = rate[1] * exp(-rate[3]) * (1-exp(-(rate[1]+ rate[2] - rate[3]))) / (rate[1]+ rate[2] - rate[3])
      tmp_p[3,1] = exp(-rate[3])
      Dtr <- array(c(tmp_p[1,1], tmp_p[2,1], 1- tmp_p[1,1], 0, tmp_p[3,1], 1 - tmp_p[3,1], 0,0,1), dim=c(S, S))
      latent_states[1,] = Dtr %*% init_state_vec
      for (t in 2:T){
        latent_states[t,] =  (Dtr %*% Mnt) %*% latent_states[t-1,]
      }
    }
    #p21_hp <- estBetaParams(true_hp[[1]][1],(true_hp[[1]][2])^2) this gives all 1
    #p31_hp <- estBetaParams(true_hp[[2]][1],(true_hp[[2]][2])^2)
    list(
      generated = lapply(seq(1:length(obs2time)), function(n){sample(c(1,2,3), 1, prob = latent_states[obs2time[n],])}),
      parameters = list(
        p21 = p21,
        p31 = p31,
        rate = rate
      )
    )
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
SBC_N = 3
SBC_M = 20 # 40
# truth-point
# workflow <- SBCWorkflow$new(DMSep, generator())
# workflow$simulate(SBC_N, true_pointest)
# truth-dist
workflow <- SBCWorkflow$new(DMSep, generator_dist())
workflow$simulate(SBC_N, dist2samp(true_hp))
stan_data$states <- unlist(posterior::draws_of(posterior::subset_draws(workflow$simulated_y, variable="y")$y))
write.csv(stan_data$states, file = "sim_y.csv")
summary(stan_data$states)
sampling_res<-DMSep$sample(data = stan_data,iter_warmup =1000,  iter_sampling = 1000, parallel_chains = 4)

workflow$fit_model(sample_iterations = SBC_M, warmup_iterations = 500, stan_data)
prior <- workflow$prior_samples
post <- workflow$posterior_samples
workflow$calculate_rank()
plot_ecdf(workflow, var = "theta")
plot_ecdf_diff(workflow, var="theta")


########
cmdgenerator_func <- custom_SBC_generator(generator_dist, true_samples)
datasets_func <- generate_datasets(cmdgenerator_func, SBC_N)
backend <- cmdstan_sample_SBC_backend(DMSep, iter_warmup = 200, iter_sampling = 200)
results <- compute_results(datasets, backend)
plot_ecdf_diff(results)
########

