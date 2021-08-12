devtools::install_github("hyunjimoon/SBC",ref="api-variant")
library(SBC)
library(mice)
library(cmdstanr)
library(dplyr)
library(posterior)
library(ConnMatTools)
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

sampling_res <- gpImputeDMSepFit(1:99)
saveRDS(sampling_res, "sample_genstate.RDS")
true_prival <- subset_draws(as_draws_rvars(sampling_res), chain = 1, iteration = seq(1, 981, by = 20))
# Three ways for prior input: truth-point bencmarking, hyperparameter-based truth-dist, sample-based truth-dist,
# The first only simulate with one set of value of parameter values while the second and the third simulate with
# different set of prior values. The second is when the modeler only knows the hyperparmeter of prior values while
# the third is when the modeler have exact n_datasets pairs of prior value in hand.

##############
# truth-dist. with custom priorval
##############
custom_nprior_generator <- function(true_prival){
  n_datasets <- niterations(true_prival)
  generated <- list()
  N = 3069
  T = 31
  S = 3
  P = 4
  obs2time = complete(generateMice(), 1)$age_ind
  initial_state = 1
  init_state_vec = rep(0, S);
  init_state_vec[initial_state] = 1
  latent_states = matrix(rep(NA, T*S),nrow = T, ncol = S)
  gen_states = array(NA, N)
  rate_vec <- rep(NA, 3)
  for(iter in 1:n_datasets){
    true_prival_i <- subset_draws(true_prival, iteration = iter)
    rate_i <- draws_of(true_prival_i$rate)
    p21_i <- draws_of(true_prival_i$p21)
    p31_i <- draws_of(true_prival_i$p31)
    # Maintenance
    Mnt = array(c(1, 0, 0, p21_i, 1- p21_i, 0, p31_i, 1-p31_i,0), dim=c(S,S))
    # Deterioration
    tmp_p <- array(rep(0, S * S), dim=c(S,S))
    tmp_p[1,1] <- exp(- rate_i[,1]-  rate_i[,2])
    tmp_p[2,1] <-  rate_i[,1] * exp(- rate_i[,3]) * (1-exp(-( rate_i[,1] +  rate_i[,2] -  rate_i[,3]))) / ( rate_i[,1] +  rate_i[,2] -  rate_i[,3])
    tmp_p[3,1] <- exp(- rate_i[,3])
    Dtr <- array(c(tmp_p[1,1], tmp_p[2,1], 1- tmp_p[1,1] - tmp_p[2,1], 0, tmp_p[3,1], 1 - tmp_p[3,1], 0,0,1), dim=c(S, S))
    latent_states[1,] <- Dtr %*% init_state_vec
    for (t in 2:T){
      latent_states[t,] <-  (Dtr %*% Mnt) %*% latent_states[t-1,]
    }
    gen_states <- unlist(lapply(seq(1:length(obs2time)), function(n){sample(c(1,2,3), 1, prob = latent_states[obs2time[n],])}))
    generated[[iter]] <- list(N = N, T = T, S = S, P = P, states = gen_states, obs2time = obs2time, initial_state = initial_state)
  }
  SBC_datasets(
    parameters = posterior::as_draws_matrix(posterior::draws_rvars(p21 = true_prival$p21, p31 = true_prival$p31, rate = true_prival$rate)),
    generated <- generated
  )
}
datasets <- custom_nprior_generator(true_prival)
backend <- cmdstan_sample_SBC_backend(DMSep, iter_warmup = 500, iter_sampling = 200)
results_new <- compute_results(datasets, backend)
plot_ecdf_diff(results_new)
