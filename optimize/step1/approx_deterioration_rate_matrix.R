if(Sys.info()['login'] == 'dashadower') {  # prevent unnecessary modifications
  setwd('/home/dashadower/git_repos/aria/regression/failure_bma')
} else setwd("C:/Users/serim/Documents/GitHub/reliability_prediction")
#########################
source(file.path(getwd(), "impute/mice_imputation.R"))
scriptDir <- getwd()
library(rstan)
model <- stan_model(file.path(getwd(), "optimize/step1/approx_deterioration_rate_matrix.stan"), verbose = FALSE) #approx_deterioration_matrix



mice_imp <- generateMice()
imputed_data <- complete(mice_imp, 1)

####################################
n_state = 5
max_allowed_state = 3
repair_state = 2
initial_state = 1

generate_maintenance_matrix <- function(n_states, max_allowed_state, repair_state){
  mat <- matrix(rep(0, len=n_states ** 2), nrow=n_states)
  for(i in 1:max_allowed_state){
    mat[i, i] <- 1
  }
  for(i in max_allowed_state+1:n_states){
    mat[i, repair_state] <- 1
  }
  
  return(mat)
}

generate_state_matrix <- function(data, n){
  state<-cut(data, breaks=c(quantile(data,seq(0,1,length.out = n+1))),labels = 1:n, include.lowest=TRUE)
  state<-as.numeric(as.character(state))
  matrix(state,nrow=31)
}

state_matrix <- generate_state_matrix(imputed_data$y_data, n_state)

options(scipen = 999)

for(engine_type in 1:5){
  states <- as.vector(t(state_matrix))[imputed_data$engine_ind == engine_type]
  onehot <- list()
  # one-hot encode per-data state to vector
  for(i in 1:length(states)){
    t_tmp <- as.vector(rep(0, len=n_state))
    t_tmp[states[i]] <- 1
    onehot[[i]] <- t_tmp
  }
  #onehot_array<-array(unlist(onehot),c(length(onehot),n_state))  # array() fills column-wise first
  onehot_array <- aperm(array(unlist(onehot),c(n_state, length(onehot))))
  opt_data <- list(N=length(onehot), n_state=n_state, state_obs=onehot_array, time_obs=imputed_data$age_ind[imputed_data$engine_ind == engine_type], max_allowed_state=max_allowed_state, repair_state=repair_state, initial_state=initial_state)
  print("start sampling")
  res <- optimizing(model, opt_data)#, init=list(rate=rep(1.0, n_state-1)))
  print("end sampling")
  res["rate[1]"]
  print(paste0("return code:", res$return_code))
  res <- res$par
  rate_matrix <- matrix(rep(0.0, n_state ** 2), nrow=n_state)
  for(i in 1:(n_state-1)){
    rate_matrix[i, i] = -as.numeric(res[paste0("rate","[",i,"]")]);
    rate_matrix[i, i+1] = as.numeric(res[paste0("rate","[",i,"]")]);
  }
  rate_matrix[n_state, n_state] = 1
  print("###########")
  print(engine_type)
  print(rate_matrix)
}

#options(scipen=0)

