source(file.path(getwd(), "impute/mice_imputation.R"))
scriptDir <- getwd()
library(rstan)
model <- stan_model(file.path(getwd(), "optimize/step1/approx_deterioration_matrix.stan"), verbose = TRUE)


mice_imp <- generateMice()
imputed_data <- complete(mice_imp, 1)

####################################
state_count = 5
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

state_matrix <- generate_state_matrix(imputed_data$y_data, state_count)
# one-hot encode to vector

options(scipen = 999)
for(engine_type in 1:5){
  states <- as.vector(t(state_matrix))[imputed_data$engine_ind == engine_type]
  onehot <- list()
  for(i in 1:length(states)){
    t_tmp <- as.vector(rep(0, len=state_count))
    t_tmp[states[i]] <- 1
    onehot[[i]] <- t_tmp
  }
  opt_data <- list(N=length(onehot), n_states=state_count, state_obs=onehot, time_obs=imputed_data$age_ind[imputed_data$engine_ind == engine_type], max_allowed_state=max_allowed_state, repair_state=repair_state, initial_state=initial_state)
  res <- optimizing(model, opt_data)
  print(paste0("return code:", res$return_code))
  res <- res$par
  mat <- matrix(rep(0.0, 25), nrow=state_count)
  for(i in 1:state_count){mat[i,i] <- 1}
  for(state in 1:(state_count-1)){
    for(col in 1:(state_count+1-state)){
      mat[state, col + (state-1)] <- as.numeric(res[paste0("state_",state,"[",col,"]")])
    }
  }
  print("###########")
  print(engine_type)
  #mat <- t(matrix(mat, nrow = 5))
  print(mat)
  print(rowSums(mat))
}
options(scipen=0)
