source(file.path(getwd(), "impute/mice_imputation.R"))
scriptDir <- getwd()
library(rstan)
library(ggplot2)
library(pracma)
model <- stan_model(file.path(getwd(), "optimize/step1/approx_deterioration_matrix.stan"), verbose = FALSE)


mice_imp <- generateMice()
imputed_data <- complete(mice_imp, 1)

####################################
state_count = 5
max_allowed_state = 3
repair_state = 2

###########################

"%^%" <- function(mat,power){
  base = mat
  out = diag(nrow(mat))
  while(power > 1){
    if(power %% 2 == 1){
      out = out %*% base
    }
    base = base %*% base
    power  = power %/% 2
  }
  out %*% base
}

generate_maintenance_matrix <- function(n_states, max_allowed_state, repair_state){
  mat <- matrix(rep(0, len=n_states ** 2), nrow=n_states)
  for(i in 1:max_allowed_state){
    mat[i, i] <- 1
  }
  for(i in (max_allowed_state+1):state_count){
    mat[i, repair_state] <- 1
  }
  
  return(mat)
}

generate_state_matrix <- function(data, n){
  state<-cut(data, breaks=c(quantile(data,seq(0,1,length.out = n+1))),labels = 1:n, include.lowest=TRUE)
  state<-as.numeric(as.character(state))
  matrix(state,nrow=31)
}

quantile_vector <- tapply(imputed_data$y_data, findInterval(imputed_data$y_data, quantile(imputed_data$y_data,seq(0,1,length.out = state_count))), mean)

###########################

optimization_horizon = 120
cm_cost = 1.0
pm_cost = 15.0

###########

generic_cost_function <- function(maintenance_interval, initial_state_vector, maintenance_matrix, deterioration_matrix, quantile_vector, horizon, cm_cost, pm_cost){
  cycle_count = floor(horizon / maintenance_interval)
  cm_count = 0
  for(i in 1:horizon){
    initial_state_vector = initial_state_vector %*% (maintenance_matrix %*% deterioration_matrix)
    cm_count = cm_count + as.double(quantile_vector %*% initial_state_vector)
  }
  return(cm_count * cm_cost + cycle_count * pm_cost)
}

minimize_max_failure <- function(maintenance_interval, initial_state_vector, maintenance_matrix, deterioration_matrix, horizon, state_threshold){
  n_states = length(initial_state_vector)
  probs = list()
  for(i in 1:horizon){
    initial_state_vector = initial_state_vector %*% (maintenance_matrix %*% deterioration_matrix)
    probs[i] <- sum(initial_state_vector[state_threshold:n_states])
  }
  return(mean(as.double(probs)))
}

#########################
plot_interval_limit = 5
initial_state = 1

####################
state_matrix <- generate_state_matrix(imputed_data$y_data, state_count)
# one-hot encode to vector

initial_state_vector <- as.vector(rep(0, state_count))
initial_state_vector[initial_state] <- 1

maintenance_matrix <- generate_maintenance_matrix(state_count, max_allowed_state, repair_state)

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
  res <- res$par
  mat <- matrix(rep(0.0, state_count ** 2), nrow=state_count)
  for(i in 1:state_count){mat[i,i] <- 1}
  for(state in 1:(state_count-1)){
    for(col in 1:(state_count+1-state)){
      mat[state, col + (state-1)] <- as.numeric(res[paste0("state_",state,"[",col,"]")])
    }
  }
  print("###########")
  print(engine_type)
  intervals = seq.int(from=1, to=plot_interval_limit)
  result = list()
  for(i in 1:plot_interval_limit){
    #result[i] <- generic_cost_function(i, initial_state_vector, maintenance_matrix, rootm(mat, i)$B, quantile_vector, optimization_horizon, cm_cost, pm_cost)[1]
    result[i] <- minimize_max_failure(i, initial_state_vector, maintenance_matrix, rootm(mat, i)$B, optimization_horizon, max_allowed_state)
    #print(result[i][1])
  }
  result <- as.double(result)
  dd <- data.frame(intervals, result)
  plt <- ggplot(dd, aes(x=intervals, y=result)) + geom_point() + ggtitle(engine_type)
  show(plt)
  #mat <- t(matrix(mat, nrow = 5))
  print(mat)
}
options(scipen=0)
