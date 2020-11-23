if(Sys.info()['login'] == 'dashadower') {  # prevent unnecessary modifications
  setwd('/home/dashadower/git_repos/aria/regression/failure_bma')
} else setwd("C:/Users/serim/Documents/GitHub/reliability_prediction")
#########################
source(file.path(getwd(), "impute/mice_imputation.R"))
scriptDir <- getwd()
library(rstan)
library(ggplot2)
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
  rate_matrix[n_state, n_state] = 0
  print("###########")
  print(engine_type)
  print(rate_matrix)
}

####### inspection time ############## 
#pick type 3
engine_type <- 3
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
rate_matrix[n_state, n_state] = 0
# if inspecition is preformed twice more
D_rate <- rate_matrix
n_state = 5
max_allowed_state = 3
repair_state = 2
pm_state = 4
cm_state = 5
initial_state = 1

model_I_t <- stan_model(file.path(getwd(), "optimize/step1/inspection_time.stan"), verbose = FALSE) #approx_deterioration_matrix

onehot_array <- aperm(array(unlist(onehot),c(n_state, length(onehot))))
opt_data_I_t <- list(N=length(onehot), n_state=n_state, state_obs=onehot_array, time_obs=imputed_data$age_ind[imputed_data$engine_ind == engine_type], max_allowed_state=max_allowed_state, repair_state=repair_state, initial_state=initial_state,
                    pm_state = pm_state, cm_state = cm_state, D_rate = D_rate)
print("start sampling")
res_I_t <- optimizing(model_I_t, opt_data_I_t)#, init=list(rate=rep(1.0, n_state-1)))
print("end sampling")


res_I_t$value # target value

xvals <- as.vector(1:30)
e_state <- vector()
for(i in 1:30){
  print("-----------")
  print(i)
  print(unlist(lapply(1:5, function(x) res_I_t$par[paste0("state_t[",i,",",x,"]")])))
  print(sum(unlist(lapply(1:5, function(x) res_I_t$par[paste0("state_t[",i,",",x,"]")]))))
  e_state[i] <- as.vector(1:5) %*% as.vector(unlist(lapply(1:5, function(x) res_I_t$par[paste0("state_t[",i,",",x,"]")])))
  #e_state[i] <- which.max(as.vector(unlist(lapply(1:5, function(x) res_I_t$par[paste0("state_t[",i,",",x,"]")]))))
}


# 결정변수를 시점이 아닌 주기로. if (i == I_t[1]|| i == I_t[2]|| i == I_t[3]) 을 era_1,2,3으로 나눠 각 구간에서 빈도를 결정변수로 바꾸는 방법시도

I_t <- sort(unlist(lapply(1:10, function(x) res_I_t$par[paste0("I_t[",x,"]")])))
ggplot() + aes(x=xvals, y=e_state) + geom_line() + geom_point(color="red") + ggtitle(paste("Expected state at time t target:", res_I_t$value)) +
  geom_vline(xintercept = I_t, linetype="dashed") + ylim(1, 5)
