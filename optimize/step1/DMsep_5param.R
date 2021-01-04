source(file.path(getwd(), "impute/mice_imputation.R"))
scriptDir <- getwd()
library(rstan)
library(ggplot2)
model_DMsep <- stan_model(file.path(getwd(), "optimize/step1/DMsep_5param.stan"), verbose = FALSE) #approx_deterioration_matrix

mice_imp <- generateMice()
imputed_data <- complete(mice_imp, 1)

#################################### DM_sep
# original policy (wihtout pm)
n_state = 3
initial_state = 1
#repair_state 
generate_state_matrix <- function(data, n){
  state<-cut(data, breaks=c(quantile(data,seq(0,1,length.out = n+1))),labels = 1:n, include.lowest=TRUE)
  state<-as.numeric(as.character(state))
  matrix(state,nrow=31)
}

state_matrix <- generate_state_matrix(imputed_data$y_data, n_state)
states <- as.vector(t(state_matrix)) #[imputed_data$engine_ind == engine_type]
onehot <- list()
# one-hot encode per-data state to vector
for(i in 1:length(states)){
  t_tmp <- as.vector(rep(0, len=n_state))
  t_tmp[states[i]] <- 1
  onehot[[i]] <- t_tmp
}

onehot_array <- aperm(array(unlist(onehot),c(n_state, length(onehot)))) # array() fills column-wise first 
opt_data <- list(N= length(imputed_data$y_data), n_state=n_state, state_obs=onehot_array, time_obs=imputed_data$age_ind, initial_state=initial_state)
print("start sampling")
res <- optimizing(model_DMsep, opt_data)#, init=list(rate=rep(1.0, n_state-1)))
print("end sampling")
print(paste0("return code:", res$return_code))
res <- res$par
print(res["p"])
rate_matrix <- matrix(rep(0.0, n_state ** 2), nrow=n_state)
for(i in 1:(n_state-1)){
  rate_matrix[i, i] = -as.numeric(res[paste0("rate","[",i,"]")]);
  rate_matrix[i, i+1] = as.numeric(res[paste0("rate","[",i,"]")]);
}
rate_matrix[n_state, n_state] = 0
print(rate_matrix)
TPM_matrix <- matrix(rep(0.0, n_state ** 2), nrow=n_state)
for(i in 1:(n_state)){
  for(j in 1:(n_state)){
  TPM_matrix[i, j] = as.numeric(res[paste0("TPM[1,",i,",", j,"]")]);
  }
}

print(TPM_matrix)
