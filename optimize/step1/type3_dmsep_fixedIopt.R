source(file.path(getwd(), "impute/mice_imputation.R"))
scriptDir <- getwd()
library(rstan)
library(ggplot2)
model_DMsep <- stan_model(file.path(getwd(), "optimize/step1/approx_deterioration_rate_matrix.stan"), verbose = FALSE) #approx_deterioration_matrix

mice_imp <- generateMice()
imputed_data <- complete(mice_imp, 1)


#################################### DM_sep
# original policy (wihtout pm)
n_state = 5
max_allowed_state = 4 # TODO change to cm_init
repair_state = 2
initial_state = 1
engine_type <- 3
generate_state_matrix <- function(data, n){
  state<-cut(data, breaks=c(quantile(data,seq(0,1,length.out = n+1))),labels = 1:n, include.lowest=TRUE)
  state<-as.numeric(as.character(state))
  matrix(state,nrow=31)
}

state_matrix <- generate_state_matrix(imputed_data$y_data, n_state)
states <- as.vector(t(state_matrix))[imputed_data$engine_ind == engine_type]
onehot <- list()
# one-hot encode per-data state to vector
for(i in 1:length(states)){
  t_tmp <- as.vector(rep(0, len=n_state))
  t_tmp[states[i]] <- 1
  onehot[[i]] <- t_tmp
}
# onehot_array<-array(unlist(onehot),c(length(onehot),n_state))  # array() fills column-wise first 
onehot_array <- aperm(array(unlist(onehot),c(n_state, length(onehot)))) #TODO choose onehot_array btw two
opt_data <- list(N=sum(imputed_data$engine_ind == engine_type), n_state=n_state, state_obs=onehot_array, time_obs=imputed_data$age_ind[imputed_data$engine_ind == engine_type], max_allowed_state=max_allowed_state, repair_state=repair_state, initial_state=initial_state)
print("start sampling")
res <- optimizing(model_DMsep, opt_data)#, init=list(rate=rep(1.0, n_state-1)))
print("end sampling")
print(paste0("return code:", res$return_code))
res <- res$par
rate_matrix <- matrix(rep(0.0, n_state ** 2), nrow=n_state)
for(i in 1:(n_state-1)){
  rate_matrix[i, i] = -as.numeric(res[paste0("rate","[",i,"]")]);
  rate_matrix[i, i+1] = as.numeric(res[paste0("rate","[",i,"]")]);
}
rate_matrix[n_state, n_state] = 0

#################################### fixedI_opt
# proposed policy (wiht pm 3,4 -> 1 and cm 5 ->2)
D_rate <- rate_matrix
n_state = 5
cm_repair = 2
pm_repair = 1
pm_init = 3
cm_init = 5
initial_state = 1

model_fixedI <- stan_model(file.path(getwd(), "optimize/step1/inspection_interval_marginalize.stan"), verbose = TRUE)


opt_data_I_t <- list(N=length(onehot), n_state=n_state, state_obs=onehot_array, time_obs=imputed_data$age_ind[imputed_data$engine_ind == engine_type], 
                     cm_init = cm_init, pm_init = pm_init, cm_repair = cm_repair, pm_repair = pm_repair, initial_state=initial_state, D_rate = D_rate)

comb_res <- rstan::sampling(model_fixedI, opt_data_I_t, algorithm="Fixed_param", chains=1, iter = 100)
sum <- summary(comb_res, c("cost"))$summary

