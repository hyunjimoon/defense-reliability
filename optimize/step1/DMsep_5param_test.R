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
  state <- cut(data, breaks=c(0, 80, 160, max(data)), labels=1:n, include.lowest = TRUE)
  state<-as.numeric(state)
  
  matrix(state,nrow=31)
}

state_matrix <- generate_state_matrix(imputed_data$y_data, n_state)
states <- as.vector(t(state_matrix)) #[imputed_data$engine_ind == engine_type] shiptype -> year
onehot <- list()
# one-hot encode per-data state to vector
for(i in 1:length(states)){
  t_tmp <- rep(0, len=n_state)
  t_tmp[states[i]] <- 1
  onehot[[i]] <- t_tmp
}

onehot_array <- aperm(array(unlist(onehot),c(n_state, length(onehot)))) # array() fills column-wise first

train_data <- onehot_array[1:(31*80), ]
test_data <- onehot_array[(31*80+1):dim(onehot_array)[1] , ]
opt_data <- list(N= dim(train_data)[1], n_state=n_state, state_obs=train_data, time_obs=imputed_data$age_ind[1:(31*80)], initial_state=initial_state)
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

pred_res <- data.frame(matrix(ncol=2, nrow=31*99 - 31*80))
colnames(pred_res) <- c("pred", "obs")
for(i in (31*80+1):length(states)){
  D_pow <- t(matrix(as.vector(unlist(lapply(1:n_state, function(row){ lapply(1:n_state, function(col){res[paste0("D_pow[",imputed_data$age_ind[i],",",row,",",col,"]")]})}))), nrow=3))
  DM_pow <- t(matrix(as.vector(unlist(lapply(1:n_state, function(row){ lapply(1:n_state, function(col){res[paste0("DM_pow[",imputed_data$age_ind[i],",",row,",",col,"]")]})}))), nrow=3))
  if(imputed_data$age_ind[i] == 1){
    pred_res[i-31*80, "pred"] <- which.max(D_pow %*% as.integer(intToBits(2 ** (initial_state-1)))[1:3])
    #print(D_pow %*% as.integer(intToBits(2 ** (initial_state-1)))[1:3])
  }
  else{
    pred_res[i-31*80, "pred"] <- which.max(DM_pow %*% as.integer(intToBits(2 ** (initial_state-1)))[1:3])
    #print(DM_pow %*% as.integer(intToBits(2 ** (initial_state-1)))[1:3])
  }
  pred_res[i-31*80, "obs"] <- which.max(test_data[i-31*80, ])#apply(test_data[i-31*80], 1, which.max)
}

library(ggplot2)
library(reshape2)

ggplot(pred_res, aes(x=as.numeric(row.names(pred_res)))) + geom_point(aes(y=pred, colour="red")) +
  geom_point(aes(y=obs, colour="green")) + ylim(1, 3)

sum(pred_res[, "pred"] == pred_res[, "obs"]) / dim(test_data)[1]
