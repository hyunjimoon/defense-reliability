#setwd("C:/Users/serim/Documents/GitHub/reliability_prediction")

source(file.path(getwd(), "impute/mice_imputation.R"))
scriptDir <- getwd()
library(rstan)
library(ggplot2)
model_DMsep <- stan_model(file.path(getwd(), "optimize/step1/DMsep_5param.stan"), verbose = FALSE) #approx_deterioration_matrix
imputed_data_gp <- read.csv("data/y_pred_5var.csv")[,-1]
imputed_data$y_data <- unlist(c(imputed_data_gp))

#################################### DM_sep
# original policy (wihtout pm)
n_state = 3
initial_state = 1

generate_state_matrix <- function(data, n){
  state <- cut(data, breaks=quantile(data,c(0,1/3,2/3,1)), labels=1:n, include.lowest = TRUE)
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
set.seed(210106)
test_ship_ind=sort(sample(1:99,5))
test_ind=c(sapply(test_ship_ind,function(i) i*31+(1:31)))
train_data <- onehot_array[-test_ind, ]
test_data <- onehot_array[test_ind, ]
opt_data <- list(N= dim(train_data)[1], n_state=n_state, state_obs=train_data, time_obs=imputed_data$age_ind[-test_ind], initial_state=initial_state)
res <- optimizing(model_DMsep, importance_resampling = TRUE,draws = 5, algorithm = "Newton",  hessian = TRUE, opt_data, verbose = TRUE)

print(res$theta_tilde[13]) # p - unlike $par, $theta_tilde is not `Named num` yet
res$par["p"]
res <- optimizing(model_DMsep, opt_data, verbose = TRUE)

print(paste0("return code:", res$return_code))
############################## 
##############################
# Debug code
for (era in 1:4){
  D <- matrix(as.vector(unlist(lapply(1:n_state, function(row){lapply(1:n_state, function(col){res$par[paste0("D[",era,",",row,",",col,"]")]})}))), nrow=3, byrow=T)
  for(j in 1:3){
    rate <- res$par[paste0("rate[", era,",", j,"]")]
    print(rate)
  }
  print(D)
}

for (era in 1:4){
  DM_pow <- matrix(as.vector(unlist(lapply(1:n_state, function(row){lapply(1:n_state, function(col){res$par[paste0("DM_pow[",era,",",row,",",col,"]")]})}))), nrow=3, byrow=T)
  print(DM_pow)
}
############################## 
##############################

pred_res <- data.frame(matrix(ncol=2, nrow=31*99 - 31*80))
colnames(pred_res) <- c("pred", "obs")
for(i in (31*80+1):length(states)){
  D_pow <- matrix(as.vector(unlist(lapply(1:n_state, function(row){ lapply(1:n_state, function(col){res[paste0("D_pow[",imputed_data$age_ind[i],",",row,",",col,"]")]})}))), nrow=3,byrow=T)
  DM_pow <- matrix(as.vector(unlist(lapply(1:n_state, function(row){ lapply(1:n_state, function(col){res[paste0("DM_pow[",imputed_data$age_ind[i],",",row,",",col,"]")]})}))), nrow=3,byrow=T)
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
