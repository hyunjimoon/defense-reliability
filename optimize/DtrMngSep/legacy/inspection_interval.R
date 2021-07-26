source(file.path(getwd(), "optimize/step1/approx_deterioration_rate_matrix.R"))

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
opt_data <- list(N=sum(imputed_data$engine_ind == engine_type), n_state=n_state, state_obs=onehot_array, time_obs=imputed_data$age_ind[imputed_data$engine_ind == engine_type], max_allowed_state=max_allowed_state, repair_state=repair_state, initial_state=initial_state)
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
max_allowed_state = 4 # = pm_state
repair_state = 2
pm_state = 3
cm_state = max_allowed_state + 1
initial_state = 1

model_I_t <- stan_model(file.path(getwd(), "optimize/step1/inspection_interval.stan"), verbose = FALSE) #approx_deterioration_matrix

onehot_array <- aperm(array(unlist(onehot),c(n_state, length(onehot))))
opt_data_I_t <- list(N=length(onehot), n_state=n_state, state_obs=onehot_array, time_obs=imputed_data$age_ind[imputed_data$engine_ind == engine_type], max_allowed_state=max_allowed_state, repair_state=repair_state, initial_state=initial_state,
                     pm_state = pm_state, cm_state = cm_state, D_rate = D_rate)
print("start sampling")
res_I_t <- optimizing(model_I_t, opt_data_I_t)#, init=list(rate=rep(1.0, n_state-1)))
print("end sampling")


res_I_t$value # target value

xvals <- as.vector(1:31)
e_state <- vector()
for(i in 1:31){
  print("-----------")
  print(i)
  print(unlist(lapply(1:n_state, function(x) res_I_t$par[paste0("state_t[",i,",",x,"]")])))
  print(sum(unlist(lapply(1:n_state, function(x) res_I_t$par[paste0("state_t[",i,",",x,"]")]))))
  e_state[i] <- as.vector(1:n_state) %*% as.vector(unlist(lapply(1:n_state, function(x) res_I_t$par[paste0("state_t[",i,",",x,"]")])))
  #e_state[i] <- which.max(as.vector(unlist(lapply(1:5, function(x) res_I_t$par[paste0("state_t[",i,",",x,"]")]))))
}

# 결정변수를 시점이 아닌 주기로. if (i == I_t[1]|| i == I_t[2]|| i == I_t[3]) 을 era_1,2,3으로 나눠 각 구간에서 빈도를 결정변수로 바꾸는 방법시도

I_t <- sort(unlist(lapply(1:10, function(x) res_I_t$par[paste0("I_t[",x,"]")])))
ggplot() + aes(x=xvals, y=e_state) + geom_line() + geom_point(color="red") + ggtitle(paste("Expected state at time t target:", res_I_t$value)) +
  geom_vline(xintercept = I_t, linetype="dashed") + ylim(1, n_state) + geom_vline(xintercept = 20, linetype="solid")

for(i in 1:2){
  print(which.max(unlist(lapply(1:10, function(x) res_I_t$par[paste0("interval[",i,",",x,"]")]))))
}

comb_model <- stan_model(file.path(getwd(), "optimize/step1/inspection_interval_marginalize.stan"), verbose = FALSE)

comb_res <- rstan::sampling(comb_model, opt_data_I_t, algorithm="Fixed_param", chains=1)

summary(comb_res, c("cost"))$summary
