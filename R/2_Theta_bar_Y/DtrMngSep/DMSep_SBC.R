library(SBC)
library(mice)
library(cmdstanr)
library(dplyr)
source(file.path(getwd(), "Y_bar_y/impute/engine_fail_imputation.R"))
source("tools/functions.r")

# Three simulators for prior, data, posterior
# 1. prior simulator
true <- as_tibble(read.csv("./data/DMSep_validation_trueparam.csv"))
# 2. data simulator
# DGP with truth-point benchmarking
p <- true$x[true$X== "p"]
q <- true$x[true$X== "q"]
rate <- matrix(nrow = n_period, ncol = n_state)
for (i in 1:n_period) for (j in 1:n_state) rate[i,j] <- true$x[true$X== sprintf("rate[%s,%s]",  i, j)]
true_pars <-list(p, q, rate)

generator <- function(){
  function(true_pars, iter){
    # Metadata
    n_data = 3069
    n_state = 3
    n_period = 4
    n_age = 31
    init_state = 1
    init_state_vec = rep(0, n_state);
    init_state_vec[init_state] = 1
    time_obs  =  complete(generateMice(), 1)$age_ind # Q.  generateMice()를 generator내에서 실행희망시 input없는거 가능?

    # Hyperparameter -> receive as input
    p = true_pars[[iter]][1]
    q = true_pars[[iter]][2]
    rate = true_pars[[iter]][3]
    #rate <- matrix(nrow = n_period, ncol = n_state)
    #for (i in 1:n_period) for (j in 1:n_state) rate[i,j] <- true$x[true$X== sprintf("rate[%s,%s]",  i, j)]

    # Maintenance
    M <- matrix(nrow = n_state, ncol = n_state)
    M[1,1]=1
    M[1,2]=p
    M[1,3]=q
    M[2,1]=0
    M[2,2]=(1-p)
    M[2,3]= (1-q)
    M[3,1]=0
    M[3,2]=0
    M[3,3]=0

    # Deterioration
    D <- array(rep(1, n_state * n_state * n_period), dim=c(n_state, n_state, n_period))
    DM_pow <- array(rep(1, n_state * n_state * n_age), dim=c(n_state, n_state, n_age))
    for(i in 1:n_period){
      D[1,1,i] = exp(-(rate[i,1]+ rate[i,2]))
      D[2,1,i] = rate[i,1] * exp(-rate[i,3]) * (1-exp(-(rate[i,1]+ rate[i,2] - rate[i,3]))) / (rate[i,1]+ rate[i,2] - rate[i,3])
      D[1,2,i] = 0
      D[3,1,i] = 1 - D[1,1,i] - D[2,1,i]
      D[2,2,i] = exp(-rate[i,3])
      D[3,2,i] = 1 - D[2,2,i]
      D[1,3,i] = 0
      D[2,3,i] = 0
      D[3,3,i] = 1
    }
    DM_pow[,,1] = D[,,1];
    for (i in 2:n_age){
      if (i <= 8){ DM_pow[,,i] = D[,,1] * M * DM_pow[,,i-1];} # colsum is 1 for right stoch.matrix
      else if (i <=20){ DM_pow[,,i] = D[,,2] * M * DM_pow[,,i-1];}
      else if (i <=26){DM_pow[,,i] = D[,,3] * M * DM_pow[,,i-1];}
      else{DM_pow[,,i] = D[,,4] * M * DM_pow[,,i-1];}
    }
  list(
    generated = lapply(seq(1:length(time_obs)), function(i){DM_pow[,,time_obs[i]] %*% init_state_vec}),
    parameters = list(
      p = p,
      q = q,
      rate = rate
    )
  )
  }
}

# 3. posterior simulator: prior and data simulator embedded in stan with supplied data
# 3-1. prior and data simulator
modelName = "DMSep"
Dir <- set_get_Dir(modelName,  "./optimize")
DMSep = cmdstanr::cmdstan_model(Dir$file)

# 3-2. simulated y data template with init as real y
n_data = 3069
n_state = 3
n_period = 4
n_age = 31
init_state = 1
init_state_vec = rep(0, n_state);
init_state_vec[init_state] = 1
mice_imp <- generateMice()
time_obs  =  complete(mice_imp, 1)$age_ind
generate_state_matrix <- function(data, n){
  state <- cut(data, breaks=quantile(data,c(0,1/3,2/3,1)), labels=1:n, include.lowest = TRUE)
  state<-as.numeric(state)
  matrix(state,nrow=31)
}
imputed_data<- complete(mice_imp, 1)
state_matrix <- generate_state_matrix(imputed_data$y_data, n_state)
states <- as.vector(t(state_matrix))
onehot <- list()
for(i in 1:n_data){
  t_tmp <- rep(0, len=n_state)
  t_tmp[states[i]] <- 1
  onehot[[i]] <- t_tmp
}
onehot_array <- aperm(array(unlist(onehot),c(n_state, length(onehot)))) # array() fills column-wise first
data_tpl <- list(n_data = dim(onehot_array)[1], n_state=n_state, n_period = n_period, n_age = n_age, state_obs=onehot_array,
     time_obs=imputed_data$age_ind, initial_state=initial_state)
N = 10
M = 10 # 40
workflow <- SBCWorkflow$new(DMSep, generator())
workflow$simulate(N, true_pars)
workflow$simulate(N, iter_tp)
workflow$fit_model(sample_iterations = M, warmup_iterations = M, data_tpl)
prior <- workflow$prior_samples
post <- workflow$posterior_samples
workflow$calculate_rank()
plot_ecdf(workflow, var = "theta")
plot_ecdf_diff(workflow, var="theta")
