library(rstan)
library(parallel)
library(tidyverse)
setwd("C:/Users/serim/Documents/academic/Bayes_Study/failure_forcast_weibull")

N<-653
N_engines <- 5
N_ships <- 99
N_ages <- 31
ship_engine_ind<-read.csv("data/engine_index.csv")$engine
count<-read.csv("data/failure_count.csv")
N_missing<-apply(count[,-1], 2, function(x) {sum(is.na(x))})
mean_failure<-apply(count[,-1], 2, function(x) {mean(x[!is.na(x)])})
mean_time<-apply(count[,-1], 2, function(x) {mean((1:31)[!is.na(x)])})


count_nona<-count
F_count<-count

for (ship in 2:(N_ships+1)){
  count_nona[,ship][is.na(count[,ship])]<-mean(count[,ship][!is.na(count[,ship])])
  F_count[,ship]<-(cumsum(count_nona[,ship])-0.3)/(sum(count_nona[,ship])+0.4)
  F_count[,ship][is.na(count[1:t,ship])]<-NA
}

y<-F_count[,-1]
y[is.na(y)]<-0

##plot

plot(y[,1],ylim=c(0,1),type='l',col=1)
for (ship in 2:(N_ships)){
  if (ship_engine_ind[ship]==5){
  lines(y[,ship],ylim=c(0,1),type='l',col=ship_engine_ind[ship])}
}

NA_ind<-apply(count[,-1], 2, function(x) ifelse(is.na(x),1,0))
NA_ind<-matrix(NA_ind,nrow=31)

weibull_data<-list(N=N,N_engines=N_engines,N_ships=N_ships,N_ages=N_ages,
                   ship_engine_ind=ship_engine_ind,
                   N_missing=N_missing,mean_failure=mean_failure,
                   mean_time=mean_time,y=y,NA_ind=NA_ind
                   )
  
weibull_model <- stan_model("models/weibull/weibull2.stan", auto_write = TRUE)
fit_weibull <- sampling(weibull_model,
                            data = weibull_data,
                            iter = 2000,
                            chains = 4)

weibull_pars=c('mu_alpha_bar','mu_sigma_bar', "alpha_bar", "sigma_bar")

print(fit_weibull, pars = weibull_pars)

plot(fit_weibull, plotfun = "trace", pars = weibull_pars)

