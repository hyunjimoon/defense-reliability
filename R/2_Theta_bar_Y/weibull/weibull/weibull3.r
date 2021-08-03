library(rstan)
library(parallel)
library(tidyverse)
setwd("C:/Users/serim/Documents/academic/Bayes_Study/failure_forcast_weibull")

N<-sum(apply(count[,-1], 2, function(x) {sum(x[!is.na(x)])})) #28013
N_engines <- 5
N_ships <- 99
N_ages <- 31
ship_engine_ind<-read.csv("data/engine_index.csv")$engine
ship_engine_ind

count<-read.csv("data/failure_count.csv")

ship_n<-apply(count[,-1], 2, function(x) {sum(x[!is.na(x)])})
ship_index<-c()
time<-c()
for (s in 1:N_ships){ ship_index<-c(ship_index,rep(s,ship_n[s]))}
for (s in 1:N_ships){
  for (t in 1:N_ages){
    if (!is.na(count[t,s+1])) {
      time<-c(time,rep(t,count[t,s+1]))
    }
  }
}

weibull3_data<-list(N=N,N_engines=N_engines,N_ships=N_ships,N_ages=N_ages,
                   ship_engine_ind=ship_engine_ind,
                   ship_n=ship_n,ship_index=ship_index,
                   time=time
)

time

weibull3_model <- stan_model("models/weibull/weibull3.stan", auto_write = TRUE)
fit_weibull3 <- sampling(weibull3_model,
                        data = weibull3_data,
                        iter = 1000,
                        chains = 4)

weibull3_pars=c('mu_alpha_bar','mu_sigma_bar', "alpha_bar", "sigma_bar")

print(fit_weibull3, pars = weibull3_pars)

plot(fit_weibull3, plotfun = "trace", pars = weibull3_pars)

fit_weibull4 <- sampling(weibull3_model,
                         data = weibull3_data,
                         iter = 1000,
                         chains = 4)

weibull3_pars=c('mu_alpha_bar','mu_sigma_bar', "alpha_bar", "sigma_bar")

print(fit_weibull4, pars = weibull3_pars)

plot(fit_weibull4, plotfun = "trace", pars = weibull3_pars)

time_pred<-summary(fit_weibull4, pars = "time_pred", probs = 0.5)$summary[,4]

alpha_sigma<-matrix(summary(fit_weibull4, pars = c("alpha_bar", "sigma_bar"), probs = 0.5)$summary[,4],ncol=2)

pred<-matrix(0,nrow=N_ages,ncol=N_ships)




par(mfrow = c(2, 3))

for (e in 1:N_engines){
  plot(seq(0,31,length=10000),dweibull(seq(0,31,length=10000),alpha_sigma[e,1],alpha_sigma[e,2]),type='l',xlab="",ylab="")
  for (i in 1:N_ships){
    if (ship_engine_ind[i]==e){
      points((1:31)[!is.na(count[,i+1])],count[,i+1][!is.na(count[,i+1])]/sum(count[,i+1][!is.na(count[,i+1])]),col='red')
    }
  }
}





