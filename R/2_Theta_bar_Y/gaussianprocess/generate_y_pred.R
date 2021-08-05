library(rstan); library(cmdstanr); library(parallel); library("tidyverse"); library(dplyr)
setwd("C:/Users/serim/Documents/academic/Bayes_Study/reliability_prediction/gaussianprocess")
set.seed(1954)
.libPaths("~/Rlib")
options(mc.cores = parallel::detectCores())
util <- new.env()
source('stan_utility.R', local=util)
source('gp_utility.R', local=util)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_trans <- c("#DCBCBC80")
c_light_highlight_trans <- c("#C7999980")
c_mid_trans <- c("#B97C7C80")
c_mid_highlight_trans <- c("#A2505080")
c_dark_trans <- c("#8F272780")
c_dark_highlight_trans <- c("#7C000080")

c_light_teal="#6B8E8E"
c_mid_teal="#487575"
c_dark_teal="#1D4F4F"

c_green_trans <- c("#00FF0080")
c_superfine <- c("#8F272705")

println <- function(msg) cat(msg); cat("\n") 
printf <- function(pattern, ...) println(sprintf(pattern, ...)) 
print_file <- function(file) println(readLines(file))


scriptDir <- getwd()
modelDir <- file.path(scriptDir, "models")
dataDir <- file.path(scriptDir, "data")
gp_fit <- function(modelName, data){
  chains <- 4
  parallel_chains <- min(chains, detectCores())
  scriptDir <- getwd()
  delivDir <- file.path(scriptDir, "deliv", modelName)
  prefit <- file.path(delivDir, paste0(modelName, ".rda"))
  stanfile <- file.path(modelDir, modelName, paste0(modelName, ".stan"))
  
  if (file.exists(prefit)){
    fit <- readRDS(prefit)
  }else{ 
    mod <- cmdstan_model(stanfile, quiet = FALSE)
    fit <- mod$sample(data, chains = chains, iter_warmup = 500, iter_sampling = 500,
                      parallel_chains = parallel_chains, save_warmup = FALSE)
    dir.create(delivDir)
    fit$save_object(file = prefit)
  }
  fit
}

div_detect <- function(stanfit){
  partition <- util$partition_div(stanfit)
  div_samples <- partition[[1]]
  nondiv_samples <- partition[[2]]
  
  par(mfrow=c(1, 3))
  plot(nondiv_samples$length_GP_engine, nondiv_samples$sigma_GP_engine, log="xy",
       col=c_dark_trans, pch=16, cex=0.8, xlab="length_GP_engine", ylab="sigma_GP_engine")
  points(div_samples$length_GP_engine, div_samples$sigma_GP_engine,
         col=c_green_trans, pch=16, cex=0.8)
  
  plot(nondiv_samples$length_GP_engine, nondiv_samples$sigma_error_ship,
       col=c_dark_trans, pch=16, cex=0.8, xlab="length_GP_engine", ylab="sigma_error_ship")
  points(div_samples$length_GP_engine, div_samples$sigma_error_ship,
         col=c_green_trans, pch=16, cex=0.8)
  
  plot(nondiv_samples$length_engine_scale, nondiv_samples$length_GP_engine_s, 
       col=c_dark_trans, pch=16, cex=0.8, xlab="length_engine_scale", ylab="length_GP_engine_s")
  points(div_samples$length_engine_scale, div_samples$length_GP_engine_s,
         col=c_green_trans, pch=16, cex=0.8)
}
```
data preparation
```{r}
N_engines <- 5 
N_ships <- 99
N_ages <- 31
N_ages_obs <- 31

ship_engine_ind <- read.csv("data/engine.csv")$engine
ship_ind <- read.csv("data/ship_index.csv")$X0
age_ind <- read.csv("data/x_age.csv", header = FALSE)[-1,1]
y <- read.csv("data/y_count_pwr.csv", header = FALSE)[,1]

data <- list(N = length(y), N_engines=N_engines,N_ships = N_ships, N_ages= N_ages, N_ages_obs = N_ages_obs, 
             ship_engine_ind =ship_engine_ind, ship_ind = ship_ind, age_ind = age_ind, y=y)


N_engines <- 5 
N_ships <- 99
N_ages <- 31
N_ages_obs <- 31

ship_engine_ind <- read.csv("data/engine.csv")$engine
ship_ind <- read.csv("data/ship_index.csv")$X0
age_ind <- read.csv("data/x_age.csv", header = FALSE)[-1,1]
y <- read.csv("data/y_count_pwr.csv", header = FALSE)[,1]

data <- list(N = length(y), N_engines=N_engines,N_ships = N_ships, N_ages= N_ages, N_ages_obs = N_ages_obs, 
             ship_engine_ind =ship_engine_ind, ship_ind = ship_ind, age_ind = age_ind, y=y)

mseNplot <- function(x, y){
  yhat<- x %>%
    dplyr::filter(str_detect(variable, "y_new_pred")) %>%
    pull(mean)
  
  yhat<- (matrix(yhat, nrow = 31, ncol = 99))
  y_hat <- rep(NA, length(y))
  
  for (i in 1:length(y)){
    y_hat[i] <- yhat[age_ind[i],ship_ind[i]]
  }
  
  par(mfrow=c(1,1))
  plot(1, type="n",xlim=c(0,31),ylim=c(-3,4),xlab="age",ylab="scaled failure count")
  for (n in 1:653){
    points(age_ind[n],y[n],col="black",pch=16)
    points(age_ind[n],y_hat[n],col="red")
  }
  
  mean((y-y_hat)^2)
}

   
data$hp_scale <- 1
data$emp_le_shape <- 18.707524
data$emp_le_scale <- 4.120272
data$emp_ls_shape <- 11.289554
data$emp_ls_scale <- 6.821292

modelName <- "gp_hp_n_5var"

sf_gp_hp_n_5var <- gp_fit(modelName, data)

sf_gp_hp_n_5var_sm<-sf_gp_hp_n_5var$summary()

yhat<- sf_gp_hp_n_5var_sm %>%
  dplyr::filter(str_detect(variable, "y_new_pred")) %>%
  pull(mean)

yhat<- (matrix(yhat, nrow = 31, ncol = 99))

write.csv(yhat,file="C:/Users/serim/Documents/academic/Bayes_Study/reliability_prediction/gaussianprocess/y_pred_5var.csv")

y_hat <- rep(NA, length(y))

for (i in 1:length(y)){
  y_hat[i] <- yhat[age_ind[i],ship_ind[i]]
}

par(mfrow=c(1,1))
plot(1, type="n",xlim=c(0,31),ylim=c(-3,4),xlab="age",ylab="scaled failure count")
for (n in 1:653){
  points(age_ind[n],y[n],col="black",pch=16)
  points(age_ind[n],y_hat[n],col="red")
}

mean((y-y_hat)^2)
