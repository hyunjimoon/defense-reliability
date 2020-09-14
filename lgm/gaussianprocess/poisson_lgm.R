library(rstan); library(cmdstanr); library(parallel); library("tidyverse")
setwd("~/Dropbox/20_paper/failure_forecast_gp")
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

scriptDir <- getwd()
modelDir <- file.path(scriptDir, "models")
dataDir <- file.path(scriptDir, "data")
outDir <- file.path(scriptDir, "deliv", modelName)

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
  
  #다음 변수와  sigma_error_ship간 Pair plot으로 pathology발견?
  #cov_engine = cov_exp_quad(ages, sigma_GP_engine, length_GP_engine);
  #cov_ship = cov_exp_quad(ages, sigma_GP_ship, length_GP_ship);
  
  par(mfrow=c(1, 3))
  plot(nondiv_samples$length_GP_engine, nondiv_samples$sigma_GP_engine, log="xy",
       col=c_dark_trans, pch=16, cex=0.8, xlab="length_GP_engine", ylab="sigma_GP_engine")
  points(div_samples$length_GP_engine, div_samples$sigma_GP_engine,
         col=c_green_trans, pch=16, cex=0.8)
  
  plot(nondiv_samples$length_engine_mu, nondiv_samples$sigma_error_ship,
       col=c_dark_trans, pch=16, cex=0.8, xlab="length_engine_mu", ylab="sigma_error_ship")
  points(div_samples$length_engine_mu, div_samples$sigma_error_ship,
         col=c_green_trans, pch=16, cex=0.8)
  
  plot(nondiv_samples$length_GP_engine, nondiv_samples$length_engine_mu,
       col=c_dark_trans, pch=16, cex=0.8, xlab="length_GP_engine", ylab="length_engine_mu")
  points(div_samples$length_GP_engine, div_samples$length_engine_mu,
         col=c_green_trans, pch=16, cex=0.8)
}
#######################################################

#######################################################
modelName <- "disease_map"
scriptDir <- getwd()
modelDir <- file.path(scriptDir, "models")
dataDir <- file.path(scriptDir, "data")
outDir <- file.path(scriptDir, "deliv", modelName)

N_engines <- 5 
N_ships <- 99
N_ages <- 31
N_ages_obs <- 31

ship_engine_ind <- read.csv("data/engine.csv")$engine
ship_ind <- read.csv("data/ship_index.csv")$X0
age_ind <- read.csv("data/x_age.csv", header = FALSE)[,1]
y_df <- read.csv("data/failure_count.csv")[,-1]
engine_mean = rep(NA, length(unique(ship_engine_ind)))
for (i in unique(ship_engine_ind)){
  engine_mean[i] <-  mean(as.matrix(y_df[c(ship_engine_ind==i)]), na.rm = TRUE) 
}
y <- as.array(as.matrix(y_df))
y <-  y[!is.na(y)]
##########################################################################
data <- list(n_obs = length(y), x = age_ind, n_covariates = 1, y = y, ye = engine_mean[ship_engine_ind][ship_ind])

data$alpha_mu_prior <- 0
data$alpha_sd_prior <- 1
data$rho_mu_prior <- 0
data$rho_sd_prior <- 1
###########################################################################
lgm_fit <- gp_fit(modelName, data)
