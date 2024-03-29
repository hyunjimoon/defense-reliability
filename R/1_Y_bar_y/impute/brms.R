library(cmdstanr); library(parallel); library("tidyverse"); library(mice); library(brms)
set.seed(1954)
.libPaths("~/Rlib")
options(mc.cores = parallel::detectCores())
set_cmdstan_path("/Users/hyunjimoon/Dropbox/20_paper/charles/code/cmdstan")

scriptDir <- getwd()
submodel <- "lgm"
modelDir <- file.path(scriptDir, submodel, "models")
dataDir <- file.path(scriptDir, "data")


gp_fit <- function(modelName, data){
  chains <- 2
  parallel_chains <- min(chains, detectCores())
  scriptDir <- getwd()
  delivDir <- file.path(scriptDir,submodel, "deliv", modelName)
  prefit <- file.path(delivDir, paste0(modelName, ".rda"))
  stanfile <- file.path(modelDir, modelName, paste0(modelName, ".stan"))
  
  if (file.exists(prefit)){
    fit <- readRDS(prefit)
  }else{ 
    mod <- cmdstan_model(stanfile, quiet = FALSE)
    fit <- mod$sample(data, chains = chains, iter_warmup = 250, iter_sampling = 250,
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
  
  plot(nondiv_samples$length_engine_mu, nondiv_samples$sigma_error_ship,
       col=c_dark_trans, pch=16, cex=0.8, xlab="length_engine_mu", ylab="sigma_error_ship")
  points(div_samples$length_engine_mu, div_samples$sigma_error_ship,
         col=c_green_trans, pch=16, cex=0.8)
  
  plot(nondiv_samples$length_GP_engine, nondiv_samples$length_engine_mu,
       col=c_dark_trans, pch=16, cex=0.8, xlab="length_GP_engine", ylab="length_engine_mu")
  points(div_samples$length_GP_engine, div_samples$length_engine_mu,
         col=c_green_trans, pch=16, cex=0.8)
}

mseNplot <- function(x, y_ext, ext_ind){
  yhat_imputed <- x %>%
    dplyr::filter(str_detect(variable, "y_pred")) %>%
    pull(mean)
  yhat_ext <- yhat_imputed[ext_ind]
  plot(1, type="n",xlim=c(0,31),ylim = c(0,300),xlab="age",ylab="failure count")
  for (n in 1:length(y_ext)){
    points(age_ind[n],y_ext[n],col="black",pch=16)
    points(age_ind[n],yhat_ext[n],col="red")
  }
  mean((y_ext-yhat_ext)^2)
}
#######################################################
# data
#######################################################
N_engines <- 5 
N_ships <- 99
N_ages <- 31
N_ages_obs <- 31

ship_engine_ind <- read.csv(paste0(dataDir,"/engine_type1to5.csv"))$engine
#ship_engine_ind <- read.csv(paste0(dataDir,"/engine_type1to4.csv"))$engine
ship_ind <- read.csv(paste0(dataDir,"/ship_index.csv"))$ship
age_ind <- read.csv(paste0(dataDir,"/x_age.csv"))[,-1]
engine_ind <- ship_engine_ind[ship_ind]

#######################################################
# BEFORE IMPUTATION
#######################################################

# existing 653 data
y_ext_df <- read.csv(paste0(dataDir,"/y_count_original.csv"))[,-1]
y_ext <- as.array(as.matrix(y_ext_df))
y_ext <-  y_ext[!is.na(y_ext)]

data = list(y_ext = y_ext,ship_ind = ship_ind, engine_ind =engine_ind, age_ind = age_ind)
fit1 <- brm(y_ext ~ ship_ind + engine_ind + age_ind, data = data,family = poisson())
fit2 <- brm(y_ext ~   age_ind + (1 + engine_ind |ship_ind), data = data,family = poisson())
fit3 <- brm(y_ext ~   ship_ind + engine_ind + age_ind + (1 + age_ind|engine_ind) + (1 + age_ind|ship_ind), data = data,family = poisson())
#engine_ind과 ship_ind의 상관성 표현
fit4 <- brm(y_ext ~   ship_ind + engine_ind + age_ind + (1 + age_ind|engine_ind) + (1 + age_ind|ship_ind), data = data,family = poisson())

#######################################################
# AFTER IMPUTATION: imputed 31*99 data
#######################################################
y_imp <- round(read.csv(paste0(dataDir,"/y_imputed9931_inverse.csv"))[,-1])
age_mean = rep(NA, length(unique(age_ind)))
for (i in unique(age_ind)){
  age_mean[i] <-  mean(as.matrix(y_df[i,]), na.rm = TRUE) 
}

data <- list(n_coordinates = 1, n_ages = length(unique(age_ind)), n_ships = length(ship_engine_ind), y = y_ext_df, ye = age_mean, ship_engine_ind = ship_engine_ind)
data$alpha_mu_prior <- 0
data$alpha_sd_prior <- 1
data$rho_mu_prior <- 0
data$rho_sd_prior <- 1

modelName <- "pois_lgm_ela2"
lgm_ela_fit <- gp_fit(modelName, data)
lgm_ela_fs <- lgm_ela_fit$summary()
ext_ind = which(!is.na(y_ext_df))
mseNplot(lgm_ela_fs, y_ext, ext_ind)