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
# Fit model with weakly-informative priors
#######################################################
modelName <- "hier_gp_tw"
#modelName <- "hier_gp_5var"

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
y <- read.csv("data/y_count_pwr.csv", header = FALSE)[,1]

data <- list(N = length(y), N_engines=N_engines,N_ships = N_ships, N_ages= N_ages, N_ages_obs = N_ages_obs, 
             ship_engine_ind =ship_engine_ind, ship_ind = ship_ind, age_ind = age_ind, y=y)

mseNplot <- function(x, y){
  yhat<- x %>%
    filter(str_detect(variable, "obs_mu")) %>%
    pull(mean)
  # yhat_fill<- hier_gp_simple %>%
  #   filter(str_detect(variable, "y_new_pred")) %>%
  #   pull(mean)
  yhat_df<- data.frame(matrix(yhat)) #, nrow = 31, ncol = 99
  yhat_df$idu <- as.numeric(row.names(yhat_df))
  yhat_df.melt <- reshape2::melt(yhat_df, id.vars="idu")
  names(yhat_df.melt) <- c("idu", "ship", "fail")
  ggplot(data = yhat_df.melt, aes(x=idu, y=fail))+
    geom_point(aes(group = ship))
  mean((y-yhat)^2)
}

########FORECAST - to rstan NOT NEEDED
fit1 <- gp_fit(modelName, data)
mse(fit1$summary(), y)

### DIVERGENCE CHECK -to rstan NEEDED
modelName <- "hier_gp_twtw"
sf_twtw_6 <- gp_fit(modelName, data)
rsf_twtw_6  <- read_stan_csv(sf_twtw_6$output_files())
util$check_all_diagnostics(rsf_twtw_6)
div_detect(rsf_twtw_6)

modelName <- "twtw_6_10"
sf_twtw_6_10 <- gp_fit(modelName, data)
rsf_twtw_6_10  <- read_stan_csv(sf_twtw_6_10$output_files())
util$check_all_diagnostics(rsf_twtw_6_10)
div_detect(rsf_twtw_6_10)

modelName <- "hier_gp_twtw"
fit_twtw <- gp_fit(modelName, data)
stanfit <- read_stan_csv(fit_twtw$output_files())
util$check_all_diagnostics(stanfit)

modelName <- "twtw_engineVar"
sf_twtw_engineVar <- gp_fit(modelName, data)
#forecast
mseNplot(sf_twtw_engineVar$summary(), y)
#divg check
rsf_twtw_engineVar  <- read_stan_csv(sf_twtw_engineVar$output_files())
util$check_all_diagnostics(rsf_twtw_engineVar)
div_detect(rsf_twtw_engineVar)


###############0810
#stan(file='models/hier_gp_weak/hier_gp_weak.stan', data=data, seed=5838298)
fit <-fail_fit
stanfit <- read_stan_csv(fit$output_files())
util$check_all_diagnostics(stanfit)

# A few extra divergences -- what's up?
# Let's look at some marginal posteriors


# diagnostics
check_all_diagnostics(fit)
pdf(file = file.path(outDir, paste(modelName,"Plots%03d.pdf", sep = "")), width = 6, height = 6, onefile = F)
parms <- c("alpha", "rho", "lp__")
mcmcHistory(fit, parms)
mcmcDensity(fit, parms)
dev.off()