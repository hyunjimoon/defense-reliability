library(rstan); library(cmdstanr); library(parallel); library("tidyverse")
setwd("~/Dropbox/20_paper/failure_forecast_gp")

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
    dir.create("test")
    fit$save_object(file = prefit)
  }
  fit
}

#######################################################
# Fit model with weakly-informative priors
#######################################################
modelName <- "hier_gp_weak"
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
 

fit1 <- gp_fit(modelName, data)
hier_gp1 <- fit1$summary()

# fit2 <- gp_fit(modelName, data)
# hier_gp2 <- fit2$summary()
# yhat <- function(modelName, y){
#   fit <- gp_fit(modelName, data)
#   hier_gp <- fit$summary()
# }
rmse <- function(x, y){
  yhat<- x %>%
    filter(str_detect(variable, "obs_mu")) %>%
    pull(mean)d
  mean((y-yhat)^2)
}
rmse(hier_gp1, y)
rmse(hier_gp2, y)


###############0810
#stan(file='models/hier_gp_weak/hier_gp_weak.stan', data=data, seed=5838298)
stanfit <- read_stan_csv(fit$output_files())
util$check_all_diagnostics(stanfit)

# A few extra divergences -- what's up?
# Let's look at some marginal posteriors
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

plot(nondiv_samples$length_GP_engine, nondiv_samples$sigma_error_ship,
     col=c_dark_trans, pch=16, cex=0.8, xlab="length_GP_engine", ylab="sigma_error_ship")
points(div_samples$length_GP_engine, div_samples$sigma_error_ship,
       col=c_green_trans, pch=16, cex=0.8)

plot(nondiv_samples$sigma_GP_engine, nondiv_samples$sigma_error_ship,
     col=c_dark_trans, pch=16, cex=0.8, xlab="sigma_GP_engine", ylab="sigma_error_ship")
points(div_samples$sigma_GP_engine, div_samples$sigma_error_ship,
       col=c_green_trans, pch=16, cex=0.8)


#hier_gp <- gp_summary("hier_gp", data)
hier_gp <- fit$summary()

yhat<- hier_gp %>%
  filter(str_detect(variable, "y_new_pred")) %>%
  pull(mean)
yhat<- data.frame(matrix(yhat, nrow = 31, ncol = 99))
yhat$idu <- as.numeric(row.names(yhat))
yhat.melt <- reshape2::melt(yhat, id.vars="idu")
names(yhat.melt) <- c("idu", "ship", "fail")
ggplot(data = yhat.melt, aes(x=idu, y=fail))+ 
  geom_point(aes(group = ship))


yhat2_df<- data.frame(matrix(yhat2, nrow = 31, ncol = 99))
yhat2$idu <- as.numeric(row.names(yhat2_df))
yhat2.melt <- reshape2::melt(yhat2, id.vars="idu")
names(yhat2.melt) <- c("idu", "ship", "fail")
ggplot(data = yhat2.melt, aes(x=idu, y=fail))+ 
  geom_point(aes(group = ship))



#TODO: compare two, ask what count is

# diagnostics
check_all_diagnostics(fit)
pdf(file = file.path(outDir, paste(modelName,"Plots%03d.pdf", sep = "")), width = 6, height = 6, onefile = F)
parms <- c("alpha", "rho", "lp__")
mcmcHistory(fit, parms)
mcmcDensity(fit, parms)
dev.off()