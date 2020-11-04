library(brms)
library(rstan)
library(loo)
library(ggplot2)
source(file.path(getwd(), "impute/mice_imputation.R"))

MSE <- function(obs, pred){
  return (mean((obs - pred)^2))
  }

scriptDir <- getwd()
###########################
# OUT OF SAMPLE TEST DATA


test_y <- read.csv(paste0(getwd(), "/impute/failure_count_test.csv"))[, 2:11]
test_y_data <- test_y[!is.na(test_y)]
test_ship <- which(!is.na(test_y), arr.ind=TRUE)[,"col"]
test_age <- which(!is.na(test_y), arr.ind=TRUE)[,"row"]
test_engine <- ship_engine_ind[test_ship]

cv_testdata <- data.frame(y_data=test_y_data, ship_ind=test_ship, engine_ind=test_engine, age_ind=test_age)


############################
# Prepare and impute data

dataDir <- dataDir <- file.path(scriptDir, "data")

ship_engine_ind <- read.csv(paste0(dataDir,"/engine_type1to4.csv"))$engine
ship_ind <- read.csv(paste0(dataDir,"/ship_index.csv"))$ship
age_ind <- read.csv(paste0(dataDir,"/x_age.csv"))[,-1]
engine_ind <- ship_engine_ind[ship_ind]

y_ground_df <- read.csv(paste0(dataDir,"/y_count_original.csv"))[,-1]
y_ground <- as.array(as.matrix(y_ground_df))

mice_imp <- generateMice()

mice_data1 <- complete(mice_imp, 1)
mice_data2 <- complete(mice_imp, 2)

ggplot(mice_data1, aes(mice_data1$ship_ind, mice_data1$age_ind, fill=mice_data1$y_data)) + scale_fill_gradientn(colours=rainbow(5)) + geom_tile() + coord_equal()
model1 <- brm(y_data ~ ship_ind + engine_ind + age_ind, data=mice_data1, family=poisson(), cores=4)
print("imputed model1 cv MSE:")
print(MSE(mice_data1$y_data, predict(model1, newdata = cv_testdata)))

model2 <- brm(y_data ~ ship_ind + engine_ind + age_ind + (1 + age_ind|engine_ind) + (1 + age_ind|ship_ind), data = mice_data1,family = poisson(), cores=4)
print("imputed model2 cv MSE:")
print(MSE(mice_data1$y_data, predict(model2, newdata = cv_testdata)))

complexity_ind = array(dim=31*99)
for(i in 1:length(complexity_ind)){
  if(i %% 99 >= 34 && i %% 99 <= 76){
    complexity_ind[i] <- 1.3
    next
  }
  else{
    complexity_ind[i] <- 1.0
  }
}

displacement <- c(5520, 3100, 170, 570, 3300)
displacement_ind = array(dim=31*99)
for(i in 1:length(displacement_ind)){
  idx <- i %% 99
  if(i %% 99 == 0){
    displacement_ind[i] <- round(3300 / 170)
  }
  else if(idx <= 7){
    displacement_ind[i] <- round(5520 / 170)
  }
  else if(idx <= 33){
    displacement_ind[i] <- round(3100 / 170)
  }
  else if(idx <= 76){
    displacement_ind[i] <- round(170 / 170)
  }
  else if(idx <= 95){
    displacement_ind[i] <- round(570 / 170)
  }
  else {
    displacement_ind[i] <- round(3300 / 170)
  }
}

# stanfit_data = list(
#   N=length(mice_data1$y_data),
#   engine_types=max(mice_data1$engine_ind),
#   y=mice_data1$y_data,
#   complexity=complexity_ind,
#   age=(mice_data1$age_ind-1)/30,
#   engine_type=mice_data1$engine_ind,
#   relative_displacement=displacement_ind,
#   engine_count=rep(2, length(mice_data1$y_data)),
#   ship_number=mice_data1$ship_ind,
#   ship_number_max=max(mice_data1$ship_ind),
#   
#   N_pred=length(mice_data1$y_data),
#   age_pred=(mice_data1$age_ind-1)/30,
#   complexity_pred=complexity_ind,
#   relative_displacement_pred=displacement_ind,
#   engine_count_pred=rep(2, length(mice_data1$y_data)),
#   engine_type_pred=mice_data1$engine_ind,
#   ship_number_pred=mice_data1$ship_ind
# )
# 
# poisson_model <- stan(paste0(getwd(), "/poisson/models/fit_data_post_pred.stan"), data=stanfit_data, cores=4)


###################################
N_engines <- 5 
N_ships <- 99
N_ages <- 31
N_ages_obs <- 31


ship_engine_ind <- read.csv(paste0(dataDir,"/engine_type1to4.csv"))$engine
gp_data <- list(N = length(mice_data1$y_data), 
             N_engines=max(ship_engine_ind),
             N_ships = max(mice_data1$ship_ind), 
             N_ages= max(mice_data1$age_ind), 
             N_ages_obs = max(mice_data1$age_ind), 
             ship_engine_ind = ship_engine_ind, 
             ship_ind = mice_data1$ship_ind, 
             age_ind = mice_data1$age_ind,
             y=mice_data1$y_data)

gp_data$hp_scale <- 3

gp_model <- stan(paste0(getwd(), "/gaussianprocess/models/gp_hp_lognorm/gp_hp_lognorm.stan"), data=gp_data, cores=4)

gp_loo <- loo(gp_model)

gp_pred = array(dim=31*99)
gp_summary <-summary(gp_model, pars=list("y_new_pred"))
for(t in 1:max(mice_data1$age_ind)){
  for(s in 1:max(mice_data1$ship_ind)){
   # print(format("y_new_pred[%d,%d]", t, s))
    gp_pred[(t-1) * 99 + s] <- gp_summary$summary[paste0("y_new_pred[", t,",", s,"]"), "mean"]
    
  }
}

###################################


# poisson_loo <- loo(poisson_model)
# model1_loo <- loo(model1)
# model2_loo = loo(model2)
# 
# #loo_compare(lost(model1_loo, poisson_loo))
# 
# rt <- loo_model_weights(list(poisson=poisson_loo, brms1=model1_loo, brms2=model2_loo))#, gp=gp_loo, brms1=model1_loo, brms2=model2_loo))
# 
# poisson_pred = rep(0, length(mice_data1$y_data))
# poisson_pred_ = extract(poisson_model, pars=list("y_post_pred"))$y_post_pred
# for(i in 1:length(mice_data1$y_data)){
#   poisson_pred[i] = mean(poisson_pred_[, i])
# }
# pooled = as.matrix(cbind(poisson=poisson_pred, brms1=predict(model1)[, "Estimate"], brms2=predict(model2)[, "Estimate"]))
# 
# MSE(mice_data1$y_data, pooled %*% as.vector(rt))
# MSE(mice_data1$y_data, poisson_pred)
# MSE(mice_data1$y_data, gp_pred)

###############################
# TEST OUT OF SAMPLE DATA


brms1_test_pred <- predict(model1, newdata = newdata)[, "Estimate"]
brms2_test_pred <- predict(model2, newdata = newdata)[, "Estimate"]

print("out of sample test data pred MSE")
print("1. brms1_test_pred")
MSE(test_y_data, brms1_test_pred)
print("2. brms2_test_pred")
MSE(test_y_data, brms2_test_pred)


gp_test_pred <- rep(0.0, length = length(test_y_data))
gp_summary <-summary(gp_model, pars=list("mu", "age_re", "engine_re", "ship_re", "GP_engine", "GP_ship"))$summary[, "mean"]
gp_mu <- gp_summary["mu"]
for(i in 1:length(test_y_data)){
  gp_test_pred[i] <- gp_mu + gp_summary[paste0("age_re[",test_age[i],"]")] + gp_summary[paste0("engine_re[",test_engine[i],"]")] + gp_summary[paste0("ship_re[",test_ship[i],"]")] + gp_summary[paste0("GP_engine[",test_age[i],",", test_engine[i], "]")] + gp_summary[paste0("GP_ship[",test_age[i],",", test_ship[i], "]")]
}

print("imputed gp model cv MSE:")
MSE(cv_testdata, gp_pred)
