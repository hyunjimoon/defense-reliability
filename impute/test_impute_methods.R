library(brms)
library(rstan)
source(file.path(getwd(), "impute/mice_imputation.R"))

MSE <- function(obs, pred){
  return (mean((obs - pred)^2))
  }

scriptDir <- getwd()
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

#model1 <- brm(y_data ~ ship_ind + engine_ind + age_ind, data=mice_data1, family=poisson(), cores=4)

#MSE(mice_data1$y_data, predict(model1))

#model2 <- brm(y_data ~ ship_ind + engine_ind + age_ind + (1 + age_ind|engine_ind) + (1 + age_ind|ship_ind), data = mice_data2,family = poisson(), cores=4)

#MSE(mice_data1$y_data, predict(model2))

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

  stanfit_data = list(
  N=length(mice_data1$y_data),
  engine_types=max(mice_data1$engine_ind),
  y=mice_data1$y_data,
  complexity=complexity_ind,
  age=(mice_data1$age_ind-1)/30,
  engine_type=mice_data1$engine_ind,
  relative_displacement=displacement_ind,
  engine_count=rep(2, length(mice_data1$y_data)),
  ship_number=mice_data1$ship_ind,
  ship_number_max=max(mice_data1$ship_ind),
  
  N_pred=length(mice_data1$y_data),
  age_pred=(mice_data1$age_ind-1)/30,
  complexity_pred=complexity_ind,
  relative_displacement_pred=displacement_ind,
  engine_count_pred=rep(2, length(mice_data1$y_data)),
  engine_type_pred=mice_data1$engine_ind,
  ship_number_pred=mice_data1$ship_ind
)

poisson_model <- stan(paste0(getwd(), "/poisson/models/fit_data_post_pred-iter2_backup(1012).stan"), data=stanfit_data, cores=1)