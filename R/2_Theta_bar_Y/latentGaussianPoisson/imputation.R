scriptDir <- "C:/Users/serim/Documents/academic/Bayes_Study/reliability_prediction"
submodel <- "lgm"
modelDir <- file.path(scriptDir, submodel, "models")
dataDir <- file.path(scriptDir, "data")

y_ext_df <- read.csv(paste0(dataDir,"/y_count_original.csv"))[,-1]
dim(y_ext_df)
y_ext <- as.vector(as.matrix(y_ext_df))

# 1. imputeTS package

#install.packages("imputeTS")
library(imputeTS)

y_ts<-ts(y_ext)

y_ts_imp<-na_interpolation(y_ts)

ggplot_na_imputations(y_ts, y_ts_imp)

y_ts_imp

y_imp<-matrix(y_ts_imp,nrow=31)

write.csv(y_imp, "C:/Users/serim/Documents/academic/Bayes_Study/reliability_prediction/data/y_imputeTS.csv")



# 2. brms package

install.packages("brms")


data = list(y_ext = y_ext,ship_ind = ship_ind, engine_ind =engine_ind, age_ind = age_ind)
fit1 <- brm(y_ext ~ ship_ind + engine_ind + age_ind, data = data,family = poisson())
fit2 <- brm(y_ext ~ ship_ind + engine_ind + age_ind, data = data,family = poisson())