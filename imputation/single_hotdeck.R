imputeFromDF <- function(full_y_df, use_row=TRUE) {
  for(i in 1:dim(missing)[1]){
    row <- missing[i, 1]
    col <- missing[i, 2]
    if(use_row){
      impute_candidate = full_y_df[row, ]
    }
    else{
      impute_candidate <- full_y_df[, colnames(full_y_df)[col]]  
    }
    impute_candidate <- unique(impute_candidate[!is.na(impute_candidate)])
    full_y_df[row, col] <- sample(impute_candidate, 1, replace=TRUE)
  }
  return(full_y_df)
}


if(interactive()){
  N_engines <- 5 
  N_ships <- 99
  N_ages <- 31
  N_ages_obs <- 31
  
  scriptDir <- getwd()
  dataDir <- dataDir <- file.path(scriptDir, "data")
  
  ship_engine_ind <- read.csv(paste0(dataDir,"/engine_type1to4.csv"))$engine
  ship_ind <- read.csv(paste0(dataDir,"/ship_index.csv"))$ship
  age_ind <- read.csv(paste0(dataDir,"/x_age.csv"))[,-1]
  engine_ind <- ship_engine_ind[ship_ind]
  
  y_ext_df <- read.csv(paste0(dataDir,"/y_count_original.csv"))[,-1]
  y_ext <- as.array(as.matrix(y_ext_df))
  n1_y_ext <-  y_ext[!is.na(y_ext)]
  imputed <- imputeFromDF(y_ext, use_row = FALSE)
  dim(imputed)
  any(is.na(imputed))
  
}


