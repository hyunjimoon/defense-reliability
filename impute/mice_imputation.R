library(mice)

generateMice <- function(){
  scriptDir <- getwd()
  dataDir <- dataDir <- file.path(scriptDir, "data")
  
  ship_engine_ind <- read.csv(paste0(dataDir,"/engine_type1to5.csv"))$engine
  ship_ind <- read.csv(paste0(dataDir,"/ship_index.csv"))$ship
  age_ind <- read.csv(paste0(dataDir,"/x_age.csv"))[,-1]
  engine_ind <- ship_engine_ind[ship_ind]
  
  y_ext_df <- read.csv(paste0(dataDir,"/y_count_original.csv"))[,-1]
  y_ext <- as.array(as.matrix(y_ext_df))
  
  
  age_ind <- array(dim=31*99)
  ship_ind <- array(dim=31*99)
  
  for(shiptype in 1:99){
    for(date in 1:31){
      ship_ind[(shiptype-1)*31 + date] <- shiptype
      age_ind[(shiptype-1)*31 + date] <- date
    }
  }
  engine_ind = array(dim=31*99)
  for(i in 1:length(engine_ind)){
    if(i %% 99 == 0){
      engine_ind[i] <- ship_engine_ind[99]
      next
    } 
    engine_ind[i] <- ship_engine_ind[i %% 99]
  }
  y_data <- as.vector(y_ext)
  data_df <- data.frame(y_data, age_ind, ship_ind, engine_ind)
  #md.pattern(data_df)
  
  mice_imputed <- mice(data_df, m=2)
  
  return(mice_imputed)
}

if(interactive()){
  #generateMice()
  
}