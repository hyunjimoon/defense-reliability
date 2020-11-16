#setwd("C:/Users/serim/Documents/GitHub/reliability_prediction")
source(file.path(getwd(), "impute/mice_imputation.R"))
scriptDir <- getwd()



mice_imp <- generateMice()
imputed_data <- complete(mice_imp, 1)

n_state = 5

generate_state_matrix <- function(data, n){
  state<-cut(data, breaks=c(quantile(data,seq(0,1,length.out = n+1))),labels = 1:n, include.lowest=TRUE)
  state<-as.numeric(as.character(state))
  matrix(state,nrow=31)
}

state_matrix <- generate_state_matrix(imputed_data$y_data, n_state)

plot_ship_state <- function(ship_n){
  ship_state <- data.frame(y=state_matrix[, ship_n], imputed=is.na(mice_imp$data[mice_imp$data$ship_ind == ship_n, ]$y_data))
  ggplot(data=ship_state, aes(x=1:nrow(ship_state), y=y, color=imputed)) + geom_line(color="gray") + ggtitle(paste("shipnum=",ship_n, "imputed=", sum(ship_state$imputed))) +
    geom_point() + scale_color_discrete(name="imputed") + labs(x="age", y="state") + coord_fixed(ratio=2) + ylim(1, 5)
}

plot_ship_state(as.integer(runif(1, 1, 99)))

