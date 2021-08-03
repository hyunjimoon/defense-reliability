#setwd("C:/Users/serim/Documents/GitHub/reliability_prediction")
source(file.path(getwd(), "Y_bar_y/impute/mice_imputation.R"))
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
  ggplot(data=ship_state, aes(x=1:nrow(ship_state), y=y, color=imputed)) + geom_line(color="gray") +
    ggtitle(paste("shipnum=",ship_n, "imputed=", sum(ship_state$imputed))) +
    geom_point() + scale_color_discrete(name="imputed") + labs(x="age", y="state") + coord_fixed(ratio=2) + ylim(1, 5)
}

plot_ship_state <- function(ship_n){
  df <- data.frame(y = numeric(0), imputed = character(0), ship = numeric(0))
  for (i in 1:ship_n){
    ship_state <- data.frame(y=state_matrix[, i], imputed=is.na(mice_imp$data[mice_imp$data$ship_ind == i, ]$y_data))
    ship_state$id = i
    df <- rbind(df, ship_state)
  }
  df
}

tot_ship_state <- plot_ship_state(99)

plot_ship_state <- function(ship_n){
  ggplot(data=tot_ship_state, aes(x=1:nrow(ship_state), y=y, color=imputed)) + geom_line(color="gray") +
    ggtitle(paste("shipnum=",ship_n, "imputed=", sum(ship_state$imputed))) +
    geom_point() + scale_color_discrete(name="imputed") + labs(x="age", y="state") + coord_fixed(ratio=2) + ylim(1, 5) + facet_wrap(~name)
}

ship_state <- pivot_longer(as_tibble(state_matrix), cols  = names(as_tibble(state_matrix)))
ship_state <- as_tibble(state_matrix)
p <- ggplot(data=ship_state, aes(x=1:nrow(state_matrix), y=value, color=imputed)) + geom_line(color="gray") +
  geom_point()

p + facet_wrap(~name, ncol = nrow(state_matrix))



+ scale_color_discrete(name="imputed") + labs(x="age", y="state") + coord_fixed(ratio=2) + ylim(1, 5)

ggtitle(paste("shipnum=",ship_n, "imputed=", sum(ship_state$imputed))) +

+ facet_wrap(~name)

plot_ship_state(as.integer(runif(1, 1, 99)))
base <- ggplot(mpg2, aes(displ, hwy)) +
  geom_blank() +
  xlab(NULL) +
  ylab(NULL)

base + facet_wrap(~class, ncol = 3)

par(mfrow=c(3,3))
for (i in 1:9){
  print(i)
  plot_ship_state(i)
}
a <- as_tibble(state_matrix)

melt(a, id.vars = n_state)
tibble(state_matrix)



