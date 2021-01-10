setwd("C:/Users/serim/Documents/GitHub/reliability_prediction")

source(file.path(getwd(), "impute/mice_imputation.R"))
scriptDir <- getwd()
library(rstan)
library(ggplot2)
library(scales)

model_DMsep <- stan_model(file.path(getwd(), "optimize/step1/DMsep_5param.stan"), verbose = FALSE) #approx_deterioration_matrix

mice_imp <- generateMice()
imputed_data<- complete(mice_imp, 1)

imputed_data_gp <- read.csv("data/y_pred_5var.csv")[,-1]
imputed_data$y_data <- unlist(c(imputed_data_gp))


#################################### DM_sep
# original policy (wihtout pm)
n_state = 3
initial_state = 1

generate_state_matrix <- function(data, n){
  state <- cut(data, breaks=quantile(data,c(0,1/3,2/3,1)), labels=1:n, include.lowest = TRUE)
  state<-as.numeric(state)
  matrix(state,nrow=31)
}

state_matrix <- generate_state_matrix(imputed_data$y_data, n_state)
states <- as.vector(t(state_matrix)) #[imputed_data$engine_ind == engine_type] shiptype -> year
onehot <- list()
# one-hot encode per-data state to vector
for(i in 1:length(states)){
  t_tmp <- rep(0, len=n_state)
  t_tmp[states[i]] <- 1
  onehot[[i]] <- t_tmp
}

onehot_array <- aperm(array(unlist(onehot),c(n_state, length(onehot)))) # array() fills column-wise first
set.seed(210106)
test_ship_ind=sort(sample(1:99,5))
test_ind=c(sapply(test_ship_ind,function(i) i*31+(1:31)))
train_data <- onehot_array[-test_ind, ]
test_data <- onehot_array[test_ind, ]
opt_data <- list(N= dim(train_data)[1], n_state=n_state, state_obs=train_data, time_obs=imputed_data$age_ind[-test_ind], initial_state=initial_state)

res <- optimizing(model_DMsep, opt_data, verbose = TRUE)
############################## 
##############################
# Debug code

rate<-matrix(0,nrow=4,ncol=3)

for (era in 1:4){
  D <- matrix(as.vector(unlist(lapply(1:n_state, function(row){lapply(1:n_state, function(col){res$par[paste0("D[",era,",",row,",",col,"]")]})}))), nrow=3, byrow=T)
  for(j in 1:3){
    rate[era,j] <- res$par[paste0("rate[", era,",", j,"]")]
  }
  print(paste0("Era:",era))
  print(D)
}
print(rate)

predicted_state<-matrix(0,nrow=31,ncol=3)
DM_pow<-array(0,dim=c(31,3,3))

for (t in 1:31){
  for (i in 1:3){
    for (j in 1:3){
      DM_pow[t,i,j]=res$par[paste0("DM_pow[",t,",",i,",", j,"]")]
    }
  }
  predicted_state[t,]=DM_pow[t,,]%*%c(1,0,0)
}

observed_count = sapply(1:3,function(i) apply(state_matrix,1,function(x) sum(x==i)))
observed_prop=t(apply(observed_count,1,function(x) x/sum(x)))
observed_scaled_state = exp(exp(exp(observed_prop)))/10
observed_scaled_state = t(apply(observed_scaled_state,1,function(x) x/sum(x)))
scaled_state=exp(exp(exp(predicted_state)))/10
scaled_state =  t(apply(scaled_state,1,function(x) x/sum(x)))

#plot(1,type="n",xlab="year",ylab="state",pch=19,cex=1,xlim=c(0,32),ylim=c(.5,3.5))


#for(i in rep(1:31,each=3)){
#  for (j in rep(1:3,31)){
#    points(i,j,pch=16,cex=4,col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
#  }
#}

ship_observed_state = data.frame(ship_id=rep(1:99,31),time=factor(rep(1:31,each=99)),state=states)
total_count_df=data.frame(time=factor(rep(1:31,3)),states=factor(rep(1:3,each=31)),observed=c(observed_count),predicted=c(predicted_state)*99)

legend_size <- c(10,8,5,1,5,8,10)


ggplot(total_count_df, aes(x = time, y = states)) +
  geom_point(aes(colour="observed", size = observed),alpha=1) +
  geom_point(aes(colour="predicted", size = predicted),alpha=0.4) +
  scale_size(breaks = c(-3:3), range = c(1,30)) +
  guides(
    size=guide_legend(override.aes = list(size = legend_size))
  ) 

ggplot(total_count_df, aes(x = time, y = states)) +
  geom_point(aes(colour="1", size = observed)) +
  scale_colour_continuous(low = "blue", high = "red", breaks = c(-3:3) ) +
  scale_size(breaks = c(-3:3), range = c(1,15)) +
  guides(
    color= guide_legend(), 
    size=guide_legend(override.aes = list(size = legend_size))
  ) 



