source(file.path(getwd(), "R/1_Y_bar_y/impute/mice_imputation.R"))
scriptDir <- getwd()
set.seed(210106)
library(rstan)
library(cmdstanr)
library(ggplot2)
library(scales)

#cmdstan: model_DMsep <-stan_m(file = (file.path(getwd(), "R/2_Theta_bar_Y/DtrMngSep/models/DMSep/DMSep.stan")))
model_DMsep <-stan_model(file = (file.path(getwd(), "R/2_Theta_bar_Y/DtrMngSep/models/DMSep/DMSep.stan")))

mice_imp <- generateMice()
imputed_data<- complete(mice_imp, 1)

imputed_data_gp <- read.csv("data/y_pred_5var.csv")[,-1]
imputed_data$y_data <- unlist(c(imputed_data_gp))
#################################### DM_sep
# original policy (wihtout pm)
n_state = 3
initial_state = 1

generate_state_matrix <- function(data, n){
  #state <- cut(data, breaks=c(0, 80, 160, max(data)), labels=1:n, include.lowest = TRUE)
  state <- cut(data, breaks=quantile(data,c(0,1/3,2/3,1)), labels=1:n, include.lowest = TRUE)
  state<-as.numeric(state)
  matrix(state,nrow=31)
}

state_matrix <- generate_state_matrix(imputed_data$y_data, n_state)
states <- as.vector(t(state_matrix))

iter=2000
MSE_df<-data.frame(index=1:2*iter,p=rep(0,2*iter),q=rep(0,2*iter),train_MSE=rep(0,2*iter),test_MSE=rep(0,2*iter))
rate_array<-data.frame(period=c(),index=c(),rate=c())
D_array <- array(0, dim=c(iter*2,4,3,3))
ship_ind_df<-matrix(0,nrow=iter*2,ncol=5)

test_ship_ind= c(17,20,24,77,82) #sort(sample(1:99,5)) #c(17,20,24,77,82)
ship_ind_df[1,]=test_ship_ind
test_ind=c(sapply(test_ship_ind,function(x) (x-1)*31+(1:31)))

train_data <- states[-test_ind]
test_data <- states[test_ind]
stan_data <- list(N= length(train_data),T = max(imputed_data$age_ind[-test_ind]), S = 3, P = 4, states=train_data, obs2time=imputed_data$age_ind[-test_ind], initial_state=initial_state)
sampling_res<-sampling(model_DMsep,stan_data, iter = 2000, cores=4)
res_df<-as.data.frame(sampling_res)
                  
# deprecated - may change DM_pow to latent_states
sample_mean<-apply(res_df,2,mean) #write.csv(sample_mean, "./data/DMSep_validation_trueparam.csv") or sample_mean[1:338]

# load sample

res_df<-readRDS("R/2_Theta_bar_Y/DtrMngSep/sample.RDS")
res_df<-as.data.frame(res_df)
summary_df<-data.frame(variable=c("lambda1","lambda2","lambda3","p21","p31"),mean=rep(0,5),sd=rep(0,5))


# rate 1,2,3

par(mfrow=c(1,3))
#png(filename=paste0("R/2_Theta_bar_Y/DtrMngSep/figure/new_lambdas.png"))
for (i in 1:3){
  ymax=max(density(res_df[,i])$y)
  xmin=min(res_df[,i])
  xmax=max(res_df[,i])
  hist(-100,ylim=c(0,ymax),xlim=c(xmin,xmax),main=paste0("Lambda ",i),xlab="",freq=FALSE)
  summary_df[i,2]=mean(res_df[,i])
  summary_df[i,3]=sd(res_df[,i])
  lines(density(res_df[,i]),lwd=2)
}
#dev.off()

# M(p12,p13)

par(mfrow=c(1,2))

for (i in c("p21","p31")){
  png(filename=paste0("R/2_Theta_bar_Y/DtrMngSep/figure/M_",i,".png"))
  ymax=max(density(res_df[,i])$y)
  xmin=min(res_df[,i])
  xmax=max(res_df[,i])
  hist(-100,ylim=c(0,ymax),xlim=c(xmin,xmax),main=paste0(i),xlab="",freq=FALSE)
  lines(density(res_df[,i]),lwd=2)
  dev.off()
}

summary_df[4,2]=mean(res_df[,"p21"])
summary_df[4,3]=sd(res_df[,"p21"])
summary_df[5,2]=mean(res_df[,"p31"])
summary_df[5,3]=sd(res_df[,"p21"])

#write.csv(summary_df,"summary.csv")

# D

Dtr_summary<-data.frame(coord=c(),mean=c(),sd=c())

png(filename=paste0("R/2_Theta_bar_Y/DtrMngSep/figure/Dtr.png"))

par(mfrow=c(3,3))
for (j in 1:3){
  for (i in 1:3){
    var=paste0("Dtr[",j,",",i,"]")
    ymax=max(density(res_df[,var])$y)
    xmin=min(res_df[,var])
    xmax=max(res_df[,var])
    mean=mean(res_df[,var])
    sd=sd(res_df[,var])
    Dtr_summary=rbind(Dtr_summary,data.frame(var,mean,sd))
    xlab=paste0("mean=",round(mean,2),", sd=",round(sd,2))
    hist(-100,ylim=c(0,ymax),xlim=c(xmin,xmax),ylab="",main=var,xlab=xlab,freq=FALSE)
    lines(density(res_df[,var]),lwd=2)
  }
}

dev.off()

write.csv(Dtr_summary,"R/2_Theta_bar_Y/DtrMngSep/Dtr_summary.csv")

# M

M_summary<-data.frame(coord=c(),mean=c(),sd=c())

png(filename=paste0("R/2_Theta_bar_Y/DtrMngSep/figure/M.png"))

par(mfrow=c(3,3))
for (j in 1:3){
  for (i in 1:3){
    var=paste0("Mnt[",j,",",i,"]")
    ymax=max(density(res_df[,var])$y)
    xmin=min(res_df[,var])
    xmax=max(res_df[,var])
    mean=mean(res_df[,var])
    sd=sd(res_df[,var])
    M_summary=rbind(M_summary,data.frame(var,mean,sd))
    xlab=paste0("mean=",round(mean,2),", sd=",round(sd,2))
    hist(-100,ylim=c(0,ymax),xlim=c(xmin,xmax),ylab="",main=var,xlab=xlab,freq=FALSE)
    lines(density(res_df[,var]),lwd=2)
  }
}

dev.off()

write.csv(M_summary,"R/2_Theta_bar_Y/DtrMngSep/M_summary.csv")


# final plot

library(ggpubr)

p1=c()
p2=c()
p3=c()
for (t in 1:31){
  p1=c(p1,res_df[,paste0("latent_states[",t,",1]")])
  p2=c(p2,res_df[,paste0("latent_states[",t,",2]")])
  p3=c(p3,res_df[,paste0("latent_states[",t,",3]")])
}
predicted_df<-data.frame(t=rep(1:31,each=4000),p1=p1,p2=p2,p3=p3)
predicted_df_2<-data.frame(t=rep(rep(1:31,each=4000),3),p=c(p1,p2,p3),state=rep(1:3,each=31*4000))
predicted_df_2[,3]=factor(predicted_df_2[,3])

ggscatter(predicted_df_2, x = "t", y = "p",
          add = "reg.line",                         # Add regression line
          conf.int = TRUE,                          # Add confidence interval
          color = "state",          # Color by groups "cyl"
          shape = "state"                             # Change point shape by groups "cyl"
)

p <- ggplot(predicted_df_2) + 
  geom_line(aes(y = p, x = t, group = state ))+
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, x =t, fill = "grey70"), alpha = 0.3)+
  scale_colour_manual("", values = "blue")+
  scale_fill_manual("", values = "grey12")


png(filename=paste0("R/2_Theta_bar_Y/DtrMngSep/figure/predicted_states.png"))

ggplot(data =predicted_df_2) + aes(x = t, y = p) +
    geom_point(alpha = 0.1, aes(color=state))


dev.off()

# y hat 

gen<-readRDS("R/2_Theta_bar_Y/DtrMngSep/sample_gen.RDS")\
gen<-as.data.frame(gen)
gen[1,(3188-31*99+1):3188]
gen[,3188]
 
# MSE

MSE_df[,2]=res_df$p21
MSE_df[,3]=res_df$p31
SSE_total = array(0,dim=c(4000,99*31))

to_onehot<-function(x){
  a=rep(0,3)
  a[x]=1
  return (a)
}
onehot_states<-sapply(states,to_onehot)

for (iter in 1:4000){
  for (ind in 1:(99*31)){
    SSE_total[iter,ind]=sum((onehot_states[,ind]-predicted_state[iter,(ind-1)%%31+1,])^2)
  }
}
MSE_df[,5]=apply(SSE_total[,test_ind],1,sum)/5
MSE_df[,4]=apply(SSE_total[,-test_ind],1,sum)/94
png(filename="/Users/choiiiiii/Documents/GitHub/defense-reliability/R/2_Theta_bar_Y/DtrMngSep/figure/train_MSE.png")
par(mfrow=c(1,1))
hist(MSE_df$test_MSE,breaks=100,main="test MSE",xlab="test MSE")
hist(MSE_df$train_MSE,breaks=100,main="train MSE",xlab="train MSE")
#hist(MSE_df$r,breaks=100,main="r",xlab="r")
dev.off()

#rate1, rate2, rate3
D_array=res_df[,15:50]

rate_array<-res_df[,1:12]

rate_vector<-array(0,dim=c(3,4,4000))

par(mfrow=c(1,1))

for(i in 1:3){
  for(era in 1:4){
    for (iter in 1:2000){
      rate_vector[i,era,iter]=rate_array[iter,(i-1)*4+era]
    }
  }
  png(filename=paste0("/Users/choiiiiii/Documents/GitHub/defense-reliability/R/2_Theta_bar_Y/DtrMngSep/figure/lambda_",i,".png"))
  ymax=max(density(rate_vector[i,1,])$y,density(rate_vector[i,2,])$y,density(rate_vector[i,3,])$y,density(rate_vector[i,4,])$y)
  xmin=min(rate_vector[i,,])
  xmax=max(rate_vector[i,,])
  hist(-100,ylim=c(0,ymax),xlim=c(xmin,xmax),main=paste0("Lambda ",i),xlab="",freq=FALSE)
  lines(density(rate_vector[i,1,]),lwd=2)
  lines(density(rate_vector[i,2,]),lwd=2,col="red")
  lines(density(rate_vector[i,3,]),lwd=2,col="blue")
  lines(density(rate_vector[i,4,]),lwd=2,col="green")
  legend("topleft", legend=c("1", "2","3","4"),
         col=c("black","red", "blue","green"), lty=1,
         title="Period", text.font=4, bg='lightblue')
  dev.off()
}


res_df[1,]

par(mfrow=c(3,3))


for(i in 1:3){
  for(j in 1:3){
    if ((j==3) || (i==1 && j==2)){
      hist(c(-0.0000001,rep(0,871)),breaks="fd",main=paste0("Histogram of M[",i,",",j,"]"),xlab=paste0("M[",i,",",j,"]"),xlim=c(-0.5,0.5))
    }
    else if (i==1 && j==1){
      hist(c(1-.0000001,rep(1,871)),breaks="fd",main=paste0("Histogram of M[",i,",",j,"]"),xlab=paste0("M[",i,",",j,"]"),xlim=c(0.5,1.5))
      print(paste0(mean(D_array[1:872,era,i,j]),"is mean of era",era,", i,j:",i,j))
      print(paste0(sd(D_array[1:872,era,i,j]),"is sd of era",era,", i,j:",i,j))
    }
    else if (i==2 && j==1){
      hist(MSE_df$p[1:872],breaks=100,main=paste0("Histogram of M[",i,",",j,"]"),xlab=paste0("M[",i,",",j,"]"),xlim=c(0,1))
      print(paste0(mean(MSE_df$p[1:872]),"is mean of p"))
      print(paste0(sd(MSE_df$p[1:872]),"is sd of p"))
    }
    else if (i==3 && j==1){
      hist(MSE_df$q[1:872],breaks=100,main=paste0("Histogram of M[",i,",",j,"]"),xlab=paste0("M[",i,",",j,"]"),xlim=c(0,1))
      print(paste0(mean(MSE_df$q[1:872]),"is mean of q"))
      print(paste0(sd(MSE_df$q[1:872]),"is sd of q"))

    }
    else if (i==2 && j==2){
      hist(1-MSE_df$p[1:872],breaks=100,main=paste0("Histogram of M[",i,",",j,"]"),xlab=paste0("M[",i,",",j,"]"),xlim=c(0,1))

    }
    else if (i==3 && j==2){
      hist(1-MSE_df$q[1:872],breaks=100,main=paste0("Histogram of M[",i,",",j,"]"),xlab=paste0("M[",i,",",j,"]"),xlim=c(0,1))
    }
  }
}

#final plot




##############################
##############################
# Debug code

rate<-matrix(0,nrow=4,ncol=3)

for (era in 1:4){
  D <- matrix(as.vector(unlist(lapply(1:n_state, function(row){lapply(1:n_state, function(col){res$par[paste0("D[",era,",",row,",",col,"]")]})}))), nrow=3, byrow=T)
  for(j in 1:3){
    rate[era,j] <- res$par[paste0("rate[", era,",", j,"]")]
  }
  h_ <- hist(rate[era,], breaks=50)
  print(h_$counts * h_$breaks)
  print(paste0("Era:",era))
  print(D)
}

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
test_count_df = data.frame(time=rep(1:31,5),observed=c(state_matrix[,test_ship_ind]),ship_ind=factor(rep(test_ship_ind,each=31)),predicted_max=rep(apply(predicted_state,1,which.max),5))

legend_size = 5

ggplot(total_count_df, aes(x = time, y = states)) +
  geom_point(aes(colour="observed", size = observed),alpha=1) +
  geom_point(aes(colour="predicted", size = predicted),alpha=0.4) +
  scale_size(breaks = c(-3:3), range = c(1,30)) +
  guides(
    size=guide_legend(override.aes = list(size = legend_size))
  )

total_df=data.frame(time=rep(1:31,3),states=rep(1:3,each=31),observed=c(observed_count),predicted=c(predicted_state)*99)

ggplot(test_count_df, aes(x = time, y = observed)) +
  geom_line(aes(colour=ship_ind),size=1.2) +
  coord_cartesian(ylim = c(0.5, 3.5)) +
  geom_point(aes(colour=ship_ind),size=3) +
  geom_point(data=total_df, aes(x=time,y=states, size = predicted),color="red",alpha=0.1) +
  scale_size(breaks = c(-3:3), range = c(1,30)) +
  guides(
    size=guide_legend(override.aes = list(size = legend_size))
  )

ggplot(test_count_df, aes(x = time, y = observed)) +
  geom_line(aes(colour=ship_ind),size=1.2) +
  coord_cartesian(ylim = c(0.5, 3.5)) +
  geom_point(aes(colour=ship_ind),size=3) +
  geom_point(aes(x=time,y=predicted_max),size=10,color="red",alpha=0.2)