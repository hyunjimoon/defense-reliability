setwd(getwd())

imp<-read.csv("data/y_pred_5var.csv") #imputed data
imp<-c(as.matrix(imp[,-1]))


n=5 #number of levels

#construct state matrix with n levels(31*99)

state_matrix <- function(n){
  state<-cut(imp, breaks=c(quantile(imp,seq(0,1,length.out = n+1))),labels = 1:n, include.lowest=TRUE)
  state<-as.numeric(as.character(state))
  matrix(state,nrow=31)
}

#construct P (n*n)

state_to_P<-function(n,statemat){
  
  cnt<-matrix(0,nrow=n,ncol=n)
  for (i in 1:n){
    for (j in 1:n){
      for (s in 1:99){
        for (a in 2:31){
          if (statemat[a-1,s]==i && statemat[a,s]==j) {cnt[i,j] = cnt[i,j]+1}
        }
      }
    }
  }
  return (apply(cnt,2,function(x) x/apply(cnt,1,sum)))
}

#steady state pi and expectation Q


steady<-function(M,P,n){
  
  r<-rep(1/n,n)
  for (i in 1:80){
    r<-r%*%(M%*%P)
  }
  r
}

Q<-function(pi,n){sum((1:n)*pi)}

### example

n=5
statemat<-state_matrix(n)
P<-state_to_P(n,statemat)

M1<-t(matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0),nrow=n))
M2<-t(matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0),nrow=n))

pi1<-steady(M1,P,n)
Q1<-Q(pi1,n)
pi2<-steady(M2,P,n)
Q2<-Q(pi2,n)

pi1
Q1
pi2
Q2
