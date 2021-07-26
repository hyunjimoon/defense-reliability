library(SBC)
true <- read.csv("./data/DMSep_validation_trueparam.csv")
true$X
true$x
modelName = "DMSep"
modDir <- set_get_Dir(modelName)$modDir
delivDir <- set_get_Dir(modelName)$delivDir
file <- set_get_Dir(modelName)$file
simple_normal = cmdstanr::cmdstan_model(file)

p

M

D
for(i in 1:4){
  D[i][1,1] = exp(-(rate[i,1]+ rate[i,2]));
  D[i][2,1] = rate[i,1] * exp(-rate[i,3]) * (1-exp(-(rate[i,1]+ rate[i,2] - rate[i,3]))) / (rate[i,1]+ rate[i,2] - rate[i,3]);
  D[i][1,2] = 0;
  D[i][3,1] = 1 - D[i][1,1] - D[i][2,1];
  D[i][2,2] = exp(-rate[i,3]);
  D[i][3,2] = 1 - D[i][2,2];
  D[i][1,3] = 0;
  D[i][2,3] = 0;
  D[i][3,3] = 1;
}

generator <- function(){  
  function(){
    M[1,1]=1;
    M[1,2]=p;
    M[1,3]=q;
    M[2,1]=0;
    M[2,2]=(1-p);
    M[2,3]= (1-q);
    M[3,1]=0;
    M[3,2]=0;
    M[3,3]= 0; 
    theta <- rnorm(1, 0, 1)
    list(
      generated = rnorm(8, theta, 1),
      parameters = list(
        theta = theta
      )
    )
  }
}
par = "theta"
D <- 8
y <- rep(0, D) # placeholder
data = list("D"= D, "y"= y, "theta_loc"= 0, "theta_scale"= 1)
nChains <- 4
parallel_chains <- min(nChains, detectCores())
N = 1000
M = 10 # 40
workflow <- SBCWorkflow$new(simple_normal, generator())
workflow$simulate(N) 
workflow$fit_model(sample_iterations = M, warmup_iterations = M, data)
prior <- workflow$prior_samples
post <- workflow$posterior_samples
workflow$calculate_rank()
plot_ecdf(workflow, var = "theta")
plot_ecdf_diff(workflow, var="theta")