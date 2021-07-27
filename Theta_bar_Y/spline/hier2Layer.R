library("dplyr")
library(splines)
data <- readRDS("~/Dropbox/21S_paper/defense-reliability/data_tb/forklift_stanInput.rds")

data_tb <- as_tibble(data)
num_knots <- 15
knot_list <- quantile(unique(data_tb$time), probs=seq(0,1,length.out=num_knots)) 
B <- bs(seq(unique(data_tb$time)), knots=knot_list[-c(1,num_knots)] ,degree=3 , intercept=TRUE )
d <- data_tb %>% select(manufacturer, sn) %>%  unique()
cat <- column_to_rownames(d, var = "sn")
standata <- list("data_tb" =data_tb$manufacturer, "K"= num_knots + 2, "N" = length(data_tb$obs),
                 "T"= max(data_tb$time), "S" = length(unique(data_tb$sn)), "E" = length(unique(data_tb$manufacturer)),
                 "age" = data_tb$time , "obj" = as.numeric(data_tb$sn), "cat" = cat$manufacturer, "B" = B, "Y"=data_tb$obs)
#     int <lower = 1> K; // number of basis functions
#     int <lower = 1> N; // number of total values
#     int <lower = 1> T; // time length of data_tb
#     int <lower = 1> S; // number of obj  e.g ship
#     int <lower = 1> E; // number of unique cat  e.g engine
#     int <lower = 1> age[N];
#     int <lower = 1> obj[N]; // obj type mapping
#     int <lower = 1> cat[S]; //obj_ID to category map
#     matrix[T,K] B; // spline values, 2d array
#     vector [N] Y; // actual values
sm = cmdstanr::cmdstan_model(stan_file = "~/Dropbox/21S_paper/defense-reliability/Theta_bar_Y/spline/models/layer2_nc.stan")
fit_mcmc <- sm$sample(data = standata, iter_warmup = 0, iter_sampling = 1) # test code
#fit_mcmc <- sm$sample(data = standata) # for inference

