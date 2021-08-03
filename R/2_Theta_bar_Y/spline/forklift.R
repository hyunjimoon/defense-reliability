library("dplyr")
library("tibble")
library("loo")
library(splines)
data <- readRDS("~/Dropbox/21S_paper/defense-reliability/data/forklift_1yrBin.rds")

data_tb <- as_tibble(data)
num_knots <- 15
knot_list <- quantile(unique(data_tb$time), probs=seq(0,1,length.out=num_knots))
B <- bs(seq(unique(data_tb$time)), knots=knot_list[-c(1,num_knots)] ,degree=3 , intercept=TRUE )

# 2. per Inspec. freq
sm_IorM = cmdstanr::cmdstan_model(stan_file = "R/2_Theta_bar_Y/spline/models/layer3_nc.stan")
cat_I <- unlist(column_to_rownames(data_tb %>% select(n_inspections, sn) %>%  unique(), var = "sn"))

standata_Inspec <- list("K"= num_knots + 2, "N" = length(data_tb$obs),
                        "T"= max(data_tb$time), "S" = length(unique(data_tb$sn)), "E" = max(data_tb$n_inspections),
                        "age" = data_tb$time , "x" = as.numeric(data_tb$sn),  "B" = B, "Y"=data_tb$obs,
                        "cat" = cat_I)

fit_mcmc_I <- sm_MorI$sample(data = standata_Inspec, iter_warmup = 500, iter_sampling = 500, parallel_chains = 4) # test code

# 3. per manufacturer
cat_M <- unlist(column_to_rownames(data_tb %>% select(manufacturer, sn) %>%  unique(), var = "sn"))

standata_Manuf <- list("K"= num_knots + 2, "N" = length(data_tb$obs),
                       "T"= max(data_tb$time), "S" = length(unique(data_tb$sn)), "E" = max(data_tb$manufacturer),
                       "age" = data_tb$time , "x" = as.numeric(data_tb$sn), "cat" = cat_M, "B" = B, "Y"=data_tb$obs)


fit_mcmc_M <- sm_MorI$sample(data = standata_Manuf, iter_warmup = 500, iter_sampling = 500) # test code


# 4. mixture
cat_IM <- unlist(aperm(array(c(cat_I, cat_M), dim=c(length(unique(data_tb$sn)),2))))
E_IM <- max(max(data_tb$manufacturer), max(data_tb$n_inspections))
standata_MI <- list("K"= num_knots + 2, "N" = length(data_tb$obs),
                    "T"= max(data_tb$time), "S" = length(unique(data_tb$sn)), "C" = 2, "E" = E_IM,
                    "age" = data_tb$time , "x" = as.numeric(data_tb$sn),"cat" = cat_IM,"B" = B, "Y"=data_tb$obs)

sm_MI = cmdstanr::cmdstan_model(stan_file = "R/2_Theta_bar_Y/spline/models/layer3_nc_mixture.stan")
fit_mcmc_MI <- sm_MI$sample(data = standata_MI, iter_warmup = 500, iter_sampling = 500, parallel_chains = 4) # test code
