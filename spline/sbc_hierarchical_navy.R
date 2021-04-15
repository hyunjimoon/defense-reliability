scriptDir <- getwd()
library(cmdstanr)
library(ggplot2)
library(SBC)

hierarchical_model = cmdstan_model("failure_bma/spline/models/layer3_nc_diffsd_parammubar.stan")

data_dir <- file.path(scriptDir, "failure_bma/data")

y_data <- read.csv(file.path(data_dir, "y_count_pwr.csv"))[, "y"]
age_data <- read.csv(file.path(data_dir, "x_age.csv"))[, "age"]
ship_data <- read.csv(file.path(data_dir, "ship_index.csv"))[, "ship"]
engine_data <- read.csv(file.path(data_dir, "engine_type1to5.csv"))[, "engine"]


basis_df <- read.csv(file.path(data_dir, "basis_df.csv"))


data_list <- list(
  K = dim(basis_df)[2],
  N = length(y_data),
  T = 31,
  S = 99,
  E = 5,
  age = age_data,
  engine = engine_data,
  ship = ship_data,
  Y = y_data,
  B = basis_df,
  N_hat = length(y_data),
  age_hat = age_data,
  ship_hat = ship_data
)

#hierarchical_model$sample(data_list)

n_datasets = 2
thin = 3

sbc_obj = SBCModel$new(name="Hierarchical", stan_model = hierarchical_model)

sbc_theta_prior = sbc_obj$sample_theta_tilde_stan(list("mu_a_bar", "mu", "mu_w_bar"), n_datasets, data=data_list)

sbc_sampled_y = sbc_obj$sample_y_tilde(sbc_theta_prior, data=data_list)


theta_post <- sbc_obj$sample_theta_bar_y(sbc_sampled_y, data=data_list, pars=list("mu_a_bar", "mu", "mu_w_bar"), fit_iter = 200)


rank <- calculate_rank(sbc_theta_prior, theta_post, thin = thin)


