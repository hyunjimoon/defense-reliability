source("tools/functions.r")
selfCalib <- function(stan_model, prior, pars, data, N, M, cnt, evolve_df, delivDir, is_param = NULL){
  data$theta_loc = mean(prior[[pars]])  #mean(prior) #summarise_draws(workflow$prior_samples)$mean
  data$theta_scale = sd(prior[[pars]]) #summarise_draws(workflow$prior_samples)$sd
  if(is_param){
    generator <- function(){
      function(theta_hp){
        theta <- rnorm(1, theta_hp)
        list(
          generated = rnorm(8, theta, 1),
          parameters = list(
            theta = theta
          )
        )
      }
    }
    theta_hp <- unlist(summarise_draws()[c("mean", "sd")])
    workflow <- SBCWorkflow$new(stan_model, generator())
    workflow$simulate(n_sbc_iterations = N, theta_hp = theta_hp)
    workflow$fit_model(sample_iterations = M, warmup_iterations = M, data)
    next_prior <- workflow$posterior_samples[pars]
  }else{
    # nonparameteric generator but parameteric stan prior
    generator_np <- function(){
      function(theta){
        list(
          generated = rnorm(8, theta, 1),
          parameters = list(
            theta = theta
          )
        )
      }
    }
    workflow <- SBCWorkflow$new(stan_model, generator_np())
    # input samples of multiple parameter
    # custom_prior: https://github.com/hyunjimoon/SBC/blob/927da3f9bc87aca19a34f4dd2061f40eae3176ea/R/util.R#L83
    workflow$simulate(n_sbc_iterations = N, param = is_param, as_draws_df(prior)[[pars]]) #custome_prior draws_of(prior[[par]])
    workflow$fit_model(sample_iterations = M, warmup_iterations = M, data)
    next_prior <- post_summ(workflow, pars, sumtype = "filtering")
  }
  t_prior <- workflow$prior_samples[pars]
  t_post <- workflow$posterior_samples[pars]
  pp_overlay_rvar(t_prior, t_post, pars, cnt)
  if(cnt == 1){
    evolve_df[row(evolve_df)==cnt] <-summarise_draws(t_prior, median, mad)[2:3]
  }else{
    d <- summarise_draws(t_prior, median, mad)[2:3]
    evolve_df <-rbind(evolve_df, as.numeric(d))
  }
  if (iter_stop(t_prior, t_post, bins)){
    csv_store(t_prior, delivDir, cnt)
    csv_store(evolve_df, delivDir, cnt,  type = "evolve")
    intv_plot_save(evolve_df)
    return (prior) # calibrated only for the target
  }
  else{
    csv_store(t_prior, delivDir, cnt)
    cnt = cnt + 1
    return (selfCalib(stan_model, next_prior, pars, data, N, M, cnt, evolve_df, delivDir, is_param = is_param))
  }
}

iter_stop <- function(prior, post, bins = 30){
  post_r_loc <- lapply(post, mean)
  post_r_scale <- lapply(post, sd)
  r_loc <- list()
  r_scale <- list()
  for (par in names(prior)){
    r_loc <- append(r_loc, E(prior[[par]]) / post_r_loc[[par]])
    r_scale <- append(r_scale, sd(prior[[par]]) / post_r_scale[[par]])
  }
  #dist_summary(prior, post, par, bins)$Minkowski_2 < 100
  return (all(r_loc > 0.9 && r_loc < 1.1 && r_scale > 0.9 && r_scale < 1.1 ))
}

# summarize NM posterior samples to N for each parameter
post_summ <- function(workflow, pars, sumtype){
  prior <- subset_draws(workflow$prior_samples, variables = pars)
  post <- subset_draws(workflow$posterior_samples, variables = pars)
  # need to work for multiple pars
  if (sumtype == "filtering"){
    if(length(draws_of(prior[[par]])) < 10){
      return (workflow$posterior_samples[names(prior_rv)])
    }else{
      for (par in pars){
        prior_par <- as_draws_df(prior)[[pars]] #draws_of(prior[[par]])
        post_par <- as_draws_df(post)[[pars]] #draws_of(post[[par]])
        q_prior <- as_draws_df(workflow$prior_samples) # to borrow the frame
        q_prior[[par]]  <- sort(prior_par)
        ar <- as_draws_array(resample_draws(q_prior, tabulate(ecdf(prior_par)(post_par) * N, nbins = length(prior_par))))
        return (as_draws_rvars(aperm(ar, c(2,1,3))))
      }
    }
  }else if(sumtype == "sample"){ # only choose the first
    return (subset_draws(post,variable = names(prior), iteration = 1))
  }
}

initDf <-function(L, summary, pars = NA){
  if(summary == "ms"){
    df <- data.frame(
                     mean = rep(NA,L),
                     sd = rep(NA,L)
    )
  }else if (summary == "q"){
    df <- data.frame(
                     q1 = rep(NA,L),
                     q3 = rep(NA,L)
    )
  }else if (summary == "pars"){
    df <- data.frame(
      median = rep(NA,L),
      mad = rep(NA,L)
      #tau_mean = rep(NA,L),
      #tau_sd = rep(NA,L)
    )
  }
  df
}
