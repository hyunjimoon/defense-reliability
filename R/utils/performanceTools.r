# performance tool: handy functions to analyze results from the computer
# experiment.

TestPerformance <- function(modelName, dim_theta = 10, corr = 1,
                            dataDir, modelDir, delivDir,
                            parametersToPlot = c("lp__", "phi", "rho",
                                                 "theta"),
                            global_pars = c("phi", "rho"),
                            localParms = FALSE,
                            nChains = 4, nPost = 500, nBurn = 500,
                            nThin = 1, plots = TRUE,
                            smart_guess = FALSE, phi = 1.0) {
  # A function to run Stan models and store results, in batches.
  # IMPORTANT: Assumes the model is already compiled.
  # Requires cmdStanTools.r and StanTools.r which should be
  # sourced before hand.
  #
  # modelName: the name of the model
  # . . .
  # smart_guess: if TRUE, solves the relevant algebraic equation
  #              of the Laplace problem beforehand, and passes the 
  #              solution to the model as the default initial guess.
  #              If not, the initial guess is a vector of mean values
  #              for theta.
  # phi: specifies the value we used when computing the initial
  #      guess. Only relevant is smart_guess = TRUE.

  dataFile <- file.path(dataDir, paste0("data_", dim_theta, ".R"))
  cat("Marker A")
  cat(dataFile)
  data <- read_rdump(dataFile)
  cat("Marker B")

  M = data$M
  if (localParms == TRUE) { init <- function() {
      list(phi = runif(1, 0.8, 1.2),
           rho = runif(1, 0, 1),
           theta = rnorm(M, sd = 2))
    }
  }
  if (localParms == FALSE) { init <- function() {
    list(phi = runif(1, 0.8, 1.2),
         rho = runif(1, 0, 1))
    }
  }

  nIter <- nPost * nThin
  nBurnin <- nBurn * nThin

  RNGkind("L'Ecuyer-CMRG")
  mc.reset.stream()

  dir.create(file.path(modelDir, modelName))
  dir.create(file.path(modelDir, modelName, "temp"))

  chains <- 1:nChains
  mclapply(chains,
           function(chain, model, data, iter, warmup, thin, init) {
             tempDir <- file.path(tempDir, chain)
             dir.create(tempDir)
             inits <- init()
             with(inits, stan_rdump(ls(inits), file = file.path(tempDir,
                                                                "init.R")))
             runModel(model = model, data = data,
                      iter = iter, warmup = warmup, thin = thin,
                      init = file.path(tempDir, "init.R"), 
                      seed = sample(1:999999, 1),
                      chain = chain, refresh = 100,
                      adapt_delta = 0.95, stepsize = 0.01)
           },
           model = file.path(modelDir, modelName),
           data = dataFile,
           init = init,
           iter = nIter, warmup = nBurnin, thin = nThin,
           mc.cores = min(nChains, detectCores()))

  fit <- read_stan_csv(file.path(modelDir, modelName, paste0(modelName, chains, ".csv")))
  dir.create(outDir)
  save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

  ###############################################################################
  ## posterior distributions of parameters
  dir.create(dirname(dirname(delivDir)))
  dir.create(dirname(delivDir))
  dir.create(delivDir)

  # gaussianParms <- c("theta[1]", "theta[2]")

  ## open graphics device
  if (plots) {
    pdf(file = file.path(delivDir, paste(modelName,"Plots%03d.pdf", sep = "")),
        width = 6, height = 6, onefile = F)

    mcmcHistory(fit, parametersToPlot)
    mcmcDensity(fit, parametersToPlot, byChain = TRUE)
    mcmcDensity(fit, parametersToPlot)
    pairs(fit, pars = parametersToPlot)
    # pairs(fit, pars = gaussianParms)

    dev.off()
  }

  ptable <- parameterTable(fit, parametersToPlot)
  write.csv(ptable, file = file.path(delivDir, paste(modelName, 
                                                     "ParameterTable.csv", sep = "")))

  ## Read out some performance summaries
  # n_eff is measured by summing the neff over the two global parameters.
  elapsed_time <- colSums(get_elapsed_time(fit))[[2]]
  summary <- summary(fit, pars = global_pars)$summary
  # n_eff <- sum(summary[, 9])
  n_eff <- summary[1, 9]
  n <- nChains * nPost
  ratio <- n_eff / n
  rho <- corr
  efficiency <- n_eff / elapsed_time

  performance <- data.frame(elapsed_time, ratio, rho, efficiency)
  write.csv(performance, file = file.path(delivDir, "performance.csv"))
}
  

read_performance <- function (model, dim_theta, corr, type_number = 3) {
  delivDir <- file.path("deliv", model, data_types[type_number])
  data_full <- data.frame()
  for (k in 1:length(dim_theta)) {
    data <- read.csv(file.path(delivDir, dim_theta[k], corr[1],
                               "performance.csv"))[, 2:4]
    for (j in 2:length(corr)) {
      data_new <- read.csv(file.path(delivDir, dim_theta[k], corr[j],
                                     "performance.csv"))[, 2:4]
      data <- rbind(data, data_new)
    }
    data$n_eff <- 4000 * data$ratio
    data$efficiency <- data$n_eff / data$elapsed_time
    data$dim_theta <- dim_theta[k]
    
    data_full <- rbind(data_full, data)
  }
  data_full
}
  