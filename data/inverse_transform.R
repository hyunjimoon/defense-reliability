inverse_yeo_johnson <- function(x, lmbda, mean, std){
  #lambda = 0.21649144
  #mean = 4.96401329
  #std = 2.3079645
  x = x * std + mean
  
  x_inv <- rep(0, length(x))
  positive_x <- which(x >= 0)
  negative_x <- which(x < 0)
  
  if(abs(lmbda) == 0){
    x_inv[positive_x] <- exp(x[positive_x]) - 1
  }
  else{
    x_inv[positive_x] <- (x[positive_x] * lmbda + 1) ** (1 / lmbda) - 1
  }
  
  if(lmbda == 2){
    x_inv[negative_x] = 1 - exp(-x[negative_x])
  }
  else{
    x_inv[negative_x] = 1 - (-(2 - lmbda) * x[negative_x] + 1) ** (1 / (2 - lmbda))
  }
  
  return(as.integer(round(x_inv)))
}

#inverse_yeo_johnson(c(-1.8, -1.31, 1, -1), 0.216,  4.96, 2.3)