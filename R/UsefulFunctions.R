##This file was created on 160722
#It contains a few useful functions, not major functions

logit <- function(x) {   #logistic transform
  if (length(which(x <= 0)) + length(which(x >= 1)) > 0) {
    return(0)
  }
  return( log( x/(1-x) ) )
}

m.which <- function(lookup.vec, x) {  #Input is (lookup vector, things you want to find)
  ind.return <- rep(0, length(x))
  for (i in 1:length(x)) {
    ind.return[i] <- which(lookup.vec == x[i])
  }
  return(ind.return)   #Output is indices of things you are looking up
}

m.grep <- function(patterns, x) {
  ind.return <- rep(0, length(patterns))
  for (i in 1:length(patterns)) {
    ind.return[i] <- grep(pattern=patterns[i], x)[1]
  }
  return(ind.return)
}

Stand.effect.1 <- function(x, y) {  #Standardized effect for the regression y ~ x
  x <- x - mean(x)
  y <- y - mean(y)
  xty <- sum(x * y)
  xtx <- sum(x * x)
  yty <- sum(y * y)
  dof <- length(x) - 2
  sigma.hat <- 1/dof * (yty - xty^2 / xtx)
  return( xty / sqrt(xtx) / sqrt(sigma.hat) )
}

Stand.effect.2 <- function(x, y, dof) {  #Standardized effect for the regression y ~ x + (additional covariates)
  xty <- sum(x * y)
  xtx <- sum(x * x)
  yty <- sum(y * y)
  sigma.hat <- 1/dof * (yty - xty^2 / xtx)
  return( xty / sqrt(xtx) / sqrt(sigma.hat) )
}



##This will compute alpha and beta for the Gamma prior by ML##

log.like.gamma <- function(n, suff, theta) {   #return -ll
  alpha <- exp(theta[1])    #alpha = e^theta_1
  beta <- exp(theta[2])    #beta = e^theta_2
  p <- length(suff)
  return( -(p * alpha * log(beta) - p * lgamma(alpha) + p * lgamma(alpha + n/2) - (alpha + n/2) * sum( log(beta + 1/2 * suff) )) )
}

d.llgamma <- function(n, suff, theta) {   #Derivative of the -ll
  theta_1 <- theta[1]; alpha <- exp(theta_1)
  theta_2 <- theta[2]; beta <- exp(theta_2)
  p <- length(suff)
  dl.d1 <- p * alpha * theta_2 - p * alpha * digamma(alpha) + p * alpha * digamma(alpha + n/2) - alpha * sum(log(beta + 1/2*suff))
  dl.d2 <- p * alpha - (alpha + n/2) * beta * sum(1/(beta + 1/2*suff))
  return(c( -dl.d1, -dl.d2 ))
}
