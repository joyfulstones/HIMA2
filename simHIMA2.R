#' Simulation Data Generator for High-dimensional Mediation Analysis
#' 
#' \code{simHIMA2} is used to generate simulation data for high-dimensional mediation analysis.
#' 
#' @param n an integer specifying sample size.
#' @param p an integer specifying the dimension of mediators.
#' @param q an integer specifying the dimension of covariates.
#' @param rou a decimal specifying the correlation of exposure.
#' @param alpha a numeric vector specifying the regression coefficients alpha (exposure --> mediators).
#' @param beta a numeric vector specifying the regression coefficients beta (mediators --> outcome).
#' @param binaryOutcome logical. Should the simulated outcome variable be binary?
#' @param seed an integer specifying a seed for random number generation.
#' 
#' @seealso see \code{\link{HIMA2}} to run HIMA2.
#' 
#' @examples
#' n <- 400  # sample size
#' p <- 1000 # the dimension of covariates
#' q <- 2 # the number of adjusted covariates
#' rou <- 0.25 # the correlation of exposure
#' 
#' # the regression coefficients alpha (exposure --> mediators)
#' alpha <- matrix(0,1,p) 
#' 
#' # the regression coefficients beta (mediators --> outcome)
#' beta <- matrix(0,1,p)
#' 
#' # the first five markers are true mediators.
#' alpha[1:5] <- c(0.20,0.25,0.15,0.30,0.35)
#' beta[1:5] <- c(0.20,0.25,0.15,0.30,0.35)
#' 
#' alpha[6] <- 0.1
#' beta[7] <- 0.1
#' 
#' # Generate simulation data
#' simdat = simHIMA2(n,p,q,rou,alpha,beta,seed=1234)
#' 
#' @export

simHIMA2<-function(n,p,q,rou,alpha,beta,binaryOutcome=FALSE,seed=123) 
{
  set.seed(seed)
  
  X <- matrix(rnorm(n, mean = 0, sd = 2),n,1) # expoure
  Z <- matrix(rnorm(n*q, mean = 0, sd = 2),n,q) # covariates
  
  mu <- matrix(0,p,1)
  sigma_e <- matrix(0,p,p)  # correlation matrix
  for (i in 1:p) {
    for (j in 1:p) {
      sigma_e[i,j]=(rou^(abs(i-j)));
    }
  }
  e <- mvrnorm(n, mu, sigma_e)  # the error terms
  
  M <- matrix(0,n,p)
  M <- X%*%(alpha) + Z%*%t(eta) + e  # mediators
  MZX <- cbind(M,Z,X)
  
  beta_gamma <- cbind(beta,delta,gamma)
  E_error <- rnorm(n, mean = 0, sd = 1)
  Y <- MZX%*%t(beta_gamma) + E_error # continuous outcome
  
  if (binaryOutcome) 
    Y <- matrix(rbinom(n, 1, 1/(1 + exp(-Y))), nrow = n)  # binary outcome
  
  return(list(Y = Y, M = M, X = X, Z = Z, n = n, p = p))
}

