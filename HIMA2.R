#' This is the main function for our proposed method for high-dimensional mediation analysis

#' High-dimensional Mediation Analysis Version 2
#' 
#' \code{HIMA2} is used to estimate and test high-dimensional mediation effects.
#' 
#' @param X a vector of exposure. 
#' @param Y a vector of outcome. Can be either continuous or binary (0-1).
#' @param M a \code{data.frame} or \code{matrix} of high-dimensional mediators. 
#' Rows represent samples, columns represent variables.
#' @param Z a \code{data.frame} or \code{matrix} of covariates dataset for testing the association M ~ X and Y ~ M.
#' 
#' @examples
#' ## Generate simulated data
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
#' simdat = simHIMA2(n,p,q,rou,alpha,beta,seed=1234)
#' 
#' ## HIMA2 output
#' hima2.fit <- HIMA2(X=simdat$X, Y=simdat$Y, M=simdat$M, Z=simdat$Z)
#' hima2.fit
#' 
#' @export

HIMA2<-function(X,Y,M,Z)
{
n <- dim(X)[1]  # number of samples
p <- dim(M)[2]  # number of mediators
d <- dim(X)[2]  # number of exposures
q <- dim(Z)[2]  # number of covariates

MZX<-cbind(M,Z,X)

#########################################################################
########################### (Step 1) SIS step ########################### 
#########################################################################
message("Step 1: Sure Independent Screening ...", "  (", Sys.time(), ")")

d_0 <- 2*round(n/log(n)) 
beta_SIS <- matrix(0,1,p) 

# Estimate the regression coefficients beta (mediators --> outcome)
for (i in 1:p){
  ID_S <- c(i, (p+1):(p+q+1))
  MZX_SIS <- MZX[,ID_S]
  fit <- lsfit(MZX_SIS,Y,intercept = TRUE)
  beta_SIS[i] <- fit$coefficients[2]
}

# Estimate the regression coefficients alpha (exposure --> mediators)
alpha_SIS <- matrix(0,1,p)
XZ <- cbind(X,Z)
for (i in 1:p){
  fit_a  <- lsfit(XZ,M[,i],intercept = TRUE)
  est_a <- matrix(coef(fit_a))[2]
  alpha_SIS[i] <- est_a
}

# Select the d_0 number of mediators with top largest effect 
ab_SIS <- alpha_SIS*beta_SIS
ID_SIS  <- which(-abs(ab_SIS) <= sort(-abs(ab_SIS))[d_0])
d <- length(ID_SIS)

#########################################################################
################### (Step 2) De-biased Lasso Estimates ##################
#########################################################################
message("Step 2: De-biased Lasso Estimates ...", "   (", Sys.time(), ")")

P_beta_SIS <- matrix(0,1,d)
beta_DLASSO_SIS_est <- matrix(0,1,d)
beta_DLASSO_SIS_SE <- matrix(0,1,d)
MZX_SIS <- MZX[,c(ID_SIS, (p+1):(p+q+1))]

DLASSO_fit <- lasso.proj(x=MZX_SIS, y=Y, family = "gaussian",Z = NULL)
beta_DLASSO_SIS_est <- DLASSO_fit$bhat[1:d]
beta_DLASSO_SIS_SE <- DLASSO_fit$se
P_beta_SIS <- t(DLASSO_fit$pval[1:d])

################### Estimate alpha ################
alpha_SIS_est <- matrix(0,1,d)
alpha_SIS_SE <- matrix(0,1,d)
P_alpha_SIS <- matrix(0,1,d)

XZ <- cbind(X,Z)
for (i in 1:d){
  fit_a  <- lsfit(XZ,M[,ID_SIS[i]],intercept = TRUE)
  est_a <- matrix(coef(fit_a))[2]
  se_a <- ls.diag(fit_a)$std.err[2]
  sd_1 <- abs(est_a)/se_a
  P_alpha_SIS[i] <- 2*(1-pnorm(sd_1,0,1))  ## the SIS for alpha
  alpha_SIS_est[i] <- est_a
  alpha_SIS_SE[i] <- se_a
}

#########################################################################
################ (step 3) The multiple-testing  procedure ###############
#########################################################################
message("Step 3: Joint significance test ...", "     (", Sys.time(), ")")

PA <- cbind(t(P_alpha_SIS),(t(P_beta_SIS)))
P_value <- apply(PA,1,max)  #The joint p-values for SIS variable

N0 <- dim(PA)[1]*dim(PA)[2]
input_pvalues <- PA + matrix(runif(N0,0,10^{-10}),dim(PA)[1],2)

# Estimate the proportions of the three component nulls
nullprop <- null_estimation(input_pvalues)

# Compute the estimated pointwise FDR for every observed p-max
fdrcut  <- fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10, nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)

ID_fdr <- which(fdrcut <= 0.05)

# Following codes extract the estimates for mediators with fdrcut<=0.05
beta_hat_est <- beta_DLASSO_SIS_est[ID_fdr]
beta_hat_SE  <- beta_DLASSO_SIS_SE[ID_fdr]

alpha_hat_est <-  alpha_SIS_est[ID_fdr]
alpha_hat_SE  <-  alpha_SIS_SE[ID_fdr]

P.value_raw <- P_value[ID_fdr]

# Indirect effect
IDE <- beta_hat_est*alpha_hat_est # mediation(indirect) effect

# Here we name the mediators as M1-Mp and extract the names of significant ones.
M<-(sprintf("M%d", 1:p))[ID_SIS[ID_fdr]]

# create a data frame with output values
output<-data.frame(cbind(M, alpha=alpha_hat_est,alpha_SE=alpha_hat_SE,beta=beta_hat_est,beta_SE=beta_hat_SE,"alpha*beta"=IDE, 
                          p_val=P.value_raw))

message("Done!", "     (", Sys.time(), ")")

return(output)
}