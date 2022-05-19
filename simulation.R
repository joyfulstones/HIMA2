rm(list=ls(all=TRUE))

simDATMed<-function(n,p,alpha,beta,binaryOutcome=FALSE,seed=123) 
{
  set.seed(seed)
  
  X <- matrix(rnorm(n, mean = 0, sd = 2),n,1) # expoure
  Z <- matrix(rnorm(n*q, mean = 0, sd = 2),n,q) # adjusting covariates
  
  mu <- matrix(0,p,1)
  e <- mvrnorm(n, mu, sigma_e)  # the error terms
  
  M <- matrix(0,n,p)
  M <- X%*%(alpha) + Z%*%t(eta) + e
  MZX <- cbind(M,Z,X)
  
  beta_gamma <- cbind(beta,delta,gamma)
  E_error <- rnorm(n, mean = 0, sd = 1)
  Y <- MZX%*%t(beta_gamma) + E_error 
  if (binaryOutcome) 
    Y <- matrix(rbinom(n, 1, 1/(1 + exp(-Y))), nrow = n)
  return(list(Y = Y, M = M, X = X, Z = Z, n = n, p = p))
}

# Standardize continuous columns of a data frame
std_fun <- function(x) {
  (x - mean(x)) / (sd(x))
}


################ Main ####################
library(survival)
library(MASS)
library(ncvreg)
library(HDMT)
library(hdi)

##############
p=1000 ## p=1000,5000 for each n.
n=300 # n=300, 600 for each p.
q=2

beta <- matrix(0,1,p)
beta[1]  <- 0.20
beta[3] <- 0.25
beta[5]  <- 0.15
beta[7]  <- 0.30
beta[9]  <- 0.35
beta[11]  <- 0.10

##
alpha <- matrix(0,1,p)
alpha[1]  <- 0.20
alpha[3] <- 0.25
alpha[5]  <- 0.15
alpha[7]  <- 0.30
alpha[9]  <- 0.35
alpha[12] <- 0.10


##
delta <- matrix(0.5,1,q)
eta <- matrix(0.3,p,q)
gamma <- matrix(0.5,1,1)

sigma_e <- matrix(0,p,p)
#correlation
rou <- 0  # Need to run for rou=0.25, 0.5, 0.75,0 
for (i in 1:p) {
  for (j in 1:p) {
    sigma_e[i,j]=(rou^(abs(i-j)));
  }
}


rand.num=1001:1500
n.rep=length(rand.num)

ab_est <- matrix(0,n.rep,p)
ab_est_MCP <- matrix(0,n.rep,p)
ab_est_HDMA <- matrix(0,n.rep,p)
P_count_HDMT <- matrix(0,n.rep,p)
P_count_MCP <- matrix(0,n.rep,p)
P_count_HDMA <- matrix(0,n.rep,p)

n_m<-8
ms=c(1,3,5,7,9,11,12,13)
IDS_Step1_HIMA2 <- matrix(0,n.rep,n_m)
IDS_Step1_HIMA <- matrix(0,n.rep,n_m)
IDS_Step1_HDMA <- matrix(0,n.rep,n_m)


for (ii in 1:n.rep)
{
  print("simulation number:")
  print(ii)
  # Generate the data
  simdat<-simDATMed(n, p,alpha,beta,binaryOutcome = FALSE, seed=rand.num[ii]) 
  
  X <-simdat$X
  M<-simdat$M
  Y<-simdat$Y
  Z<-simdat$Z
  
  n <- dim(X)[1]
  p<-dim(M)[2]
  d <-dim(X)[2]
  q<-dim(Z)[2]
  
  #Standardize the exposure 
  #X_exp<- (X - mean(X)) / (sd(X))
  #X<-as.matrix(X_exp)
  
  
  #standardize all the continuous columns and merge with other columns
  # M<-as.data.frame(M)
  # M<-as.data.frame(lapply(M[c(1:p)], std_fun))
  # M<-as.matrix(M)
  
  
  ######################################################################
  ###############step 1
  ######################################################################
  MZX<-cbind(M,Z,X)
  
  ## the new method
  ######  Step 1--the  SIS  step
  
  d_0 <- 2*round(n/log(n))
  beta_SIS <- matrix(0,1,p)
  for (i in 1:p){
    ID_S <- c(i, (p+1):(p+q+1))
    MZX_SIS <- MZX[,ID_S]
    fit <- lsfit(MZX_SIS,Y,intercept = TRUE)
    beta_SIS[i] <- fit$coefficients[2]
  }
  # beta_SIS<-std_fun(beta_SIS)    # standardize the coefficient estimate
  
  ##
  alpha_SIS <- matrix(0,1,p)
  XZ <- cbind(X,Z)
  for (i in 1:p){
    fit_a  <- lsfit(XZ,M[,i],intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    alpha_SIS[i] <- est_a
  }
  # alpha_SIS<-std_fun(alpha_SIS)  # standardize the coefficient estimate
  
  ab_SIS <- alpha_SIS*beta_SIS
  ID_SIS  <- which(-abs(ab_SIS) <= sort(-abs(ab_SIS))[d_0])
  d <- length(ID_SIS)
  
  
  #calculate the frequency of each mediator being selected at step 1 for HIMA2 
  IDS_Step1_HIMA2[ii,]= ms %in% ID_SIS

  
  ######## \Step 2--the following is the code for DLASSO
  P_beta_SIS <- matrix(0,1,d)
  beta_DLASSO_SIS_est <- matrix(0,1,d)
  beta_DLASSO_SIS_SE <- matrix(0,1,d)
  MZX_SIS <- MZX[,c(ID_SIS, (p+1):(p+q+1))]
  
  DLASSO_fit <- lasso.proj(x=MZX_SIS, y=Y, family = "gaussian",Z = NULL)
  beta_DLASSO_SIS_est <- DLASSO_fit$bhat[1:d]
  beta_DLASSO_SIS_SE <- DLASSO_fit$se
  P_beta_SIS <- t(DLASSO_fit$pval[1:d])
  # 
  
  
  ####################  the following is the code for the estimation of alpha
  
  
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
  
  ################### the joint test
  
  PA <- cbind(t(P_alpha_SIS),(t(P_beta_SIS)))
  P_value <- apply(PA,1,max)  # the joint p-values for SIS variable
  
  
  ################  \step 3  the multiple-testing  procedure ####
  
  N0 <- dim(PA)[1]*dim(PA)[2]
  input_pvalues <- PA + matrix(runif(N0,0,10^{-10}),dim(PA)[1],2)
  nullprop <- null_estimation(input_pvalues,lambda=0.5)
  fdrcut  <- fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10, nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
  ID_fdr <- which(fdrcut <= 0.05)
  #print(ID_fdr)
  
  
  if(!is.null(ID_fdr))
  {
    P_count_HDMT[ii, ID_SIS[ID_fdr]] <- 1
  }
  #print(P_count_HDMT[ii,])
  
  ab_est[ii, ID_SIS] <- beta_DLASSO_SIS_est*alpha_SIS_est
  
  #########################   the HIMA method ####################
  
  d <- 2*round(n/log(n))
  ID_SIS <- which(-abs(beta_SIS) <= sort(-abs(beta_SIS))[d])
  

  #calculate the frequency of each mediator being selected at step 1 in HIMA
  IDS_Step1_HIMA[ii,]= ms %in% ID_SIS
  
  
  ####
  MZX_SIS <- MZX[,c(ID_SIS, (p+1):(p+q+1))]
  fit <- ncvreg(MZX_SIS,Y,family = "gaussian",penalty="MCP",alpha=1)
  lam <- fit$lambda[which.min(BIC(fit))]
  
  beta_MCP <- coef(fit, lambda=lam)[2:(d+1)]
  # 
  id_non <- ID_SIS[which(beta_MCP != 0)] # the ID of non-zero by MCP method
  ##
  MZX_MCP <- MZX[,c(id_non, (p+1):(p+q+1))]
  fit <- lsfit(MZX_MCP,Y,intercept = TRUE)
  beta_est_cox <- fit$coefficients[2:(dim(MZX_MCP)[2] + 1)]
  beta_SE_cox <- ls.diag(fit)$std.err[2:(dim(MZX_MCP)[2] + 1)]
  P_val_beta <- 2*(1-pnorm(abs(beta_est_cox)/beta_SE_cox,0,1))
  P_val_beta  <- matrix(P_val_beta, 1, length(id_non))
  ####################
  XZ <- cbind(X,Z)
  alpha_est_MCP <- matrix(0,1,length(id_non))
  alpha_SE_MCP <- matrix(0,1,length(id_non))
  P_alpha_MCP <- matrix(0,1,length(id_non))
  for ( i in 1:length(id_non)){
    fit_a  <- lsfit(XZ,M[,id_non[i]],intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    se_a <- ls.diag(fit_a)$std.err[2]
    sd_1 <- abs(est_a)/se_a
    P_alpha_MCP[i] <- 2*(1-pnorm(sd_1,0,1))  ## the SIS for alpha
    alpha_est_MCP[i] <- est_a
    alpha_SE_MCP[i] <- se_a
  }
  #######  FWER
  # P_JS_MCP <- apply(rbind(P_alpha_MCP, P_val_beta),2,max)*length(id_non)
  # 
  # ID_MCP <- id_non[which(P_JS_MCP <= 0.05)]
  
  ############### B-H for FDR
  d_0 <- length(id_non)
  P_BH <- (1:d_0)*(0.05/d_0)
  P_adj_MCP <- apply(rbind(P_alpha_MCP, P_val_beta),2,max)
  P_sort <- sort(P_adj_MCP)
  SN <- sum(as.numeric(P_sort < P_BH))
  ID_BH <- order(P_adj_MCP)[1:SN]
  ID_MCP <- id_non[ID_BH]
  
  #############
  if(!is.null(ID_MCP))
  {
    P_count_MCP[ii,ID_MCP] <- 1
  }
  
  ab_est_MCP[ii, id_non] <- alpha_est_MCP*beta_MCP[which(beta_MCP != 0)]
  
  #print(ii) 
  
  #########################   the HDMA method ####################
  ####\step 1---This use the screening step of HIMA
  d <- 2*round(n/log(n))
  ID_SIS <- which(-abs(beta_SIS) <= sort(-abs(beta_SIS))[d])
  d <- length(ID_SIS)
  
  #calculate the frequency of each mediator being selected at step 1 in HDMA
  IDS_Step1_HDMA[ii,]= ms %in% ID_SIS

  
  ######## \Step 2--the following is the code for DLASSO
  P_beta_SIS <- matrix(0,1,d)
  beta_DLASSO_SIS_est <- matrix(0,1,d)
  beta_DLASSO_SIS_SE <- matrix(0,1,d)
  MZX_SIS <- MZX[,c(ID_SIS, (p+1):(p+q+1))]
  
  DLASSO_fit <- lasso.proj(x=MZX_SIS, y=Y, family = "gaussian",Z = NULL)
  beta_DLASSO_SIS_est <- DLASSO_fit$bhat[1:d]
  beta_DLASSO_SIS_SE <- DLASSO_fit$se
  P_beta_SIS <- t(DLASSO_fit$pval[1:d])
  # 
  
  ####################  the following is the code for the estimation of alpha
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
  
  ############Step 3: Here we just extract IDs with p_values<=0.05
  P_Val <- apply(rbind(P_alpha_SIS, P_beta_SIS),2,max)
  ID_H <- which(P_Val<=0.05)
  ID_HDMA <- ID_SIS[ID_H]
  
  
  if(!is.null(ID_H))
  {
    P_count_HDMA[ii, ID_SIS[ID_H]] <- 1
  }
  #print(P_count_HDMA[ii,])
  
  ab_est_HDMA[ii, ID_SIS] <- beta_DLASSO_SIS_est*alpha_SIS_est
  
} # the end of n.rep


# save(IDS_Step1_HIMA2,file="IDS_Step1_HIMA2_075.RData")
# save(IDS_Step1_HIMA,file="IDS_Step1_HIMA_075.RData")
# save(IDS_Step1_HDMA,file="IDS_Step1_HDMA_075.RData")

imp_m<-c("M1","M2","M3","M4","M5","M6","M7","M8")
freq_count_step1<-matrix(0,8,4)
freq_count_step1=as.data.frame(freq_count_step1)
colnames(freq_count_step1)=c("M","HIMA2","HIMA","HDMA")
freq_count_step1$M<-imp_m
freq_count_step1[,2]=colSums(IDS_Step1_HIMA2)
freq_count_step1[,3]=colSums(IDS_Step1_HIMA)
freq_count_step1[,4]=colSums(IDS_Step1_HDMA)
freq_count_step1

prob_count_step1=freq_count_step1[,c(2:4)]/n.rep
prob_count_step1

##############################

# save(ab_est, file="ab_est_075.RData")
# save(ab_est_MCP, file="ab_est_MCP_075.RData")
# save(ab_est_HDMA, file="ab_est_HDMA_075.RData")
# save(P_count_HDMT, file="P_count_HDMT_075.RData")
# save(P_count_MCP, file="P_count_MCP_075.RData")
# save(P_count_HDMA, file="P_count_HDMA_075.RData")


#S <- 1:10
S<-c(1,3,5,7,9)
ab_t<-matrix(0,n.rep,p)

for(i in 1:n.rep)
{
  ab_t[i,]<-alpha*beta
}

#####For HDMT##################################
EST_ab_HDMT<-apply(ab_est,2,mean)
Bias_ab_HDMT <- apply(ab_est,2,mean) - alpha*beta
SSE_ab_HDMT <- apply(ab_est,2,sd)
MSE_ab_HDMT <- colMeans((ab_est-ab_t)^2)

#
FDP_our <- matrix(0,1,n.rep)
for (i in 1:n.rep){
  FDP_our[i] <- sum(P_count_HDMT[i,c(2,4,6,8,10:p)])/max(1,sum(P_count_HDMT[i,]))
}
#
FDR_HDMT <- mean(FDP_our)

###
TPR_our_count <- matrix(0,1,n.rep)

for (i in 1:n.rep){
  TPR_our_count[i] <- mean(P_count_HDMT[i,S])
}
Power_HDMT <- mean(TPR_our_count)

TPR_our_count2 <- matrix(0,n.rep,length(S))

for (i in 1:n.rep){
  for(j in 1:length(S))
  {
    TPR_our_count2[i,j] <- mean(P_count_HDMT[i,S[j]])
  }
}
Power_HDMT_NEW <- colMeans(TPR_our_count2)
#####For HIMA##################################
EST_ab_MCP<-apply(ab_est_MCP,2,mean)
Bias_ab_MCP <- apply(ab_est_MCP,2,mean) - alpha*beta
SSE_ab_MCP <- apply(ab_est_MCP,2,sd)
MSE_ab_MCP <- colMeans((ab_est_MCP-ab_t)^2)

FDP_MCP <- matrix(0,1,n.rep)
for (i in 1:n.rep){
  FDP_MCP[i] <- sum(P_count_MCP[i,c(2,4,6,8,10:p)])/max(1,sum(P_count_MCP[i,]))
}
#
FDR_MCP <- mean(FDP_MCP)

TPR_MCP_count <- matrix(0,1,n.rep)
for (i in 1:n.rep){
  TPR_MCP_count[i] <- mean(P_count_MCP[i,S])
}

Power_MCP <- mean(TPR_MCP_count)

TPR_MCP_count2 <- matrix(0,n.rep,length(S))

for (i in 1:n.rep){
  for(j in 1:length(S))
  {
    TPR_MCP_count2[i,j] <- mean(P_count_MCP[i,S[j]])
  }
}
Power_MCP_NEW <- colMeans(TPR_MCP_count2)

#####For HDMA##################################
EST_ab_HDMA<-apply(ab_est_HDMA,2,mean)
Bias_ab_HDMA <- apply(ab_est_HDMA,2,mean) - alpha*beta
SSE_ab_HDMA <- apply(ab_est_HDMA,2,sd)
MSE_ab_HDMA <- colMeans((ab_est_HDMA-ab_t)^2)

#
FDP_HDMA <- matrix(0,1,n.rep)
for (i in 1:n.rep){
  FDP_HDMA[i] <- sum(P_count_HDMA[i,c(2,4,6,8,10:p)])/max(1,sum(P_count_HDMA[i,]))
}
#
FDR_HDMA <- mean(FDP_HDMA)

###
TPR_HDMA_count <- matrix(0,1,n.rep)

for (i in 1:n.rep){
  TPR_HDMA_count[i] <- mean(P_count_HDMA[i,S])
}
Power_HDMA <- mean(TPR_HDMA_count)

TPR_HDMA_count2 <- matrix(0,n.rep,length(S))

for (i in 1:n.rep){
  for(j in 1:length(S))
  {
    TPR_HDMA_count2[i,j] <- mean(P_count_HDMA[i,S[j]])
  }
}
Power_HDMA_NEW <- colMeans(TPR_HDMA_count2)

#############
imp_m<-c("M1","M2","M3","M4","M5","M6","M7","M8")
Final_output_HDMT<-data.frame(CpG=imp_m,"Est_alpha*beta"=EST_ab_HDMT[c(1,3,5,7,9,11,12,13)],"Bias_alpha*beta"=Bias_ab_HDMT[c(1,3,5,7,9,11,12,13)]
                              ,"SSE_alpha*beta"=SSE_ab_HDMT[c(1,3,5,7,9,11,12,13)],"MSE_alpha*beta"=MSE_ab_HDMT[c(1,3,5,7,9,11,12,13)])
Final_output_HDMT


Final_output_MCP<-data.frame(CpG=imp_m,"Est_alpha*beta"=EST_ab_MCP[c(1,3,5,7,9,11,12,13)],"Bias_alpha*beta"=Bias_ab_MCP[c(1,3,5,7,9,11,12,13)]
                             ,"SSE_alpha*beta"=SSE_ab_MCP[c(1,3,5,7,9,11,12,13)],"MSE_alpha*beta"=MSE_ab_MCP[c(1,3,5,7,9,11,12,13)])
Final_output_MCP


Final_output_HDMA<-data.frame(CpG=imp_m,"Est_alpha*beta"=EST_ab_HDMA[c(1,3,5,7,9,11,12,13)],"Bias_alpha*beta"=Bias_ab_HDMA[c(1,3,5,7,9,11,12,13)]
                              ,"SSE_alpha*beta"=SSE_ab_HDMA[c(1,3,5,7,9,11,12,13)],"MSE_alpha*beta"=MSE_ab_HDMA[c(1,3,5,7,9,11,12,13)])
Final_output_HDMA


Power<-c(Power_HDMT,Power_MCP,Power_HDMA)
FDR<-c(FDR_HDMT,FDR_MCP,FDR_HDMA)


Final<-data.frame(CpG=imp_m,"Bias_HIMA2_alpha*beta"=Bias_ab_HDMT[c(1,3,5,7,9,11,12,13)],"Bias_HIMA_alpha*beta"=Bias_ab_MCP[c(1,3,5,7,9,11,12,13)]
                  ,"Bias_HDMA_alpha*beta"=Bias_ab_HDMA[c(1,3,5,7,9,11,12,13)],"SSE_HIMA2_alpha*beta"=SSE_ab_HDMT[c(1,3,5,7,9,11,12,13)],
                  "SSE_HIMA_alpha*beta"=SSE_ab_MCP[c(1,3,5,7,9,11,12,13)],"SSE_HDMA_alpha*beta"=SSE_ab_HDMA[c(1,3,5,7,9,11,12,13)],
                  "MSE_HIMA2_alpha*beta"=MSE_ab_HDMT[c(1,3,5,7,9,11,12,13)],"MSE_HIMA_alpha*beta"=MSE_ab_MCP[c(1,3,5,7,9,11,12,13)],
                  "MSE_HDMA_alpha*beta"=MSE_ab_HDMA[c(1,3,5,7,9,11,12,13)])
Final

Method<-c("HIMA2","HIMA","HDMA")
data.frame(Method,Power,FDR)


###This will results the mediator wise power for each method.
data.frame(cbind(True_mediator=c(1:5),Power_HIMA2=Power_HDMT_NEW,Power_HIMA=Power_MCP_NEW,Power_HDMA=Power_HDMA_NEW))






