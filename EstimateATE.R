# -------------------------------------
# Function for Estimating Average Treatment Effect
# -------------------------------------
library(ranger)
library(ggplot2)
library(MASS)
library(tidyverse)
library(gbm)

source("EstimateAlphat.R")

EstimateATE <- function(X.train,Tr.train,Y.train,A_mat.train,
                        X.val=NULL,Tr.val=NULL,Y.val=NULL,A_mat.val=NULL,
                        learner=NULL,CV=FALSE){
  
  
  if(CV==FALSE){
    X.val <- X.train
    Tr.val <- Tr.train
    Y.val <- Y.train
    A_mat.val <- A_mat.train
  }
  
  est.tr <- EstimateAlphat(X.train,Tr.train,Y.train,A_mat.train,
                           X.val,Tr.val,Y.val,A_mat.val,
                           learner = learner, tr_arm = "treatment",CV=CV)
  est.control <- EstimateAlphat(X.train,Tr.train,Y.train,A_mat.train,
                                X.val,Tr.val,Y.val,A_mat.val,
                                learner = learner, tr_arm = "control",CV=CV)
  
  p1 <- est.tr$raw_est
  mu1 <- est.tr$subgroup_est
  e1 <- est.tr$e1
  
  p0 <- est.control$raw_est
  mu0 <- est.control$subgroup_est
  
  
  
  subgroup_ATE <-  mu1-mu0
  
  # Efficient influence function
  P_A_j <- colSums(A_mat.val)/nrow(A_mat.val)
  IC_mat <- A_mat.val/P_A_j*(Tr.val/e1*(Y.val-p1) +p1 - mean(p1))-
    A_mat.val/P_A_j*((1-Tr.val)/(1-e1)*(Y.val-p0)+p0-mean(p0))
  
  # correlation matrix
  cor.mat <- cor(IC_mat)
  
  Z <- mvrnorm( n = length(Y), mu = rep(0,ncol(IC_mat)), Sigma = cor.mat )
  
  kappa<- quantile(apply(Z,1, max),0.975)
  
  cov.mat <- cov(IC_mat)
  
  CI.length <- NULL
  
  for(k in  1:ncol(cov.mat)){
    CI.length[k] <- kappa * sqrt(cov.mat[k,k])/sqrt(length(Y))
  }
  
  CI.upper <- subgroup_ATE + CI.length
  CI.lower <- subgroup_ATE - CI.length
  
  
  return(list("subgroup_ATE" = subgroup_ATE, "Upper_CI" = CI.upper, "Lower_CI" = CI.lower))
}
