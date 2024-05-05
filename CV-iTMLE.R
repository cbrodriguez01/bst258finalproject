
library(ranger)
library(ggplot2)
library(MASS)
library(tidyverse)
library(caret)

source("EstimateAlphat.R")
source("EstimateRR.R")
source("EstimateOR.R")
source("EstimateATE.R")

CViTMLE <- function(X,Tr,Y,A_mat, learner=NULL, parameter=NULL, n.folds=3){
  
  ## create fold list
  fold.list <- createFolds(Y, k = n.folds)
  
  cv.subgroup.est <- matrix(NA,nrow = n.folds, ncol = ncol(A_mat))
  
  for(j in 1:n.folds){
    
    Y.train <- Y[-fold.list[[j]]]
    Tr.train <- Tr[-fold.list[[j]]]
    X.train <- X[-fold.list[[j]],]
    A_mat.train <- A_mat[-fold.list[[j]],]
    
    Y.val <- Y[fold.list[[j]]]
    Tr.val <- Tr[fold.list[[j]]]
    X.val <- X[fold.list[[j]],]
    A_mat.val <- A_mat[fold.list[[j]],]
    
    if(parameter=="Alpha_1"){
      res.alpha1 <- EstimateAlphat(X.train,Tr.train,Y.train,A_mat.train,
                     X.val,Tr.val,Y.val,A_mat.val,
                     learner=learner, tr_arm="treatment",CV=TRUE)
      
      cv.subgroup.est[j,] <- res.alpha1$subgroup_est
      
    }else if(parameter == "Alpha_0"){
      res.alpha0 <- EstimateAlphat(X.train,Tr.train,Y.train,A_mat.train,
                                   X.val,Tr.val,Y.val,A_mat.val,
                                   learner=learner, tr_arm="control",CV=TRUE)
      
      cv.subgroup.est[j,] <- res.alpha0$subgroup_est
      
      }else if(parameter=="RR"){
      res.rr <- EstimateRR(X.train,Tr.train,Y.train,A_mat.train,
                               X.val,Tr.val,Y.val,A_mat.val,
                               learner=learner, CV=TRUE)
      
      cv.subgroup.est[j,] <- res.rr$subgroup_RR
      
    }else if(parameter=="OR"){
      res.or <- EstimateOR(X.train,Tr.train,Y.train,A_mat.train,
                           X.val,Tr.val,Y.val,A_mat.val,
                           learner=learner,CV=TRUE)
      
      cv.subgroup.est[j,] <- res.or$subgroup_OR
      
      
    }else if(parameter=="ATE"){
      res.ate <- EstimateATE(X.train,Tr.train,Y.train,A_mat.train,
                           X.val,Tr.val,Y.val,A_mat.val,
                           learner=learner, CV=TRUE)
      
      cv.subgroup.est[j,] <- res.ate$subgroup_ATE
      
    }else{
      return(NULL)
    }
  }
  
  final.est <- colMeans(cv.subgroup.est,na.rm = TRUE)
   
  return(list("CV_subgroup_est" =  final.est))
}