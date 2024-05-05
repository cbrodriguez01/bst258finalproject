# -------------------------------------
# Function for Estimating alpha_t
# -------------------------------------
library(ranger)
library(ggplot2)
library(MASS)
library(tidyverse)
library(gbm)

EstimateAlphat <- function(X.train,Tr.train,Y.train,A_mat.train,
                           X.val=NULL,Tr.val=NULL,Y.val=NULL,A_mat.val=NULL,
                           learner=NULL, tr_arm="treatment",CV=FALSE){
  if(CV==FALSE){
    X.val <- X.train
    Tr.val <- Tr.train
    Y.val <- Y.train
    A_mat.val <- A_mat.train
  }
  
  ## organize data frames
  train.TX <- cbind(Tr.train,X.train)
  colnames(train.TX)[1] <-"Tr"
  val.TX <- cbind(Tr.val,X.val)
  colnames(val.TX)[1] <-"Tr"
 
  train.YTX <- cbind(Y.train,Tr.train,X.train)
  colnames(train.YTX)[1:2] <-c("Y","Tr")
 
  
  
  
  if(tr_arm=="treatment"){
    print("Estimating subgroup conditional risk under the treatment arm >>>>>>>>>>")
    
    val.mat1 <- cbind(Tr.val,X.val)
    val.mat1[,1] <- rep(1,nrow(val.mat1))
    colnames(val.mat1)[1] <-"Tr"
  
  if(learner =="RF"){
    
    # propensity score model
    psModel <- ranger(x=X.train,y=Tr.train)
    e1 <- predict(psModel,data = as.data.frame(X.val),type="response")
    e1 <- e1$predictions
    
    
    muModel <- ranger(x=train.TX,y=Y.train)
    mu <- predict(muModel,data =val.TX,type="response")
    mu <- mu$predictions
    
    # predict conditional risk under the treatment arm
    mu1 <- predict(muModel,data = val.mat1,type="response")
    mu1 <- mu1$predictions
    
  }else if(learner=="GB"){
    
    # propensity score model
    psModel <- gbm(Tr.train~.,data = data.frame(X.train),distribution = "bernoulli")
    e1 <- predict(psModel,newdata = as.data.frame(X.val),type="response")
    
    # conditional outcome model
    muModel <- gbm(Y~.,data = data.frame(train.YTX),distribution = "bernoulli")
    mu <- predict(muModel,newdata = as.data.frame(val.TX),type="response")
    
    # predict conditional risk under the treatment arm
    mu1 <- predict(muModel,newdata = as.data.frame(val.mat1),type="response")
    
  }else{
    
    # propensity score model
    psModel <- glm(Tr.train~X.train,family = "binomial")
    e1 <- predict(psModel,data = as.data.frame(X.val),type="response")
    
    # conditional outcome model
    muModel <- glm(Y.train~ train.TX,family = "binomial")
    mu <- predict(muModel,data = as.data.frame(val.TX),type="response")
    
    # predict conditional risk under the treatment arm
    mu1 <- predict(muModel,data = as.data.frame(val.mat1),type="response")
  }
  
  # avoid positivity violation
  mu1[mu1==1] <- e1[e1==1] <- mu[mu==1] <- 0.9999
  
  mu1[mu1==0] <- mu[mu==0] <- e1[e1==0] <- 0.0001
  
  
  P_A_j <- colSums(A_mat.val)/nrow(A_mat.val)
  
  if(any(P_A_j==0)){
    return(rep(NA,ncol(A_mat)))
  }
  

  # prep for iterative update-------------------------------
  S_mat <- A_mat.val/P_A_j * Tr.val/e1
  
  mu1.tmle <- mu1
  
  Iter <- 300
  
  neg.log.like.old <- -100
  
  delta0 <- delta <- 0.00001
  
  print("Performing iTMLE update >>>>>>>>>>")
  
  for (j in 1:Iter){
    
    phi <- A_mat.val/P_A_j * Tr.val/e1 * (Y.val-mu1.tmle)
    
    correction<-sqrt(sum((colMeans(phi))^2))
    
    S <- rowSums(S_mat * colMeans(phi))/correction
    
    neg.log.like.new <- sum(Y.val* (qlogis(mu1.tmle) + delta * S) -
                              log(1+exp(qlogis(mu1.tmle) + delta *S)))
    
    if(neg.log.like.new > neg.log.like.old){
      #mu1.tmle <- plogis(logit.mu1 + epsilon1*S)
      mu1.tmle <- plogis(qlogis(mu1.tmle) + delta*S)
      #delta.old <- delta.new
      neg.log.like.old <- neg.log.like.new
    }else{
      delta <- -delta0
    }
    
  }
  

  mustar <- mu1.tmle
  
  subgroup_est <- colSums(A_mat.val*mustar,na.rm = TRUE)/colSums(A_mat.val)
  
 
  
  }else{
    
    print("Estimating subgroup conditional risk under the controlled arm >>>>>>>>>>")
    
    val.mat0 <- data.frame(cbind(Tr.val,X.val))
    val.mat0[,1] <- rep(0,nrow(val.mat0))
    colnames(val.mat0)[1] <- "Tr"
    
    if(learner =="RF"){
      
      # propensity score model
      psModel <- ranger(x=X.train,y=Tr.train)
      e1 <- predict(psModel,data = as.data.frame(X.val),type="response")
      e1 <- e1$predictions
      e0 <- 1-e1
      
      # conditional outcome model
      muModel <- ranger(x=train.TX,y=Y.train)
      mu <- predict(muModel,data = val.TX,type="response")
      mu <- mu$predictions
      
      # predict conditional risk under the control arm
      mu0 <- predict(muModel,data = val.mat0,type="response")
      mu0 <- mu0$predictions
      
    }else if(learner=="GB"){
      
      # propensity score model
      psModel <- gbm(Tr.train~.,data = data.frame(X.train),distribution = "bernoulli")
      e1 <- predict(psModel,newdata = data.frame(X.val),type="response")
      e0 <- 1-e1
      
      # conditional outcome model
      muModel <- gbm(Y~.,data = data.frame(train.YTX),distribution = "bernoulli")
      mu <- predict(muModel,newdata = data.frame(val.TX),type="response")
      
      # predict conditional risk under the control arm
      mu0 <- predict(muModel,newdata = as.data.frame(val.mat0),type="response")
      
    }else{
      
      # propensity score model
      psModel <- glm(Tr.train~X.train,family = "binomial")
      e1 <- predict(psModel,data = as.data.frame(X.val),type="response")
      e0 <- 1-e1
      
      # conditional outcome model
      muModel <- glm(Y.train~ train.TX,family = "binomial")
      mu <- predict(muModel,data = as.data.frame(val.TX),type="response")
      
      # predict conditional risk under the control arm
      mu0 <- predict(muModel,data = as.data.frame(val.mat0),type="response")
    }
    
    # avoid positivity violation
    mu0[mu0==1] <- mu[mu==1] <- e0[e0==1] <- 0.9999
    
    mu0[mu0==0] <- mu[mu==0] <- e0[e0==0] <- 0.0001
    
    P_A_j <- colSums(A_mat.val)/nrow(A_mat.val)
    
    if(any(P_A_j==0)){
      return(rep(NA,ncol(A_mat.val)))
    }
    
    
    # prep for iterative update-------------------------------
    S_mat <- A_mat.val/P_A_j * (1-Tr.val)/e0
    
    mu0.tmle <- mu0
    
    Iter <- 300
    
    neg.log.like.old <- -100
    
    delta0 <- delta <- 0.00001
    
    print("Performing iTMLE update >>>>>>>>>>")
    
    for (j in 1:Iter){
      
      phi <- A_mat.val/P_A_j * (1-Tr.val)/e0 * (Y.val-mu0.tmle)
      
      correction<-sqrt(sum((colMeans(phi))^2))
      
      S <- rowSums(S_mat * colMeans(phi))/correction
      
      neg.log.like.new <- sum(Y.val* (qlogis(mu0.tmle) + delta * S) -
                                log(1+exp(qlogis(mu0.tmle) + delta *S)))
      
      if(neg.log.like.new > neg.log.like.old){
       
        mu0.tmle <- plogis(qlogis(mu0.tmle) + delta*S)

        neg.log.like.old <- neg.log.like.new
      }else{
        delta <- -delta0
      }
      
    }
    
    
    mustar <- mu0.tmle
    
    subgroup_est <- colSums(A_mat.val*mustar,na.rm = TRUE)/colSums(A_mat.val)
    
    
    
  }
  
  
  print("Finishing estimation >>>>>>>>>>")
  
  return(list("raw_est" = mustar,"subgroup_est"=subgroup_est, "e1"=e1))
}
