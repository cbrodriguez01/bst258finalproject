###WE WILL NONT SH


library(ranger)
library(xgboost)
library(ggplot2)
library(MASS)
library(latex2exp)
library(ipw)
library(cowplot)
library(tidyverse)
library(gbm)
library(tidyverse)
library(DoubleML)
library(mlr3)
library(mlr3learners)



################################################################################
# iTMLE-- i will add some code here. This method is a bit different from the others
################################################################################
set.seed(1994)
# simulate ground truth
n <- 50000
# generate sigma matrix
d <- 4
Sigma <- matrix(0,d,d)
for (i in 1:d){
  for(j in 1:d){
    Sigma[i,j] <- 0.5^(abs(i-j))
  }
}
L <- mvrnorm(n = n, mu = rep(0,d), Sigma = Sigma)
colnames(L) <- c("L1","L2","L3","L4")
#Treatment  & outcome
A <- rbinom(n,1,plogis( L[,1] - 0.5*L[,2] +0.25*L[,3] + 0.1*L[,4])) 
EY <- 21 + A + 27.4*L[,1] + 13.7*L[,2] + 13.7*L[,3] + 13.7*L[,4]
EY1 <- 21 + 1 + 27.4*L[,1] + 13.7*L[,2] + 13.7*L[,3] + 13.7*L[,4]
EY0 <- 21 + 0 + 27.4*L[,1] + 13.7*L[,2] + 13.7*L[,3] + 13.7*L[,4]
Y <- rbinom(n,1, plogis(EY))
Y1 <- rbinom(n,1, plogis(EY1))
Y0 <- rbinom(n,1, plogis(EY0))


source("EstimateAlphat.R")
source("EstimateRR.R")
source("CV-iTMLE.R")
## Preprate Xj -- these have to be separated by levels of the pre-specified groups
### No overlap between subgroups

# generate subgroups
# generate subgroups indices
X1<- L[,1] > quantile(L[,1],0.1)
E.X1 <- sum(Y1[X1])/sum(X1)
X2 <- (L[,2] > quantile(L[,2],0.1)) &(L[,2] < quantile(L[,2],0.9))
E.X2 <- sum(Y1[X2])/sum(X2)
X3 <- (L[,3] + L[,4]) > -2
E.X3 <- sum(Y1[X3])/sum(X3)
X4 <- L[,4] > -1
E.X4 <- sum(Y1[X4])/sum(X4)
X_mat <- cbind(X1, X2, X3, X4)

truePsi <- c(E.X1,E.X2,E.X3,E.X4)
truePsi


###-Risk Ratio iTMLE
# subgroup conditional risk under the treatment arm 
### DRAWBACK-- CANNOT HANDLE MISSING DATA!!
res_1 <- EstimateAlphat(L,A,Y,X_mat,
                        learner="RF", 
                        tr_arm="treatment",
                        CV=FALSE)

cat("The subgroup conditional risks under the treatment arm are:\n ",res_1$subgroup_est)

# subgroup conditional risk under the control arm 
res_0 <- EstimateAlphat(L,A,Y,X_mat,
                        learner="RF",
                        tr_arm="control",
                        CV=FALSE)
cat("The subgroup conditional risks under the control arm are:\n ",res_0$subgroup_est)

### Estimate Risk Ratio-- we use random forest to generate the initial estimate.
# subgroup conditional risk under the treatment arm 
res_RR <- EstimateRR(L,A,Y,X_mat,
                     learner="RF",CV=F)
cat("The subgroup relative risks are:\n ",res_RR$subgroup_RR)
res_RR$Standard_error


###- CV- iTMLE
#the cross-validated version of iTMLE. First, we start with estimating 
#the subgroup conditional risk under the treatment arm.
res_cvitmle<-CViTMLE(L,A,Y,X_mat,  # data
                     learner="RF", # learner for generating initial estimate
                     parameter="RR", # parameter of interest
                     n.folds=3)

res_alpha1$CV_subgroup_est


