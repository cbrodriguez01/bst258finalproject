

library(data.table)

library(mgcv) #GAM
library(coefplot) #plot coefs--will use for iTMLE
library(BayesTree) #BART
library(rlearner)
library(SuperLearner)
library(causalTree)
library(grf) #Generalized random forests
library(FindIt)
library(htetree)
library(rpart)
library(expss)
library(causal)
library(MASS)
library(devtools)


### Seems like there is an issue with this package from GitHub. 
# devtools::install_github("forestry-labs/causalToolbox") ## this package has the metalearners but i wasnt able to download



################################################################################
# Meta learners
#Looking at different simulation settings from Kunzel et al.
#In this case Y is not binary, but we can try that as well
################################################################################

#Setting 1: No treatment effect

set.seed(1994)
# simulate ground truth
n <- 10000
# generate sigma matrix
d <- 5
Sigma <- matrix(0,d,d)
for (i in 1:d){
  for(j in 1:d){
    Sigma[i,j] <- 0.5^(abs(i-j))
  }
}
L <- mvrnorm(n = n, mu = rep(0,d), Sigma = Sigma)
colnames(L) <- c("cov1", "cov2", "cov3", "cov4","cov5")
my_data <- data.table(L)

#propensity score
pi <- 0.5
# treatment
A <- rbinom(n = n, size = 1, prob = pi)

my_data[, A:= A]

#random noise epsilon
e <- rnorm(n)

#mu_a = E[Y|L]-- no treatment effect so these are equal
beta <- matrix(runif(d,min = 1, max=30),ncol = 1,nrow=d)
mu_0 <- L%*%beta
mu_1 <- mu_0
Y1 <- mu_1 + e
Y0 <- mu_0 + e
my_data[, Y := A*Y1 + (1-A)*Y0]



#Setting 2: Unbalanced treatment assignment
d <- 20
Sigma <- matrix(0,d,d)
for (i in 1:d){
  for(j in 1:d){
    Sigma[i,j] <- 0.5^(abs(i-j))
  }
}
L <- mvrnorm( n = n, mu = rep(0,d), Sigma = Sigma)
#propensity score
pi<-0.01
# treatment
A <- rbinom(n = n,size = 1,prob = pi)
#mu_a = E[Y|L]
beta<-matrix(runif(d, -5, 5),ncol = 1,nrow=d)
mu_0<-L%*%beta + 5*(L[,1] > 0.5)
mu_1<- mu_0 + 8*(L[,2] > 0.1)
Y1<-mu_1 + e
Y0<-mu_0 + e
Y<- A*Y1 + (1-A)*Y0




################################################################################
# iTMLE-- i will add some code here. This method is a bit different from the others
################################################################################



