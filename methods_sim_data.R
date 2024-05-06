### Running meta-learners and superLearner

library(tidyhte)
library(tidyverse)
library(data.table)
library(MASS)
library(SuperLearner)
## data generation -------------------------------------------------------------


set.seed(100)

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

### No Treatment Effect --------------------------------------------------------

my_data_no_hte <- data.table(L)

#propensity score
pi <- 0.5
# treatment
A <- rbinom(n = n, size = 1, prob = pi)

my_data_no_hte[, A:= A]

#random noise epsilon
e <- rnorm(n)

#mu_a = E[Y|L]-- no treatment effect so these are equal
beta <- matrix(runif(d,min = 1, max=30),ncol = 1, nrow=d)
mu_0 <- L%*%beta
mu_1 <- mu_0
Y1 <- mu_1 + e
Y0 <- mu_0 + e
my_data_no_hte[, Y := A*Y1 + (1-A)*Y0]
my_data_no_hte[, uid := 1:n]
my_data_no_hte[, ps := pi]

my_data_no_hte0 <- my_data_no_hte[A == 0,]
my_data_no_hte1 <- my_data_no_hte[A == 1,]

mean(mu_1-mu_0)
### Unbalanced tx --------------------------------------------------------------

my_data_un_tx <- data.table(L)

#propensity score
pi <- 0.01
# treatment
A <- rbinom(n = n,size = 1,prob = pi)
my_data_un_tx[, A:= A]

#mu_a = E[Y|L]
e <- rnorm(n)
beta<-matrix(runif(d, -5, 5),ncol = 1,nrow=d)
mu_0 <-L%*%beta + 5*(L[,1] > 0.5)
mu_1 <- mu_0 + 8*(L[,2] > 0.1)
true_effect <- mean(mu_1 - mu_0)
Y1<-mu_1 + e
Y0<-mu_0 + e
Y<- A*Y1 + (1-A)*Y0
my_data_un_tx[, Y := A*Y1 + (1-A)*Y0]
my_data_un_tx[, uid := 1:n]
my_data_un_tx[, ps := pi]

my_data_un_tx0 <- my_data_un_tx[A == 0,]
my_data_un_tx1 <- my_data_un_tx[A == 1,]

### Balanced & no confounding --------------------------------------------------

my_data_bal_tx <- data.table(L)
#propensity score
pi <- 0.5
# treatment
A <- rbinom(n = n,size = 1,prob = pi)
my_data_bal_tx[, A:= A]
beta1 <-matrix(runif(d, 1, 30),ncol = 1,nrow=d)
beta0 <-matrix(runif(d, 1, 30),ncol = 1,nrow=d)
mu_1 <- L %*% beta1
mu_0 <- L %*% beta0
mean(mu_1 - mu_0)

Y1 <- mu_1 + e
Y0 <- mu_0 + e
my_data_bal_tx[, Y := A*Y1 + (1-A)*Y0]
my_data_bal_tx[, uid := 1:n]
my_data_bal_tx[, ps := pi]

my_data_bal_tx0 <- my_data_bal_tx[A == 0,]
my_data_bal_tx1 <- my_data_bal_tx[A == 1,]

### SuperLearner & Random forest for CATE estimation  --------------------------

cate_estimation <- function(my_data, my_data0, my_data1, covariates, tx){

hte_superLearner0 <- SuperLearner(Y = my_data0$Y, X = my_data0[, ..covariates], 
                                     family = gaussian(), 
                                     SL.library = c("SL.mean","SL.glmnet", 
                                                    "SL.ranger", 
                                                    "SL.gam"))

hte_superLearner1 <- SuperLearner(Y = my_data1$Y, X = my_data1[, ..covariates], 
                                     family = gaussian(), 
                                     SL.library = c("SL.mean","SL.glmnet", 
                                                    "SL.ranger", 
                                                    "SL.gam"))

all_cov <- c(covariates,tx)
hte_superLearner <- SuperLearner(Y = my_data$Y, X = my_data[,..all_cov], 
                                    family = gaussian(), 
                                    SL.library = c("SL.mean","SL.glmnet", "SL.ranger", 
                                                   "SL.gam"))

hte_m_superLearner <- SuperLearner(Y = my_data$Y, X = my_data[,..covariates], 
                                      family = gaussian(), 
                                      SL.library = c("SL.mean", "SL.glmnet",
                                                     "SL.ranger", 
                                                     "SL.gam"))

## T-Learner 

mu0_hat <- rep(0, n)
mu0_hat[my_data$A==0] <- hte_superLearner0$SL.predict
mu0_hat[my_data$A==1] <- predict(hte_superLearner0, my_data1, 
                                 onlySL = TRUE)$pred




mu1_hat <- rep(0, n)
mu1_hat[my_data$A==1] <- hte_superLearner1$SL.predict
mu1_hat[my_data$A==0] <- predict(hte_superLearner1, my_data0,
                                 onlySL = TRUE)$pred


cate_t_learner <- mu1_hat-mu0_hat


## S-learner

my_dataTMP <- my_data
my_dataTMP$A <- 0 


mu0_hat_s <- rep(0, n) 

mu0_hat_s[my_data$A==0] <- hte_superLearner$SL.predict[my_data$A==0]
mu0_hat_s[my_data$A==1] <- predict(hte_superLearner, my_dataTMP,
                                   onlySL = TRUE)$pred[my_data$A==1]


my_dataTMP$A <- 1 


mu1_hat_s <- rep(0, n) 

mu1_hat_s[my_data$A==1] <- hte_superLearner$SL.predict[my_data$A==1]
mu1_hat_s[my_data$A==0] <- predict(hte_superLearner, my_dataTMP, 
                                   onlySL = TRUE)$pred[my_data$A==0]

cate_s_learner <- mu1_hat_s - mu0_hat_s


## X - Learner

# Compute pseudo-outcome
psi_x_0 <- predict(hte_superLearner1, my_data0, onlySL=T)$pred - my_data0$Y
psi_x_1 <- my_data1$Y - predict(hte_superLearner0, my_data1, onlySL=T)$pred 

# Fit regressions of pseudo outcome against covariates

tau_x_0_fit <- SuperLearner(Y = psi_x_0, X = my_data0[, ..covariates], 
                            family = gaussian(), 
                            SL.library = c("SL.mean", "SL.glmnet", 
                                           "SL.gam"))

tau_x_1_fit <- SuperLearner(Y = psi_x_1, X = my_data1[, ..covariates], 
                            family = gaussian(), 
                            SL.library = c("SL.mean","SL.glmnet", 
                                           "SL.gam"))

# predict TX effects
tau_x_0_hat <- rep(0,n)
tau_x_0_hat[my_data$A==0] <-tau_x_0_fit$SL.predict
tau_x_0_hat[my_data$A==1] <- predict(tau_x_0_fit, my_data1, onlySL=T)$pred 

tau_x_1_hat <- rep(0,n)
tau_x_1_hat[my_data$A==1] <-tau_x_1_fit$SL.predict
tau_x_1_hat[my_data$A==0] <- predict(tau_x_1_fit, my_data0, onlySL=T)$pred 

# Estimate the propensity score: 

pi_fit <- glm(A ~ ., data = my_data[, ..all_cov], family = binomial())

pi_hat <- predict(pi_fit, my_data[, ..all_cov], type="response")


# Ensure positivity
epsilon <- .01 

pi_hat <- ifelse(pi_hat < epsilon,
                  epsilon,
                  ifelse(pi_hat >1-epsilon, 
                         1-epsilon, pi_hat))

cate_x_learner <- pi_hat*tau_x_0_hat +(1-pi_hat)* tau_x_1_hat


## DR Learner

# pseudo-outcome
psi_dr <- mu1_hat - mu0_hat +
  ((1/pi_hat) * (my_data$A * (my_data$Y - mu1_hat)) - 
     (1/(1-pi_hat)) * ((1-my_data$A)*(my_data$Y - mu0_hat)))
# regress pseudo-outcome against covariates 

tau_dr_fit <- SuperLearner(Y = psi_dr, X = my_data[, ..covariates],
                           family = gaussian(), 
                           SL.library = c("SL.mean","SL.glmnet", 
                                          "SL.gam", "SL.ranger"))

# CATE DR-Learner
cate_dr_learner <- c(predict(tau_dr_fit, my_data[, ..covariates], 
                             onlySL = TRUE)$pred)


### R-Learner

# model for conditional mean outcome 

m_hat <- hte_m_superLearner$SL.predict

# pseudo-outcome

residual_tx <- my_data$A - pi_hat
residual_outcome <- my_data$Y - m_hat
psi_r <- residual_outcome/residual_tx

# weights
w = residual_tx^2

# regress pseudo-outcome (randomForest)

tau_dr_fit <- ranger::ranger(y = psi_dr, x =
                        my_data[, ..covariates], keep.inbag =
                        TRUE)

cate_r_learner <- tau_dr_fit$predictions

return(cbind("S-Learner" = cate_s_learner, 
         "T-Learner" = cate_t_learner, 
         "X-Learner" = cate_x_learner,
         "Dr-Learner" = as.numeric(cate_dr_learner),
         "R-Leaner" = cate_r_learner))
}


### NO HTE ---------------------------------------------------------------------
cate_no_hte <- cate_estimation(my_data = my_data_no_hte,
                     my_data0 = my_data_no_hte0,
                     my_data1 = my_data_no_hte1,
                     covariates = colnames(my_data_no_hte[, c(1:5)]),
                     tx = colnames(my_data_no_hte[, c(6)]))

psych::pairs.panels(cate_no_hte, hist.col = "skyblue")

save(cate_no_hte, file = "cate_no_hte.rda")

### Unbalanced tx --------------------------------------------------------------

cate_un_tx <- cate_estimation(my_data = my_data_un_tx,
                              my_data0 = my_data_un_tx0,
                              my_data1 = my_data_un_tx1,
                              covariates = colnames(my_data_un_tx[, c(1:5)]),
                              tx = colnames(my_data_un_tx[, c(6)]))

psych::pairs.panels(cate_un_tx,  hist.col = "plum4")
save(cate_un_tx, file = "cate_un_tx.rda")

### balanced tx & no conf ------------------------------------------------------

cate_bal_tx <- cate_estimation(my_data = my_data_bal_tx,
                              my_data0 = my_data_bal_tx0,
                              my_data1 = my_data_bal_tx1,
                              covariates = colnames(my_data_bal_tx[, c(1:5)]),
                              tx = colnames(my_data_bal_tx[, c(6)]))

psych::pairs.panels(cate_bal_tx,  hist.col = "aquamarine4")
save(cate_bal_tx, file = "cate_bal_tx.rda")
