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
my_data[, uid := 1:n]
my_data[, ps := pi]

my_data0 <- my_data[A == 0,]
my_data1 <- my_data[A == 1,]
### SuperLearner fit -----------------------------------------------------------

no_hte_superLearner0 <- SuperLearner(Y = my_data0$Y, X = my_data0[,c(1:5)], 
                                     family = gaussian(), 
                                     SL.library = c("SL.mean","SL.glmnet", 
                                                    "SL.ranger", 
                                                    "SL.gam"))

table(simplify2array(no_hte_superLearner0$whichDiscreteSL))

no_hte_superLearner1 <- SuperLearner(Y = my_data1$Y, X = my_data1[,c(1:5)], 
                                     family = gaussian(), 
                                     SL.library = c("SL.mean","SL.glmnet", 
                                                    "SL.ranger", 
                                                    "SL.gam"))

table(simplify2array(no_hte_superLearner1$whichDiscreteSL))



no_hte_superLearner <- SuperLearner(Y = my_data$Y, X = my_data[,c(1:6)], 
                                    family = gaussian(), 
                                    SL.library = c("SL.mean","SL.glmnet", "SL.ranger", 
                                                   "SL.gam"))

## T-Learner 

mu0_hat <- rep(0, n)
mu0_hat[my_data$A==0] <- no_hte_superLearner0$SL.predict
mu0_hat[my_data$A==1] <- predict(no_hte_superLearner0, my_data1, 
                                 onlySL = TRUE)$pred




mu1_hat <- rep(0, n)
mu1_hat[my_data$A==1] <- no_hte_superLearner1$SL.predict
mu1_hat[my_data$A==0] <- predict(no_hte_superLearner1, my_data0,
                                 onlySL = TRUE)$pred


cate_t_learner <- mu0_hat-mu1_hat
hist(cate_t_learner)

mean(cate_t_learner)

## S-learner

my_dataTMP <- my_data
my_dataTMP$A <- 0 


mu0_hat_s <- rep(0, n) 

mu0_hat_s[my_data$A==0] <- no_hte_superLearner$SL.predict[my_data$A==0]
mu0_hat_s[my_data$A==1] <- predict(no_hte_superLearner, my_dataTMP,
                                   onlySL = TRUE)$pred[my_data$A==1]


my_dataTMP$A <- 1 


mu1_hat_s <- rep(0, n) 

mu1_hat_s[my_data$A==1] <- no_hte_superLearner$SL.predict[my_data$A==1]
mu1_hat_s[my_data$A==0] <- predict(no_hte_superLearner, my_dataTMP, 
                                   onlySL = TRUE)$pred[my_data$A==0]

cate_s_learner <- mu1_hat_s - mu0_hat_s
hist(cate_s_learner)
mean(cate_s_learner)

## X - Learner

# Compute pseudo-outcome
psi_x_0 <- predict(no_hte_superLearner1, my_data0, onlySL=T)$pred - my_data0$Y
psi_x_1 <- my_data1$Y - predict(no_hte_superLearner0, my_data1, onlySL=T)$pred 

# Fit regressions of pseudo outcome against covariates

tau_x_0_fit <- SuperLearner(Y = psi_x_0, X = my_data0[, c(1:5)], 
                            family = gaussian(), 
                            SL.library = c("SL.mean", "SL.glmnet", 
                                           "SL.gam"))

tau_x_1_fit <- SuperLearner(Y = psi_x_1, X = my_data1[, c(1:5)], 
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

pi_fit <- glm(A ~ ., data = my_data[, c(1:6)], family = binomial())

pi_hat <- predict(pi_fit, my_data[, c(1:6)], type="response")
pi_hat

# Ensure positivity
epsilon <- .01 

pi_hat <- ifelse(pi_hat < epsilon,
                  epsilon,
                  ifelse(pi_hat >1-epsilon, 
                         1-epsilon, pi_hat))

cate_x_learner <- pi_hat*tau_x_0_hat +(1-pi_hat)* tau_x_1_hat
hist(cate_x_learner)
mean(cate_x_learner)

## DR Learner

# pseudo-outcome
psi_dr <- mu1_hat - mu0_hat +
  ((1/pi_hat) * (my_data$A * (my_data$Y - mu1_hat)) - 
     (1/(1-pi_hat)) * ((1-my_data$A)*(my_data$Y - mu0_hat)))
# regress pseudo-outcome against covariates 

tau_dr_fit <- SuperLearner(Y = psi_dr, X = my_data[, c(1:5)],
                           family = gaussian(), 
                           SL.library = c("SL.mean","SL.glmnet", 
                                          "SL.gam"))
tau_dr_fit
# CATE DR-Learner
cate_dr_learner <- predict(tau_dr_fit, my_data[, c(1:5)], onlySL = TRUE)$pred

hist(cate_dr_learner)

mean(cate_dr_learner)

### R-Learner

# model for conditional mean outcome 
no_hte_m_superLearner <- SuperLearner(Y = my_data$Y, X = my_data[,c(1:5)], 
                                      family = gaussian(), 
                                      SL.library = c("SL.mean", "SL.glmnet", "SL.ranger", 
                                                     "SL.gam"))


m_hat <- no_hte_m_superLearner$SL.predict

# pseudo-outcome

residual_tx <- my_data$A - pi_hat
residual_outcome <- my_data$Y - m_hat
psi_r <- residual_outcome/residual_tx

# weights
w = residual_tx^2

# regress pseudo-outcome (randomForest)

tau_dr_fit <- ranger::ranger(y = psi_dr, x =
                        my_data[, c(1:5)], keep.inbag =
                        TRUE)

  # other option is to use a superLearner
  # SuperLearner(Y = psi_r, X = my_data[, c(1:5)],
  #                          obsWeights = w,
  #           family = gaussian(), 
  #           SL.library = c("SL.mean","SL.ranger"))

cate_r_learner <- tau_dr_fit$predictions
hist(cate_r_learner)
mean(cate_r_learner)
