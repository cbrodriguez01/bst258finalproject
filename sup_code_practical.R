### Supplementary R code using SuperLearner R package for main supplement of Phillips et al. 2023, Practical considerations for specifying a super learner

############################## load R packages #################################

### How to load SuperLearner R package (assuming it's already installed):
library(SuperLearner) # Version: 2.0-28

### At the end of this file, the R session information is provided, which 
### includes the version numbers for all packages used in this file.

##############################################################################
### EXAMPLE 1
### Super learner (SL) specification based on guidelines
##############################################################################

############################### simulate data ##################################
set.seed(59)
n <- 50
X1 <- runif(n = n, min = -5, max = 5)
X2 <- runif(n = n, min = -5, max = 5)
X3 <- runif(n = n, min = -5, max = 5)
X4 <- runif(n = n, min = -5, max = 5)
Y <- 6 + 0.4*X1 - 0.36*(X2^2) + 0.1*(X1>0)*(X3^3) + rnorm(n)
d <- data.frame(Y, X1, X2, X3, X4)
### We will use data "d" in examples 1-4

################################# Part a #######################################
### Use discrete SL (dSL) to predict outcome (Y) from predictors (X) 
library_ex1 <- c("SL.mean", "SL.nnls", "SL.glm", "SL.glm.interaction", 
                 "SL.bayesglm", "SL.glmnet", "SL.ridge", "SL.gam", "SL.earth")
set.seed(924)
sl_fit_ex1 <- SuperLearner(Y = d$Y, X = d[,-1], SL.library = library_ex1, 
                           cvControl = list(V = 20))

### results 
sl_fit_ex1

### ensemble SL (eSL) is the default SL in the SuperLearner R package, with 
### meta-learner "method.NNLS", and this eSL coefficients are column "Coef"
### in the table above that's output when sl_fit_1a is called. (We check its
### CV risk in example 1b.)

### the dSL is the learner with the lowest CV risk. 
### As shown in the table above, that's SL.earth_All, which is learner 
### SL.earth with all covariates (i.e., SL.earth is not coupled with a screener)
### predictions for the dSL can be obtained from library.predict:
pred_dSL_ex1 <- sl_fit_ex1$library.predict[, "SL.earth_All"]


################################## Part b ######################################
### Use dSL to predict Y from X, considering the ensemble SL (eSL) as an 
### additional candidate for the dSL to select by examining its CV risk. 
### NOTE: This requires the specifications in SuperLearner (above) and 
### CV.SuperLearner (below) to be the same.
set.seed(924)
cvsl_fit_ex1 <- CV.SuperLearner(Y = d$Y, X = d[,-1], SL.library = library_ex1, 
                                cvControl = list(V = 20),
                                innerCvControl = list(list(V = 20))) ### CV for essemble learning
summary(cvsl_fit_ex1)
# Risk is based on: Mean Squared Error
# 
# All risk estimates are based on V =  20 
# 
### Column "Ave" is the average CV risk estimate across all folds. The dSL 
### selects the learner with the smallest CV risk. As shown above, that's 
### algorithm "SL.earth_All". 

### Row with "Super Learner" algorithm is the CV performance of the eSL that 
### was fit in example 1a, whose coefficients are provided in example 1a table. 
### If the eSL, "Super Learner", had smallest "Ave" then it would be selected 
### by the dSL and pred_dSL_ex1 could be changed as shown below: 
### pred_dSL_ex1 <- sl_fit_ex1$SL.predict

### Note: 
### - The slight difference in CV risk results in examples 1a and 1b are 
### due to variations in the cross-validation folds. 
### - Algorithm "Discrete SL" in the table above can be ignored. It's not 
### necessary to cross-validate the dSL; the dSL is just a copy of a 
### learner that has already been cross-validated.

##############################################################################
### EXAMPLE 2
### Poor SL library specification
##############################################################################
library_ex2 <- c("SL.mean", "SL.nnls", "SL.glm", "SL.glm.interaction", 
                 "SL.bayesglm", "SL.glmnet", "SL.ridge")
set.seed(924)
sl_fit_ex2 <- SuperLearner(Y = d$Y, X = d[,-1], SL.library = library_ex2,
                           cvControl = list(V = 20))
sl_fit_ex2
set.seed(924)
cvsl_fit_ex2 <- CV.SuperLearner(Y = d$Y, X = d[,-1], SL.library = library_ex2,
                                cvControl = list(V = 20),
                                innerCvControl = list(list(V = 20)))
summary(cvsl_fit_ex2)

##############################################################################
### EXAMPLE 3
### Poor cross-validation specification
##############################################################################
set.seed(924)
sl_fit_ex3 <- SuperLearner(Y = d$Y, X = d[,-1], SL.library = library_ex1,
                           cvControl = list(V = 2))
sl_fit_ex3
set.seed(924)
cvsl_fit_ex3 <- CV.SuperLearner(Y = d$Y, X = d[,-1], SL.library = library_ex1,
                                cvControl = list(V = 2),
                                innerCvControl = list(list(V = 2)))
summary(cvsl_fit_ex3)
##############################################################################
### EXAMPLE 4
### Super learner (SL) specification for rare binary outcome based on guidelines
##############################################################################

############################### simulate data ##################################
### we'll draw on the same simulation as examples 1-3 but this time we will 
### simulate a rare binary outcome. We do this by specifying an arbitrary 
### cutoff, -8.5, that is close to the minimum value of Y.
set.seed(59)
n <- 5000
X1 <- runif(n = n, min = -5, max = 5)
X2 <- runif(n = n, min = -5, max = 5)
X3 <- runif(n = n, min = -5, max = 5)
X4 <- runif(n = n, min = -5, max = 5)
Y <- as.numeric(6 + 0.4*X1 - 0.36*(X2^2) + 0.1*(X1>0)*(X3^3) + rnorm(n) < -8.75)
d_binaryY <- data.frame(Y, X1, X2, X3, X4)

library_ex4 <- c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.bayesglm", 
                 "SL.glmnet", "SL.gam", "SL.earth")
set.seed(924)
sl_fit_ex4 <- SuperLearner(
  Y = d_binaryY[,1], X = d_binaryY[,-1], SL.library = library_ex4, 
  method = "method.NNloglik", family = binomial(),
  cvControl = list(V = 20, stratifyCV = TRUE)
)
sl_fit_ex4

set.seed(924)
cvsl_fit_ex4 <- CV.SuperLearner(
  Y = d_binaryY[,1], X = d_binaryY[,-1], SL.library = library_ex4,
  method = "method.NNloglik", family = binomial(),
  cvControl = list(V = 20, stratifyCV = TRUE),
  innerCvControl = list(list(V = 20, stratifyCV = TRUE))
)
summary(cvsl_fit_ex4)








##############################################################################
### EXAMPLE 5
### Super learner (SL) specification for rare binary outcome with poor 
### cross-validation specification
##############################################################################
set.seed(924)
sl_fit_ex5 <- SuperLearner(
  Y = d_binaryY[,1], X = d_binaryY[,-1], SL.library = library_ex4, 
  method = "method.NNloglik", family = binomial(),
  cvControl = list(V = 5)
)
sl_fit_ex5
set.seed(924)
cvsl_fit_ex5 <- CV.SuperLearner(
  Y = d_binaryY[,1], X = d_binaryY[,-1], SL.library = library_ex4,
  method = "method.NNloglik", family = binomial(),
  cvControl = list(V = 5),
  innerCvControl = list(list(V = 5))
)
summary(cvsl_fit_ex5)
























