##### SANDBOX ####

### tidyhte --------------------------------------------------------------------

## Using vignette as provided by authors

# Initial part find nuisance parameter estimation; 
basic_config() %>%
  add_known_propensity_score("ps") %>% # in this case ps is known
  add_outcome_model("SL.gam") %>%
  add_outcome_model("SL.glm.interaction") %>%
  add_outcome_model("SL.glmnet", alpha = c(0, 1)) %>%
  add_outcome_model("SL.glmnet.interaction", alpha = c(0, 1)) %>%
  add_outcome_diagnostic("RROC") %>%
  add_effect_model("SL.gam") %>%
  add_effect_model("SL.glm.interaction") %>%
  add_effect_model("SL.glmnet", alpha = c(0, 1)) %>%
  add_effect_model("SL.glmnet.interaction", alpha = c(0, 1)) %>%
  add_effect_diagnostic("RROC") %>%
  add_moderator("KernelSmooth", cov1,cov2, cov3, cov4, cov5) %>%
  add_vimp(sample_splitting = FALSE) ->
  hte_cfg

# Find quantities of interest based on nuisance estimation and pseudo-outcomes

my_data %>%
  attach_config(hte_cfg) %>%
  make_splits(uid, .num_splits = 3) %>%
  produce_plugin_estimates(
    Y,
    A,
    cov1, cov2, cov3,cov4, cov5,
  ) %>%
  construct_pseudo_outcomes(Y, A) -> prepped_data

prepped_data %>%
  estimate_QoI(cov1, cov2, cov3,cov4, cov5) -> results

results


filter(results, grepl("SATE|PATE", estimand))


filter(results, grepl("SL risk", estimand)) %>%
  mutate(
    level = factor(level, levels = c("Control Response", "Treatment Response", 
                                     "Effect Surface"))
  ) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_pointrange(
    aes(
      x = reorder(term, -estimate),
      y = estimate,
      ymin = estimate - 1.96 * std_error,
      ymax = estimate + 1.96 * std_error)
  ) +
  expand_limits(y = 0) +
  scale_x_discrete("Model name") +
  scale_y_continuous("CV Risk in SuperLearner Ensemble") +
  facet_wrap(~level, scales = "free_x") +
  coord_flip() +
  ggtitle("Submodel Risk Estimates") +
  theme_minimal()


filter(results, grepl("SL coefficient", estimand)) %>%
  mutate(level = factor(level, levels = c("Control Response", "Treatment Response"))) %>%
  ggplot(aes(
    x = reorder(term, estimate),
    y = estimate,
    ymin = estimate - 1.96 * std_error,
    ymax = estimate + 1.96 * std_error
  )) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_pointrange() +
  expand_limits(y = 0) +
  scale_x_discrete("Model name") +
  scale_y_continuous("Coefficient in SuperLearner Ensemble") +
  facet_wrap(~level) +
  coord_flip() +
  ggtitle("SuperLearner Ensemble") +
  theme_minimal()


filter(results, grepl("RROC", estimand)) %>%
  mutate(
    level = factor(level, levels = c("Control Response", 
                                     "Treatment Response", "Effect Surface"))
  ) %>%
  ggplot() +
  geom_line(
    aes(
      x = value,
      y = estimate
    )
  ) +
  geom_point(
    aes(x = value, y = estimate),
    data = filter(results, grepl("RROC", estimand)) %>% group_by(level) %>% slice_head(n = 1)
  ) +
  expand_limits(y = 0) +
  scale_x_continuous("Over-estimation") +
  scale_y_continuous("Under-estimation") +
  facet_wrap(~level, scales = "free_x") +
  coord_flip() +
  ggtitle("Regression ROC Curves") +
  theme_minimal()


ggplot(filter(results, estimand == "VIMP")) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_pointrange(
    aes(
      x = term,
      y = estimate,
      ymin = estimate - 1.96 * std_error,
      ymax = estimate + 1.96 * std_error
    )
  ) +
  expand_limits(y = 0) +
  scale_x_discrete("Covariate") +
  scale_y_continuous("Reduction in R² from full model") +
  coord_flip() +
  ggtitle("Covariate Importance") +
  theme_minimal()


for (cov in c("cov1", "cov2", "cov3", "cov4","cov5")) {
  ggplot(filter(results, estimand == "MCATE", term == cov)) +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
    geom_ribbon(
      aes(
        x = value,
        ymin = estimate - 1.96 * std_error,
        ymax = estimate + 1.96 * std_error
      ),
      alpha = 0.75
    ) +
    geom_line(
      aes(x = value, y = estimate)
    ) +
    expand_limits(y = 0) +
    scale_x_continuous("Covariate level") +
    scale_y_continuous("CATE") +
    ggtitle(paste("Marginal effects across", cov)) +
    theme_minimal() -> gp
  print(gp)
}

## Other essemble combinations are possible. 
## I tried running forest based methods but takes long to run. 

### R-Learner ------------------------------------------------------------------

library(devtools) 
install_github("xnie/rlearner")
library(rlearner)
library(glmnet)

## Here we need to choose a base learner; ie. kernel, lasso, oracle, boost
n.train = 0.7*nrow(my_data)
n.test = 0.2*nrow(my_data)
n.holdout = nrow(my_data) - n.train - n.test

X.all <- as.matrix(my_data[,c(1:5)])
W.all <- my_data$A
Y.all <- my_data$Y

# Train datanrow()# Train data
X = X.all[1:n.train,]
W = W.all[1:n.train]
Y = Y.all[1:n.train]

# Test data
X.test = X.all[n.train + 1:n.test,]
W.test = W.all[n.train + 1:n.test]
Y.test = Y.all[n.train + 1:n.test]

# Holdout data
X.holdout = X.all[n.train + n.test + 1:n.holdout,]
W.holdout = W.all[n.train + n.test + 1:n.holdout]
Y.holdout = Y.all[n.train + n.test + 1:n.holdout]

#
# fit propensity model
#

W.boost = cvboost(X, W, objective = "binary:logistic", nthread = 4)
W.hat.boost = predict(W.boost)


W.lasso = cv.glmnet(X, W, family = "binomial", keep = TRUE)
W.hat.lasso =
  W.lasso$fit.preval[,!is.na(colSums(W.lasso$fit.preval))][, W.lasso$lambda == W.lasso$lambda.min]

# boosting wins
print(round(c(-mean(W * log(W.hat.boost) + (1 - W) * log(1 - W.hat.boost)),
              -mean(W * log(W.hat.lasso) + (1 - W) * log(1 - W.hat.lasso), na.rm = T)), 4))
W.hat = W.hat.boost


#
# fit marginal response model
#

Y.boost = cvboost(X, Y, objective = "reg:squarederror", nthread = 4)
Y.hat.boost = predict(Y.boost)
proc.time()

Y.lasso = cv.glmnet(X, Y, keep = TRUE, family = "gaussian")
Y.hat.lasso =
  Y.lasso$fit.preval[,!is.na(colSums(Y.lasso$fit.preval))][, Y.lasso$lambda == Y.lasso$lambda.min]

# boosting wins
print(round(c(-mean(Y * log(Y.hat.boost) + (1 - Y) * log(1 - Y.hat.boost)),
              -mean(Y * log(Y.hat.lasso) + (1 - Y) * log(1 - Y.hat.lasso))), 4))
Y.hat = Y.hat.boost

#
# fit R-learner given chosen nuisance components
#


tau.boost = rboost(X, W, Y, p_hat = W.hat, m_hat = Y.hat, nthread = 4)
tau.hat.boost = predict(tau.boost)
proc.time()

tau.lasso = rlasso(X, W, Y, p_hat = W.hat, m_hat = Y.hat)
tau.hat.lasso = predict(tau.lasso)

# who wins on CV error?
print(round(c(mean((Y - Y.hat - tau.hat.boost * (W - W.hat))^2),
              mean((Y - Y.hat - tau.hat.lasso * (W - W.hat))^2)), 4))

# lasso wins on holdout
tau.hat.boost.holdout = predict(tau.boost, X.holdout)
tau.hat.lasso.holdout = predict(tau.lasso, X.holdout)


tau.hat.lasso.test = predict(tau.lasso, X.test)

hist(tau.hat.lasso.test)


### t-learner ### Two learners (one model for tx and one model for control)


library(ranger)
# Create separate data frames for the control and the treatment group:
dfs <- my_data
dfs0 <- dfs[dfs$A == 0, ] 
dfs1 <- dfs[dfs$A == 1, ]

# Train a random forest (could be any other learner) for the control group
mu0_fit <- ranger(y = dfs0$Y, x = dfs0 [, c(1:5)], keep.inbag = TRUE)
# Obtain predictions for mu_0
mu0_hat <- rep(0, n)
mu0_hat[dfs$A==0] <- mu0_fit$predictions 
# OOB predictions
mu0_hat[dfs$A==1] <- predict(mu0_fit, dfs1)$predictions

# Train a random forest (could be any other learner) for the treatment group
mu1_fit <- ranger(y = dfs1$Y, x = dfs1 [, c(1:5)], keep.inbag = TRUE)
# Obtain predictions predictions where 
mu1_hat <- rep(0,n)
mu1_hat[dfs$A==1] <- mu1_fit$predictions 
# OOB 
mu1_hat[dfs$A==0] <- predict(mu1_fit, dfs0)$predictions

cate_t <- mu1_hat - mu0_hat

hist(cate_t) ## centered at zero --- what we expect to see


### S-Leaner

## Fit the whole data 

mu_fit <- ranger(y = dfs$Y, x = dfs[, c (1:6)], keep.inbag
                 = TRUE)

dfsTMP <- dfs
dfsTMP$A <- 0
mu0_hat_s <- rep(0, n) 
mu0_hat_s[dfs$A==0] <- mu_fit$predictions[dfs$A==0] 
mu0_hat_s[dfs$A==1] <- predict(mu_fit,dfsTMP)$predictions[dfs$A==1]
# Predict mu_1 by setting A = 1 for all persons, use OOB predictions (see

dfsTMP$A <- 1
mu1_hat_s <- rep(0, n)
mu1_hat_s[dfs$A==1] <- mu_fit$predictions[dfs$A==1] 
mu1_hat_s[dfs$A==0] <- predict(mu_fit,dfsTMP)$predictions[dfs$A==0]
# Compute the CATE as the difference between the predictions by treatment status (see Equation 5):
cate_s <- mu1_hat_s - mu0_hat_s
hist(cate_s)
