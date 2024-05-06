# bst258finalproject
BIOSTAT 258 final project: A Survey of  Heterogeneous Treatment Effects by  Monica Robles, Carmen Rodriguez and George Sawyer
- **methods_sim_data.R** contains simulation results from using meta-learners to estimate the CATE for 3 scenarios.
  
- **iTMLE.R** contains a small simulation to show how to implement the method proposed by Wei et al.. 2023 to estimate HTE across multiple subgroups using TMLE.
       - This file calls on EstimateATE.R, EstimateRR.R, EstimateOR.R, and CV-iTMLE.R, which are from the author's repository. See their paper.

- The remaining files  are some examples we found online that provide a guide for using the SuperLearner package and the and tidyHTE.

