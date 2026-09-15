
# Cross-Validation for ToM & Pain Network Activity

# Note: This script generates a synthetic dummy dataset matching the 
# expected structure of clinical & neuroimaging data.


# Load packages
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("boot", quietly = TRUE)) install.packages("boot")
if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")

library(dplyr)
library(boot)
library(tibble)

# Generate dummy dataset

set.seed(42) # For exact reproducibility
n_sim <- 150 # Exemplary sample size

dummy_data <- tibble(
  ID = sprintf("sub-%03d", 1:n_sim),
  
  # Target variables (e.g., mean BOLD contrast estimates/beta weights)
  ToM_Network_Activity  = rnorm(n_sim, mean = 0.35, sd = 0.25),
  Pain_Network_Activity = rnorm(n_sim, mean = 0.28, sd = 0.20),
  
  # ASR DSM-oriented scales (eg t-scores or raw scores)
  ASR_DSM_Depressive_Problems             = round(rnorm(n_sim, mean = 55, sd = 10)),
  ASR_DSM_Anxiety_Problems                = round(rnorm(n_sim, mean = 54, sd = 9)),
  ASR_DSM_Somatic_Problems                = round(rnorm(n_sim, mean = 52, sd = 8)),
  ASR_DSM_Avoidant_Personality_Problems   = round(rnorm(n_sim, mean = 53, sd = 10)),
  ASR_DSM_Inattention                     = round(rnorm(n_sim, mean = 56, sd = 11)),
  ASR_DSM_Hyperactivity                   = round(rnorm(n_sim, mean = 51, sd = 8)),
  ASR_DSM_Antisocial_Personality_Problems = round(rnorm(n_sim, mean = 50, sd = 7)),
  
  # ASR broad-band overscales
  ASR_Internalizing    = round(rnorm(n_sim, mean = 58, sd = 12)),
  ASR_Externalizing    = round(rnorm(n_sim, mean = 52, sd = 11)),
  ASR_Thought_Problems = round(rnorm(n_sim, mean = 53, sd = 9))
)

# When using empirical data, load it here:
# my_data <- read.csv("path/to/real_data.csv")
my_data <- dummy_data


# function for LOOCV and out-of-sample R²

calc_loocv_metrics <- function(model, data, y_var_name) {
  n <- nrow(data)
  cv_res <- cv.glm(data = data, glmfit = model, K = n)
  
  # delta[1] = raw LOOCV MSE
  mse_raw  <- cv_res$delta[1]
  var_tot  <- var(data[[y_var_name]], na.rm = TRUE)
  
  # Out-of-sample R² (can be negative under overfitting/pure noise)
  r2_loocv <- 1 - (mse_raw / var_tot)
  
  return(list(
    mse_raw  = mse_raw,
    r2_loocv = r2_loocv
  ))
}


# Data prep 

analysis_data <- my_data %>%
  select(
    ToM_Network_Activity,
    Pain_Network_Activity,
    ASR_DSM_Depressive_Problems,
    ASR_DSM_Anxiety_Problems,
    ASR_DSM_Somatic_Problems,
    ASR_DSM_Avoidant_Personality_Problems,
    ASR_DSM_Inattention,
    ASR_DSM_Hyperactivity,
    ASR_DSM_Antisocial_Personality_Problems,
    ASR_Internalizing,
    ASR_Externalizing,
    ASR_Thought_Problems
  ) %>%
  na.omit()


# GLM modeling
# Note: ASR_DSM_ADH_Problems is excluded to avoid multicollinearity/rank deficiency

# ToM
glm_tom_dsm <- glm(
  ToM_Network_Activity ~ ASR_DSM_Depressive_Problems + ASR_DSM_Anxiety_Problems +
    ASR_DSM_Somatic_Problems + ASR_DSM_Avoidant_Personality_Problems +
    ASR_DSM_Inattention + ASR_DSM_Hyperactivity +
    ASR_DSM_Antisocial_Personality_Problems,
  data = analysis_data, family = gaussian()
)

glm_tom_over <- glm(
  ToM_Network_Activity ~ ASR_Internalizing + ASR_Externalizing + ASR_Thought_Problems,
  data = analysis_data, family = gaussian()
)

# Pain
glm_pain_dsm <- glm(
  Pain_Network_Activity ~ ASR_DSM_Depressive_Problems + ASR_DSM_Anxiety_Problems +
    ASR_DSM_Somatic_Problems + ASR_DSM_Avoidant_Personality_Problems +
    ASR_DSM_Inattention + ASR_DSM_Hyperactivity +
    ASR_DSM_Antisocial_Personality_Problems,
  data = analysis_data, family = gaussian()
)

glm_pain_over <- glm(
  Pain_Network_Activity ~ ASR_Internalizing + ASR_Externalizing + ASR_Thought_Problems,
  data = analysis_data, family = gaussian()
)


# LOOCV
res_tom_dsm   <- calc_loocv_metrics(glm_tom_dsm,   analysis_data, "ToM_Network_Activity")
res_tom_over  <- calc_loocv_metrics(glm_tom_over,  analysis_data, "ToM_Network_Activity")
res_pain_dsm  <- calc_loocv_metrics(glm_pain_dsm,  analysis_data, "Pain_Network_Activity")
res_pain_over <- calc_loocv_metrics(glm_pain_over, analysis_data, "Pain_Network_Activity")


# Results summary table

loocv_summary <- tibble(
  Target_Network = c("ToM", "ToM", "Pain", "Pain"),
  Predictors     = c("ASR DSM-Scales", "ASR Overscales", "ASR DSM-Scales", "ASR Overscales"),
  N              = nrow(analysis_data),
  LOOCV_MSE      = c(res_tom_dsm$mse_raw, res_tom_over$mse_raw, res_pain_dsm$mse_raw, res_pain_over$mse_raw),
  LOOCV_R2       = c(res_tom_dsm$r2_loocv, res_tom_over$r2_loocv, res_pain_dsm$r2_loocv, res_pain_over$r2_loocv)
)

print(loocv_summary)

