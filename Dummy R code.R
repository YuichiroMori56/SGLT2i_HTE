# --- Annotated R Code for Reproducibility ---
# This script provides a simplified, annotated example of the methodology used in the paper:
# "Heterogeneous Effect of Sodium-Glucose Cotransporter 2 Inhibitors on Cardiovascular Events
# among People with Type 2 Diabetes: Machine-Learning Causal Forest Method and Target
# Trial Emulation"

# --- 0. Load Required Libraries and Create Dummy Data ---

# Load necessary packages
# install.packages(c("tidyverse", "survival", "MatchIt", "grf", "ggplot2", "lmtest", "sandwich"))
library(tidyverse)
library(survival)
library(MatchIt)
library(grf)
library(ggplot2)
library(lmtest)
library(sandwich)

# Generate a comprehensive dummy dataset including covariates and outcomes
set.seed(1234)
n_patients <- 10000
all_patients <- data.frame(
  id = 1:n_patients,
  age = rnorm(n_patients, mean = 55, sd = 8),
  sex = factor(rbinom(n_patients, 1, 0.6), labels = c("Female", "Male")),
  bmi = rnorm(n_patients, mean = 27, sd = 4),
  sbp = rnorm(n_patients, mean = 125, sd = 17),
  hba1c = rnorm(n_patients, mean = 7.6, sd = 1.5),
  egfr = rnorm(n_patients, mean = 80, sd = 18),
  hist_cvd = rbinom(n_patients, 1, 0.15),
  # Simulate cohort entry and SGLT2i initiation dates in months
  cohort_entry_month = sample(-12:12, n_patients, replace = TRUE),
  sglt2_initiation_month = round(runif(n_patients, 0, 12))
) 
# Some patients never initiate SGLT2i (these will serve as potential controls)
all_patients$sglt2_initiation_month[sample(1:n_patients, size = n_patients * 0.7)] <- NA

# Simulate patient outcomes with built-in Heterogeneous Treatment Effects (HTE)
all_patients <- all_patients %>%
  mutate(
    ever_treated = ifelse(is.na(sglt2_initiation_month), 0, 1),
    
    # Define a treatment effect that varies with patient characteristics (HTE)
    # Benefit (hazard reduction) increases with higher SBP, BMI, and HbA1c
    # Benefit slightly decreases with older age
    treatment_effect_log_hr = ever_treated * (
      -0.30                              # Average treatment effect (log Hazard Ratio)
      - 0.090 * (sbp - 112)              # Effect modifier: SBP
      - 0.065 * (bmi - 27)               # Effect modifier: BMI
      - 0.080 * (hba1c - 7.6)            # Effect modifier: HbA1c
      + 0.015 * (age - 55)               # Effect modifier: Age
      + 0.002 * (age - 55) * (bmi - 27)  # Complex interactions
      + 0.003 * (age - 55) * (sbp - 112)
    ),
    
    # Define baseline risk (prognostic effect)
    prognostic_log_hr = 0.03 * (age - 55) + 0.01,
    
    # Calculate final hazard incorporating baseline risk and HTE
    hazard = 0.005 * exp(prognostic_log_hr + treatment_effect_log_hr),
    
    time_to_event = rexp(n(), rate = hazard),
    # Follow-up extends up to 84 months
    time_to_censor = runif(n(), 0, 84),
    survival_time = pmin(time_to_event, time_to_censor),
    event = as.numeric(time_to_event < time_to_censor)
  )


# --- 1. Target Trial Emulation via Sequential Propensity Score Matching ---

# This loop emulates the target trial by creating monthly exposure datasets
# where new users of SGLT2i are matched to eligible non-users chronologically.

# NOTE (1): For demonstration, this loop is simplified to perform matching only
# during the first 13 months (0 to 12). The original study matched for the entire period.

# NOTE (2): The original study used time-varying covariates updated at each
# matching month. This sample code uses fixed baseline covariates for simplicity.

monthly_matched_cohorts <- list()
matched_patient_ids <- c()

# Loop through a simplified calendar period (first 13 months)
for (index_month in 0:12) {
  
  # Identify SGLT2i new users for the current month
  sglt2_new_users <- all_patients %>%
    filter(
      !id %in% matched_patient_ids,
      cohort_entry_month <= index_month,
      sglt2_initiation_month == index_month
    )
  
  # Identify eligible controls for the current month
  eligible_controls <- all_patients %>%
    filter(
      !id %in% matched_patient_ids,
      cohort_entry_month <= index_month,
      is.na(sglt2_initiation_month) | sglt2_initiation_month > index_month
    )
  
  if (nrow(sglt2_new_users) == 0 || nrow(eligible_controls) == 0) {
    next
  }
  
  # Create the dataset for this month's matching
  monthly_matching_pool <- bind_rows(
    sglt2_new_users %>% mutate(treatment = 1),
    eligible_controls %>% mutate(treatment = 0)
  )
  
  # Perform propensity score matching
  match_obj_month <- matchit(
    treatment ~ age + sex + bmi + sbp + hba1c + egfr + hist_cvd,
    data = monthly_matching_pool,
    method = "nearest", ratio = 1
  )
  
  # Store the matched data and update the list of matched IDs
    monthly_matched_data <- match.data(match_obj_month) %>%
      mutate(match_month = index_month)
    monthly_matched_cohorts[[as.character(index_month)]] <- monthly_matched_data
    matched_patient_ids <- c(matched_patient_ids, monthly_matched_data$id)
}

# Combine all monthly matched cohorts into the final analytic dataset
emulated_trial_data <- bind_rows(monthly_matched_cohorts)


# --- Post-Matching Analyses (Steps 2-8) ---
# The following analyses are performed on the `emulated_trial_data`.

# --- 2. Outcome Modelling via Adjusted Cox Regression (Average Treatment Effect) ---

# Estimate the Average Treatment Effect (ATE) using a Cox proportional hazards model.
# After matching, the model is further adjusted for covariates to address residual confounding.
cox_model <- coxph(
  Surv(survival_time, event) ~ treatment + age + sex + bmi + sbp + hba1c + egfr + hist_cvd + cluster(subclass),
  data = emulated_trial_data
)
summary(cox_model)

# --- 3. Causal Forest Training (for Heterogeneous Treatment Effects) ---

# Define the 3-year binary outcome and censoring status.
emulated_trial_data <- emulated_trial_data %>%
  mutate(
    outcome_3yr = ifelse(survival_time < 36 & event == 1, 1, 0),
    censored_3yr = ifelse(survival_time < 36 & event == 0, 1, 0)
  )

# For the main analysis, we use only subjects with complete 3-year outcome data.
analysis_data_complete <- emulated_trial_data %>% filter(censored_3yr == 0)

# Define components for the causal forest model using the complete-case data
X <- analysis_data_complete %>% select(age, sex, bmi, sbp, hba1c, egfr, hist_cvd)
Y <- analysis_data_complete$outcome_3yr
W <- analysis_data_complete$treatment

# Define folds for cross-fitting
num_folds <- 10
set.seed(1234)
folds <- sample(1:num_folds, size = nrow(X), replace = TRUE)

# Train the causal forest. `tune.parameters = "all"` automatically tunes key hyperparameters.
causal_forest_model <- causal_forest(
  X = data.matrix(X), Y = Y, W = W,
  num.trees = 5000,
  honesty = TRUE,
  tune.parameters = "all", # Enable automatic hyperparameter tuning
  clusters = folds,
  seed = 1234
)


# --- 4. Estimation of Individualised Treatment Effects (ITEs) ---

# Predict ITEs by retrieving the out-of-fold predictions from the trained model.
# This is a result of the cross-fitting (`clusters` argument) and is crucial for honest estimation.
tau.hat <- predict(causal_forest_model)$predictions
analysis_data_complete$ite <- tau.hat


# --- 5. Calibration and Validation of ITE Estimates ---

calibration_results <- test_calibration(causal_forest_model)
print(calibration_results)

# --- 6. Variable Importance Ranking ---

var_importance <- variable_importance(causal_forest_model)
var_importance_df <- data.frame(
  Variable = colnames(X),
  Importance = var_importance
) %>% arrange(desc(Importance))
ggplot(var_importance_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity") + coord_flip() + theme_minimal()

# --- 7. Stratification of ITEs and CATEs using a Doubly-Robust Estimator ---

num.ranking = 5

ranking <- rep(NA, nrow(analysis_data_complete))
for(fold in seq(num_folds)) {
  tau.hat.quantiles <- quantile(tau.hat[folds == fold], probs = seq(0,1, by = 1/num.ranking))
  ranking[folds == fold] <- cut(tau.hat[folds == fold], tau.hat.quantiles, include.lowest = TRUE, labels = seq(num.ranking))
}


# Stratify the population into quintiles based on their predicted ITE.
analysis_data_complete$ranking <- as.factor(ranking)

# Manually compute the doubly-robust AIPW (Augmented Inverse Propensity Weighting) scores.
# This provides a robust estimate of the treatment effect for each subgroup.

# Extract nuisance predictions from the model object
e_hat <- causal_forest_model$W.hat # Propensity score: P(W=1|X)
m_hat <- causal_forest_model$Y.hat # Conditional mean outcome: E(Y|X)

# Estimate potential outcomes mu_0(X) and mu_1(X)
mu_hat_0 <- m_hat - e_hat * analysis_data_complete$ite
mu_hat_1 <- m_hat + (1 - e_hat) * analysis_data_complete$ite

# Compute the AIPW scores
aipw_scores <- analysis_data_complete$ite +
  (analysis_data_complete$treatment / e_hat) * (analysis_data_complete$outcome_3yr - mu_hat_1) -
  ((1 - analysis_data_complete$treatment) / (1 - e_hat)) * (analysis_data_complete$outcome_3yr - mu_hat_0)

ols <- lm(aipw_scores ~ 0 + factor(analysis_data_complete$ranking))
forest.ate <- data.frame("aipw", paste0("Q", seq(num.ranking)), coeftest(ols, vcov = vcovHC(ols, "HC2"))[,1:2])
colnames(forest.ate) <- c("method", "ranking", "estimate", "std.err")
rownames(forest.ate) <- NULL

res <- rbind(forest.ate)
res

ggplot(res %>% mutate(estimate = estimate*100, 
                      std.err = std.err*100,
                      ymin = estimate -2 * std.err, 
                      ymax = estimate +2 * std.err)) + 
  aes(x = ranking, y = estimate) + 
  geom_point() + 
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .2) + 
  ylab("Risk difference (percentage point)") + xlab("") + 
  ggtitle("Average risk difference within CATE ranking") + 
  theme_minimal() + 
  geom_hline(yintercept = 0) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  ggsci::scale_color_jama()


# --- 7.5 Visualizing Treatment Effect Heterogeneity (Partial Dependence Plots) ---

# Visualize the relationship between key covariates and the predicted ITEs.
# A negative ITE from the model indicates risk reduction. We multiply by -100
# to express this as "Risk Reduction in percentage points (pp)".

plot_data <- analysis_data_complete %>%
  mutate(
    benefit_pp = ite * -100,
    has_benefit = benefit_pp > 0
  )

# Plot for SBP
p_sbp <- ggplot(plot_data, aes(x = sbp, y = benefit_pp, color = has_benefit)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_hline(yintercept = 0, color = "#AE1022") +
  scale_color_manual(values = c("TRUE" = "#09B6B3" , "FALSE" = "#EF7918")) +
  labs(
    title = "ITE vs. Systolic Blood Pressure",
    x = "Systolic Blood Pressure (mmHg)",
    y = "Predicted Risk Reduction (pp)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

p_sbp

# --- 8. IPCW (Inverse Probability of Censoring Weighting) Sensitivity Analysis ---

# To correct for potential selection bias from censoring, we use IPCW.
# Step 8a: Model the probability of being censored using the full matched cohort (before filtering).
censor_model <- glm(
  censored_3yr ~ age + sex + bmi + sbp + hba1c + egfr + hist_cvd,
  data = emulated_trial_data, # Use the full dataset before filtering out censored cases
  family = "binomial"
)

# Step 8b: Predict the probability of NOT being censored for the complete cases.
prob_not_censored <- 1 - predict(censor_model, newdata = analysis_data_complete, type = "response")

# Step 8c: Calculate the IPCW weights as the inverse of this probability.
analysis_data_complete$ipcw_weights <- 1 / prob_not_censored

# Step 8d: Re-train the causal forest on the complete cases, applying the IPCW weights.
# The components X, Y, W, and folds are the same as in the main analysis.
causal_forest_ipcw <- causal_forest(
  X = data.matrix(X), Y = Y, W = W,
  sample.weights = analysis_data_complete$ipcw_weights, # Apply IPCW weights
  num.trees = 5000,
  honesty = TRUE,
  tune.parameters = "all", # Also tune parameters for the weighted model
  clusters = folds,
  seed = 1234
)
# repeat the same analysis