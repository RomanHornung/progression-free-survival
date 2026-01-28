## ---------- Models selected using AIC and BIC ----------

# 0. Import data, libraries and functions----
## 0.1 Loading data----
load(file = "./Data/pseudo_data.RData")

df$ldh <- df$ldh/1000

## 0.2 Loading functions----
source("./R_Code/functions.R")

## 0.3 libraries----
source("./R_Code/libraries.R")


# Determine the number of complete cases:
sum(complete.cases(df))

# Events per variable:
sum(df$status=="progressive disease")
sum(df$status=="progressive disease")/(ncol(df)-2)


# Use reverse Kaplan-Meier to estimate the median follow-up time:

rev_km <- survfit(Surv(time, 1-(as.numeric(status)-1)) ~ 1, data = df)
summary(rev_km)$table["median"]/30.44

df <- as.data.frame(df)

# 1. Imputation of the data----
# simultaneous imputation and rephrasing from list to dataframe
set.seed(12345)
df_imp <- do.call("rbind", imputedata(df, m = 20))

# 2. Model
## 2.1 Calculate model without covariates and saturated----
cox_without <- coxph(Surv(time, as.numeric(status)-1) ~ 1, data = df_imp)
formula <- "Surv(time, as.numeric(status)-1) ~ leuc + lymph + n_gran + mono + e_gran + crp + albumin + protein +
ldh + mg"
cox_saturated <- coxph(as.formula(formula), data = df_imp, ties = "efron", singular.ok = FALSE)

## 2.2 Calculating step AIC----
cox_aic <- step(cox_without, direction="forward",
                scope=list(lower=cox_without, upper=cox_saturated),
                k=2*length(df_imp)/mean(!is.na(df[,!c(names(df) %in% c("status", "time"))])))

## 2.3 Calculating step BIC----
cox_bic <- step(cox_without, direction="forward",
                scope = list(lower = cox_without, upper = cox_saturated),
                k=log(nrow(df_imp)) * length(df_imp) / mean(!is.na(df[,!c(names(df) %in% c("status", "time"))])))


round(exp(cox_aic$coefficients), 3)
round(exp(cox_bic$coefficients), 3)

# 3. Saving ----
## 3.1 Models----
save(cox_aic, cox_bic, file = "./Data/models.Rda")

## 3.2 Imputed dataset----
save(df_imp, file = "./Data/imputed_dataset.Rda")






## ---------- LASSO with simple K-fold CV and train/test MI separation ----------

library(glmnet)
library(survival)

# --- Settings ---
K <- 5                  # number of CV folds
m <- 20                 # number of imputations
set.seed(12345)

# Ensure an id column exists (optional but useful)
# if (!"id" %in% names(df)) df$id <- seq_len(nrow(df))

df <- as.data.frame(df)

# There are two observations with time == 0 --> Lasso won't work.
# Solution: Add 1e-6 to all observations in 'time':
df$time <- df$time + 1e-6


# Define outcome coding once (helper to add numeric event)
.status_to_event <- function(x) as.integer(x == "progressive disease")

# Define the predictor formula for model.matrix (no intercept)
pred_formula <- ~ leuc + lymph + n_gran + mono + e_gran + crp + albumin + protein + ldh + mg

# Helper to compute C-index given a data.frame with time/event and a numeric risk score vector
cindex_from_score <- function(time, event, score) {
  # survival::survConcordance returns a list with $concordance
  survival::survConcordance(Surv(time, event) ~ score)$concordance
}

# Create CV folds on the ORIGINAL (uninmputed) data (similar to your makeCVdiv)
foldind <- makeCVdiv(data = df, yname = "status", ncv = K)

# Storage for CV performance across lambda values (filled after we learn the lambda grid)
lambda_grid <- NULL
cv_cindex_mat <- NULL  # rows = folds, cols = lambdas

# ---- Outer loop over folds ----
for (k in 1:K) {
  cat(sprintf("CV fold %d / %d\n", k, K))
  
  # --- Separate imputations on test and train ---
  # Test: list of m imputed datasets
  test_list  <- imputedata(data = df[foldind == k, ], m = m)
  # Train: stacked m imputations (single data.frame)
  train_stack <- do.call("rbind", imputedata(data = df[foldind != k, ], m = m))
  
  # Add numeric event to train
  train_stack$event <- .status_to_event(train_stack$status)
  
  # Build design matrix for glmnet (no intercept)
  x_train <- model.matrix(pred_formula, data = train_stack)[, -1]
  y_train <- Surv(train_stack$time, train_stack$event)
  
  # Weights so that each person contributes ~1 across the m stacked copies
  w_train <- rep(1 / m, nrow(train_stack))
  
  # --- Obtain (and fix) the lambda grid like glmnet would do ---
  # We take the grid from the FIRST fold's training data and reuse it for all folds
  if (is.null(lambda_grid)) {
    fit_path <- glmnet(x_train, y_train, family = "cox", alpha = 1, standardize = TRUE)
    lambda_grid <- fit_path$lambda
    # Initialize matrix to collect per-fold performance per lambda
    cv_cindex_mat <- matrix(NA_real_, nrow = K, ncol = length(lambda_grid))
    colnames(cv_cindex_mat) <- signif(lambda_grid, 5)
  }
  
  # --- Fit the full path on this fold's training data using the FIXED lambda grid ---
  fit_path_k <- glmnet(
    x_train, y_train,
    family = "cox",
    alpha  = 1,
    standardize = TRUE,
    lambda = lambda_grid,
    weights = w_train
  )
  
  # --- Evaluate on each of the m test imputations, average C-index over m ---
  # For each test imputed dataset:
  mean_cindex_per_lambda <- rep(NA_real_, length(lambda_grid))
  for (i in seq_len(m)) {
    test_i <- test_list[[i]]
    test_i$event <- .status_to_event(test_i$status)
    x_test_i <- model.matrix(pred_formula, data = test_i)[, -1]
    y_time_i <- test_i$time
    y_event_i <- test_i$event
    
    # Predict risk scores (linear predictors) for all lambdas in one go
    # predict(..., s = lambda_grid, type = "link") returns a matrix (n x length(lambda_grid))
    eta_mat <- predict(fit_path_k, newx = x_test_i, s = lambda_grid, type = "link")
    # If only one lambda, coerce to matrix
    if (is.null(dim(eta_mat))) eta_mat <- matrix(eta_mat, ncol = 1)
    
    # Compute C-index for each lambda and accumulate
    ci_vec <- apply(eta_mat, 2, function(score) cindex_from_score(y_time_i, y_event_i, as.numeric(score)))
    
    if (i == 1) {
      mean_cindex_per_lambda <- ci_vec
    } else {
      mean_cindex_per_lambda <- mean_cindex_per_lambda + ci_vec
    }
  }
  mean_cindex_per_lambda <- mean_cindex_per_lambda / m
  
  # Store per-fold results
  cv_cindex_mat[k, ] <- mean_cindex_per_lambda
}

# ---- Aggregate across folds and choose best lambda ----
mean_cv_cindex <- colMeans(cv_cindex_mat, na.rm = TRUE)
best_lambda <- lambda_grid[ which.max(mean_cv_cindex) ]
cat("Best lambda (by mean CV C-index):", best_lambda, "\n")


# ---- Fit FINAL LASSO model on the full data (stacked MI on the entire dataset) ----
# Impute full data m times and stack (consistent with your workflow)
full_stack <- do.call("rbind", imputedata(data = df, m = m))
full_stack$event <- .status_to_event(full_stack$status)
x_full <- model.matrix(pred_formula, data = full_stack)[, -1]
y_full <- Surv(full_stack$time, full_stack$event)
w_full <- rep(1 / m, nrow(full_stack))

lasso_final <- glmnet(
  x_full, y_full,
  family = "cox",
  alpha = 1,
  standardize = TRUE,
  lambda = best_lambda,
  weights = w_full
)

# Extract selected variables
coef_final <- as.matrix(coef(lasso_final))
sel_vars <- rownames(coef_final)[as.numeric(coef_final) != 0]

round(exp(coef_final[sel_vars,]), 3)



# Hazard Ratios from both models
hr_lasso <- round(exp(coef_final[sel_vars, ]), 3)
hr_aic   <- round(exp(cox_aic$coefficients), 3)

# Merge both into a single data frame
comparison <- merge(
  data.frame(variable = names(hr_lasso), HR_LASSO = as.numeric(hr_lasso)),
  data.frame(variable = names(hr_aic), HR_AIC = as.numeric(hr_aic)),
  by = "variable",
  all = TRUE
)

# Optionally sort by variable name
comparison <- comparison[order(comparison$variable), ]

# Display
comparison



# --- Report and save ---
print(list(
  best_lambda = best_lambda,
  selected_variables = sel_vars,
  mean_cv_cindex = mean(mean_cv_cindex),
  cv_cindex_path = mean_cv_cindex
))

save(best_lambda, lasso_final, sel_vars, mean_cv_cindex, cv_cindex_mat, lambda_grid,
     file = "./Data/lasso_model.Rda")
