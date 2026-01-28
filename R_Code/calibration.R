# ============================================================
# Model calibration (AIC, BIC, and LASSO)
# ============================================================

# 0) Load data and models ----
load("./Data/imputed_dataset.Rda")   # df_imp
load("./Data/models.Rda")            # cox_aic, cox_bic
load("./Data/lasso_model.Rda")       # lasso_final, best_lambda, sel_vars

source("./R_Code/libraries.R")
source("./R_Code/graphic_parameters.R")

library(survival)
library(ggplot2)

# Weight each observation by 1/M since data are stacked imputations
M <- 20
w_imp <- rep(1 / M, nrow(df_imp))

# ------------------------------------------------------------
# 1. Calibration slopes
# ------------------------------------------------------------

# Add linear predictors
df_imp$lp_aic   <- predict(cox_aic, type = "lp")
df_imp$lp_bic   <- predict(cox_bic, type = "lp")

# Build LP for LASSO (design matrix as during fitting)
pred_formula <- ~ leuc + lymph + n_gran + mono + e_gran + crp + albumin + protein + ldh + mg
X_imp <- model.matrix(pred_formula, data = df_imp)[, -1, drop = FALSE]
df_imp$lp_lasso <- as.numeric(predict(lasso_final, newx = X_imp, type = "link"))

# Fit slope models (slope ≈ 1 = well calibrated)
cal_aic   <- coxph(Surv(time, as.numeric(status) - 1) ~ lp_aic,
                   data = df_imp, weights = w_imp)
cal_bic   <- coxph(Surv(time, as.numeric(status) - 1) ~ lp_bic,
                   data = df_imp, weights = w_imp)
cal_lasso <- coxph(Surv(time, as.numeric(status) - 1) ~ lp_lasso,
                   data = df_imp, weights = w_imp)

summary(cal_aic)
summary(cal_bic)
summary(cal_lasso)


# ------------------------------------------------------------
# 2. Calibration plots (6, 12, 24 months)
# ------------------------------------------------------------

cal_times <- c(6, 12, 24) * 30.44  # months -> days
ng <- 5  # number of risk groups

calib_df <- function(fit, data, lp_col, timepoint, model_label) {
  # Predict survival at given time using baseline survival
  base <- survfit(fit)
  t_idx <- max(which(base$time <= timepoint))
  base_surv_t <- base$surv[t_idx]
  
  # predicted individual survival and risk
  S_pred <- base_surv_t ^ exp(data[[lp_col]])
  risk_pred <- 1 - S_pred
  data$risk_pred <- risk_pred
  
  # Group by predicted risk
  data$group <- cut(risk_pred,
                    breaks = quantile(risk_pred, probs = seq(0, 1, length.out = ng + 1)),
                    include.lowest = TRUE, labels = FALSE)
  
  # Observed risk per group from KM
  obs <- tapply(1:nrow(data), data$group, function(idx) {
    sf <- survfit(Surv(time, as.numeric(status) - 1) ~ 1, data = data[idx, ], weights = w_imp[idx])
    idx_t <- max(which(sf$time <= timepoint))
    obs_risk <- 1 - ifelse(length(idx_t) == 0, 1, sf$surv[idx_t])
    obs_risk
  })
  
  data.frame(
    model = model_label,
    time = timepoint / 30.44,
    group = 1:ng,
    pred = tapply(risk_pred, data$group, mean),
    obs = as.numeric(obs)
  )
}

# Combine data for all models and times
df_cal_all <- do.call(rbind, c(
  lapply(cal_times, function(t) calib_df(cox_aic, df_imp, "lp_aic", t, "AIC model")),
  lapply(cal_times, function(t) calib_df(cox_bic, df_imp, "lp_bic", t, "BIC model")),
  lapply(cal_times, function(t) calib_df(coxph(Surv(time, as.numeric(status)-1) ~ offset(lp_lasso),
                                               data = df_imp, weights = w_imp, ties = "efron"),
                                         df_imp, "lp_lasso", t, "LASSO model"))
))

# Calibration plot
p_cal <- ggplot(df_cal_all, aes(x = pred, y = obs, color = model, shape = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray40") +
  geom_point(size = 2) +
  facet_wrap(
    ~ time,
    nrow = 1,
    labeller = as_labeller(function(x)
      paste0(x, " Months"))
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Predicted risk of progression/death",
       y = "Observed risk (Kaplan–Meier estimate)",
       color = "Model", shape = "Model") +
  theme_bw()

print(p_cal)

# Save
ggsave("./Plots/calibration_cox_models.png", p_cal, width = 11, height = 3.5, dpi = 300)




# how broad are the predicted risks really?
check_spread <- function(x, name) {
  cat("\n", name, "\n")
  print(quantile(x, c(.05, .25, .5, .75, .95)))
  cat("range:", range(x), "\n")
}

# AIC
base_aic <- survfit(cox_aic)
t12 <- 12 * 30.44
s0_12 <- base_aic$surv[max(which(base_aic$time <= t12))]
risk_aic <- 1 - s0_12 ^ exp(df_imp$lp_aic)
check_spread(risk_aic, "AIC")

# BIC
base_bic <- survfit(cox_bic)
s0_12b <- base_bic$surv[max(which(base_bic$time <= t12))]
risk_bic <- 1 - s0_12b ^ exp(df_imp$lp_bic)
check_spread(risk_bic, "BIC")

# LASSO
fit_lasso_cal <- coxph(
  Surv(time, as.numeric(status)-1) ~ offset(lp_lasso),
  data = df_imp, weights = w_imp, ties = "efron"
)
base_lasso <- survfit(fit_lasso_cal)
s0_12l <- base_lasso$surv[max(which(base_lasso$time <= t12))]
risk_lasso <- 1 - s0_12l ^ exp(df_imp$lp_lasso)
check_spread(risk_lasso, "LASSO")
