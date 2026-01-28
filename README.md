# R Code for a Statistical Analysis for Prediction of Progression-Free Survival under Immune Checkpoint Inhibition

This repository contains the R code used for the development and validation of Cox proportional hazards models to predict progression-free survival (PFS) in patients treated with immune checkpoint inhibitors (ICI).  
The analysis includes multiple imputation for missing data, model building via forward selection based on AIC and BIC criteria, and performance evaluation using repeated cross-validation, the concordance index, and the Brier score. In addition, a LASSO-penalized Cox regression model was fitted as a sensitivity analysis.  
Model diagnostics were performed using plots of scaled Schoenfeld residuals and calibration plots.

The results generated with this code are part of a retrospective study published in the following paper:

> **Dorman, K., Breitenwieser, K., Fischer, L. et al.**  
> *Real-world analysis of immune checkpoint inhibitor efficacy and response predictors in patients treated at the CCCMunichLMU outpatient clinic.*  
> **Scientific Reports**, 15, 43269 (2025).  
> https://doi.org/10.1038/s41598-025-30220-0

The original patient data cannot be included in this repository.  
Instead, pseudo data with similar statistical properties to the original data are provided. Consequently, numerical results obtained using the pseudo data will differ from those reported in the publication, while the analysis workflow remains reproducible.

The structure of the data and the code is as follows:

- **R Code**
  - **calibration.R**: Required data: *imputed_dataset.Rda*, *models.Rda*, *lasso_model.Rda*. Calculation of the calibration plots. 
  - **functions.R**: Required data: none. Functions for imputing data and index formation for cross-validation.
  - **generation_of_pseudo_data.R**: Required data: original data (not included here). Generation of a pseudo data set from the original data
  - **graphic_parameters.R**: Required data: none. Colors, heights and widths for generating graphics.
  - **libraries.R**: Required data: none. All packages used, for fast loading.
  - **model.R**: Required data: *pseudo_data.RData*. Calculation of the models with step-wise selection according to AIC and BIC and the LASSO-penalized Cox regression model.
  - **preparation.R**: Required data: none. Preparation of the original data (not included here).
  - **schoenfeld.R**: Required data: *imputed_dataset.Rda*, *models.Rda*. Calculation of the Schoenfeld residual plots.
  - **validation.R**: Required data: *pseudo_data.RData*. Validation of the models with step-wise selection according to AIC and BIC.
  - **validation_lasso.R**: Required data: *pseudo_data.RData*. Validation of the LASSO-penalized Cox regression model.

- **Data**:
  - **imputed dataset.Rda**: (*df_imp*) 20-fold stacked, imputed data set from the pseudo data
  - **models.Rda**: Original Coxph models, selected using AIC (*cox_aic*) and BIC (*cox_bic*)
  - **lasso_model.Rda**: the LASSO-penalized Cox regression model
  - **prep.RData**: prepared dataset from pseudo data (*df*)
  - **pseudo_data.RData**: generated pseudo data set from the original data (*df*)
  - **session_info.Rda**: Version information of R and packages
  - **validation.Rda**: From original data for cross-validation: optimal cutpoint for the AIC (*opt_cut_aic*) and BIC (*opt_cut_bic*) model, survival plots for the cutpoints for the AIC (*surv_val_aic*; plot is included in Plots/Cutoffs_validated_AICmodel.pdf) and BIC (*surv_val_bic*; plot is included in Plots/Cutoffs_validated_BICmodel.pdf) model, Brier and Cindex with respect to the AIC and BIC model (*performance_plot*; plot is included in Plots/Performance_validated_cox.pdf) and the corresponding values (res).
  - **lasso_repeatedCV_performance.Rda**: C-index and Brier score values from the repeated cross-validation for estimating the predictive performance of the LASSO-penalized Cox regression model

- **Plots**
  - **calibration_cox_models.png**: Calibration plots for the AIC model, the BIC model, and the LASSO-penalized Cox regression model (Graphic is calculated in *R Code/calibration.R*)
  - **schoenfeldres_AICmodel.png**: Schoenfeld residuals for the original AIC model (Graphic is calculated in *R Code/ schoenfeld.R*)
  - **schoenfeldres_BICmodel.png**: Schoenfeld residuals for the original BIC model (Graphic is calculated in *R Code/schoenfeld.R*)
