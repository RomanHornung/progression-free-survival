# 0. Data, libraries and functions----
## 0.1 Loading data----
#### imputed Data
load(file = "./Data/imputed_dataset.Rda")

#### Models
load(file = "./Data/models.Rda")

## 0.2 Loading libraries----
source("./R_Code/libraries.R")

## 0.3 Loading graphic parameters----
source("./R_Code/graphic_parameters.R")

# 1. Diagnostics----
## 1.1 Proportional hazard assumption----
### 1.1.1 Test----
#### Testing the proportional hazard assumptions is consistent.
#### All co-variables should have a non-significant test.
cox.zph(cox_aic)
cox.zph(cox_bic)
#### In both cases, except for the monocytes,
#### all covariates are significant (in some cases even strongly significant),
#### so that the proportional hazard assumption may not apply.

### 1.1.2 Schoenfeld Residuals over time for seperate covariates----
#### The Schoenfeld residuals should scatter around 0 over time
#### (covariates are independent of time when proportional hazards are assumed)
#### and should have a uniform variability.

png("./Plots/schoenfeldres_AIC_model.png", width = 1800, height = 1300, res = 150)
zph <- cox.zph(cox_aic)
dimnames(zph$y)[[2]] <- c("albumin", "lactate dehydrogenase", "lymphocytes", "monocytes")
ggcoxzph(zph, point.col = gray3[1])
dev.off()


# pdf(width = pdf_w_h[1], height = pdf_w_h[2],file = "./Plots/schoenfeldres_BICmodel.pdf")
png("./Plots/schoenfeldres_BIC_model.png", width = 1800, height = 1300, res = 150)
zph <- cox.zph(cox_bic)
dimnames(zph$y)[[2]] <- c("albumin", "lactate dehydrogenase")
ggcoxzph(zph, point.col = gray3[1])
dev.off()

#### In contrast to the tests from 1.1.1, the coefficients in both models are well distributed,
#### lying at 0 on the one hand and scattering fairly evenly (except for the edges) on the other.

