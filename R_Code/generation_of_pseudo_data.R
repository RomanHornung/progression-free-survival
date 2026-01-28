# 0. Data, libraries----
## 0.1 Loading Data----
load(file = "./Data/prep.RData")

## 0.2 Loading libraries----
source("./R_Code/libraries.R")

# 1. Generation of pseudo data----
## 1.1 metric parameters----
set.seed(12345)
df_pseudo <- sim_df(df, n = nrow(df))[-1]

#### Gammadistribution for time
temp <- egamma(df$time)$parameters
df_pseudo$time <- summary(rgamma(nrow(df), shape = temp[1], scale = temp[2]))
rm(temp)

## 1.2 Categorical parameters----
df_pseudo$status <- sample(c("survival", "progressive disease"), nrow(df), replace = T)
df_pseudo$status <- factor(df_pseudo$status)

## 1.3 Production of NAs in pseudo data----
df_pseudo$leuc[which(is.na(df$leuc))] <- NA
df_pseudo$lymph[which(is.na(df$lymph))] <- NA
df_pseudo$n_gran[which(is.na(df$n_gran))] <- NA
df_pseudo$mono[which(is.na(df$mono))] <- NA
df_pseudo$e_gran[which(is.na(df$e_gran))] <- NA
df_pseudo$crp[which(is.na(df$crp))] <- NA
df_pseudo$albumin[which(is.na(df$albumin))] <- NA
df_pseudo$protein[which(is.na(df$protein))] <- NA
df_pseudo$ldh[which(is.na(df$ldh))] <- NA
df_pseudo$mg[which(is.na(df$mg))] <- NA

#### removing df and renaming df_pseudo to simply go through the entire code with the pseudo data set
rm(df)
df <- df_pseudo

## 2. Saving----
save(df, file = "./Data/pseudo_data.RData")
