setwd("/nfsmb/koll/hornung/Projects/SideProjects/ProjectDormanetal/Revision/GitHub_nach_Revision/progression-free-survival")

## ===== Repeated nested CV for LASSO–Cox with MI-separation (C-index + IBS) =====

# 0. Import data, libraries and functions----
## 0.1 Loading data----
load(file = "./Data/pseudo_data.RData")

df$ldh <- df$ldh/1000

## 0.2 Loading functions----
source("./R_Code/functions.R")

## 0.3 libraries----
source("./R_Code/libraries.R")

library(survival)
library(glmnet)
library(pec)

## --- Settings ---
K_outer <- 5          # outer CV folds
K_inner <- 5          # inner CV folds for lambda tuning
R <- 10               # repetitions of the entire outer CV
m <- 20               # number of imputations
set.seed(12345)

## --- Data prep ---
df <- as.data.frame(df)
if (!"id" %in% names(df)) df$id <- seq_len(nrow(df))  # stable id for alignment
df$status <- factor(df$status, levels = c("survival","progressive disease"))
df$time <- df$time + 1e-6     # ensure strictly positive times

.status_to_event <- function(x) as.integer(x == "progressive disease")

## --- Helpers ---

# Robust C-index from linear predictor (checks lengths and finiteness)
cindex_from_score <- function(time, event, score) {
  n <- length(time)
  if (length(event) != n || length(score) != n) {
    stop(sprintf("Length mismatch: time=%d, event=%d, score=%d", n, length(event), length(score)))
  }
  idx <- which(is.finite(time) & is.finite(event) & is.finite(score))
  survival::survConcordance(Surv(time[idx], event[idx]) ~ score[idx])$concordance
}

# Build model matrix with fixed columns (align test to train), no row dropping
build_X <- function(data, terms_obj, ref_colnames = NULL) {
  mf <- model.frame(terms_obj, data = data, na.action = na.pass)  # do not drop rows
  X <- model.matrix(terms_obj, data = mf)[, -1, drop = FALSE]     # drop intercept
  # Align to reference columns if provided
  if (!is.null(ref_colnames)) {
    missing <- setdiff(ref_colnames, colnames(X))
    if (length(missing)) {
      Xmiss <- matrix(0, nrow = nrow(X), ncol = length(missing))
      colnames(Xmiss) <- missing
      X <- cbind(X, Xmiss)
    }
    extra <- setdiff(colnames(X), ref_colnames)
    if (length(extra)) X <- X[, setdiff(colnames(X), extra), drop = FALSE]
    X <- X[, ref_colnames, drop = FALSE]
  }
  stopifnot(nrow(X) == nrow(data))  # critical: ensure no row dropping
  X
}

# Predictor specification
pred_formula <- ~ leuc + lymph + n_gran + mono + e_gran + crp + albumin + protein + ldh + mg
mm_terms <- terms(pred_formula)

## --- Storage for repeated CV results ---
cindex_lasso <- brier_lasso <- numeric(0)

## --- Repeated outer CV ---
for (rep_ix in seq_len(R)) {
  message(sprintf("=== Repetition %d / %d ===", rep_ix, R))
  
  # Fresh outer folds each repetition (stratified by status)
  set.seed(1e6 + rep_ix)
  fold_outer <- makeCVdiv(data = df, yname = "status", ncv = K_outer)
  
  # Temporary collectors for this repetition
  cindex_outer <- brier_outer <- rep(NA_real_, K_outer)
  
  for (k in seq_len(K_outer)) {
    message(sprintf("  Outer fold %d / %d", k, K_outer))
    
    # Split original data
    df_train_outer <- df[fold_outer != k, , drop = FALSE]
    df_test_outer  <- df[fold_outer == k, , drop = FALSE]
    
    # --- MI for outer train/test ---
    # Train: stack m imputations
    train_stack_outer <- do.call("rbind", imputedata(data = df_train_outer, m = m))
    train_stack_outer <- train_stack_outer[order(train_stack_outer$id), , drop = FALSE]
    train_stack_outer$event <- .status_to_event(train_stack_outer$status)
    
    # Test: list of m imputations
    test_list_outer <- imputedata(data = df_test_outer, m = m)
    test_list_outer <- lapply(test_list_outer, function(d) d[order(d$id), , drop = FALSE])
    
    # Build X/y for outer training
    x_train_outer <- build_X(train_stack_outer, mm_terms)
    ref_cols <- colnames(x_train_outer)
    y_train_outer <- Surv(train_stack_outer$time, train_stack_outer$event)
    w_train_outer <- rep(1 / m, nrow(train_stack_outer))
    
    # --- INNER CV for lambda tuning (on df_train_outer only) ---
    set.seed(2e6 + rep_ix*10 + k)
    fold_inner <- makeCVdiv(data = df_train_outer, yname = "status", ncv = K_inner)
    
    lambda_grid <- NULL
    cvmat_inner <- NULL  # rows = inner folds, cols = lambdas (C-index)
    
    for (kk in seq_len(K_inner)) {
      # Split inner train/test (still on original outer-train)
      dtr_in  <- df_train_outer[fold_inner != kk, , drop = FALSE]
      dte_in  <- df_train_outer[fold_inner == kk, , drop = FALSE]
      
      # MI separation for inner train/test
      train_stack_in <- do.call("rbind", imputedata(data = dtr_in, m = m))
      test_list_in   <- imputedata(data = dte_in, m = m)
      
      train_stack_in <- train_stack_in[order(train_stack_in$id), , drop = FALSE]
      test_list_in   <- lapply(test_list_in, function(d) d[order(d$id), , drop = FALSE])
      
      train_stack_in$event <- .status_to_event(train_stack_in$status)
      
      x_in <- build_X(train_stack_in, mm_terms)
      ref_cols_in <- colnames(x_in)
      y_in <- Surv(train_stack_in$time, train_stack_in$event)
      w_in <- rep(1 / m, nrow(train_stack_in))
      
      # Define lambda grid on the first inner fold; reuse for others
      if (is.null(lambda_grid)) {
        fit_path <- glmnet(x_in, y_in, family = "cox", alpha = 1, standardize = TRUE)
        lambda_grid <- fit_path$lambda
        cvmat_inner <- matrix(NA_real_, nrow = K_inner, ncol = length(lambda_grid))
        colnames(cvmat_inner) <- signif(lambda_grid, 5)
      }
      
      # Fit glmnet path for fixed lambda grid on inner-train
      fit_path_in <- glmnet(x_in, y_in, family = "cox", alpha = 1,
                            standardize = TRUE, lambda = lambda_grid, weights = w_in)
      
      # Evaluate on each of the m inner-test imputations and average C-index
      mean_cindex_lambda <- rep(0, length(lambda_grid))
      for (i in seq_len(m)) {
        dte_i <- test_list_in[[i]]
        dte_i$event <- .status_to_event(dte_i$status)
        x_te_i <- build_X(dte_i, mm_terms, ref_colnames = ref_cols_in)
        stopifnot(nrow(x_te_i) == nrow(dte_i))  # critical check
        
        # Predict linear predictors for all lambdas
        eta_mat <- predict(fit_path_in, newx = x_te_i, s = lambda_grid, type = "link")
        if (is.null(dim(eta_mat))) eta_mat <- matrix(eta_mat, ncol = 1)
        
        ci_vec <- apply(eta_mat, 2, function(sc)
          cindex_from_score(dte_i$time, dte_i$event, as.numeric(sc))
        )
        mean_cindex_lambda <- mean_cindex_lambda + ci_vec
      }
      mean_cindex_lambda <- mean_cindex_lambda / m
      
      cvmat_inner[kk, ] <- mean_cindex_lambda
    } # end inner folds
    
    # Choose best lambda by mean inner-C-index
    mean_inner <- colMeans(cvmat_inner, na.rm = TRUE)
    lambda_best <- lambda_grid[ which.max(mean_inner) ]
    
    # --- Fit final outer-train model at tuned lambda ---
    fit_lasso_outer <- glmnet(
      x_train_outer, y_train_outer,
      family = "cox", alpha = 1, standardize = TRUE,
      lambda = lambda_best, weights = w_train_outer
    )
    
    # Prepare an offset-based coxph on outer-train for IBS computation
    lp_train_outer <- as.numeric(predict(fit_lasso_outer, newx = x_train_outer, type = "link"))
    stopifnot(length(lp_train_outer) == nrow(train_stack_outer))
    train_stack_outer$lp <- lp_train_outer
    cox_offset_outer <- coxph(Surv(time, event) ~ offset(lp),
                              data = train_stack_outer, ties = "efron", x = TRUE)
    
    # --- Evaluate on outer-test (average over m imputations) ---
    cidx_vec <- brier_vec <- numeric(m)
    
    for (i in seq_len(m)) {
      dte_o <- test_list_outer[[i]]
      dte_o$event <- .status_to_event(dte_o$status)
      
      x_te <- build_X(dte_o, mm_terms, ref_colnames = ref_cols)
      stopifnot(nrow(x_te) == nrow(dte_o))  # critical check
      lp_te <- as.numeric(predict(fit_lasso_outer, newx = x_te, type = "link"))
      stopifnot(length(lp_te) == nrow(dte_o))
      
      cidx_vec[i] <- cindex_from_score(dte_o$time, dte_o$event, lp_te)
      
      # IBS via pec on the offset-cox trained in outer-train:
      dte_o$lp <- lp_te
      pe <- suppressMessages(suppressWarnings(
        pec::pec(cox_offset_outer,
                 formula = Surv(time, event) ~ 1,
                 data = dte_o, exact = FALSE)
      ))
      brier_vec[i] <- as.numeric(pec::crps(pe))[2]
    }
    
    cindex_outer[k] <- mean(cidx_vec)
    brier_outer[k]  <- mean(brier_vec)
  } # end outer folds
  
  # Collect repetition results
  cindex_lasso <- c(cindex_lasso, cindex_outer)
  brier_lasso  <- c(brier_lasso,  brier_outer)
}

## --- Summary of repeated CV performance ---
perf_lasso <- data.frame(
  cindex = cindex_lasso,
  ibrier = brier_lasso
)

print(summary(perf_lasso$cindex))
print(summary(perf_lasso$ibrier))

save(perf_lasso, file = "./Data/lasso_repeatedCV_performance.Rda")
