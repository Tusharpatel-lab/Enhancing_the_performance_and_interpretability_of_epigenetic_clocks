# Author: Tushar Patel (Tushar.Patel@leibniz-fli.de)
# Date: 29-01-2026
# Description: Run TFMethyl clock
# Tested with: R version 4.4.0; data.table 1.17.6; glmnet 4.1-9

predict_age_from_beta <- function(
    beta_matrix,
    model,
    cluster_info,
    aggregate_fun = median
) {
  ## ------------------------------------------------------------
  ## Args:
  ##   beta_matrix   : data.frame or data.table
  ##                   Rows = CpGs, columns = samples
  ##   model         : trained model object with predict() method
  ##   cluster_info  : data.frame with:
  ##                   - rownames = CpG IDs
  ##                   - column 'cluster' defining CpG clusters
  ##   aggregate_fun : function used to aggregate CpGs per cluster
  ##                   (default = median)
  ##
  ## Returns:
  ##   Numeric vector of predicted chronological age (years),
  ##   named by sample ID
  ## ------------------------------------------------------------
  
  ## Coerce to data.frame for safety
  beta_matrix <- as.data.frame(beta_matrix)
  
  ## -----------------------------
  ## Input validation
  ## -----------------------------
  if (is.null(rownames(beta_matrix))) {
    stop("beta_matrix must have CpG IDs as rownames.")
  }
  
  if (is.null(rownames(cluster_info))) {
    stop("cluster_info must have CpG IDs as rownames.")
  }
  
  if (!"cluster" %in% colnames(cluster_info)) {
    stop("cluster_info must contain a column named 'cluster'.")
  }
  
  #finding CpGs with non-zero coefficients
  coefficients <- coef(model, s = model$lambda.1se)
  coefficients <- as.data.frame(as.matrix(coefficients))
  coefficients <- coefficients[coefficients[, 1] != 0, , drop = FALSE]
  coefficients_names <- gsub("^cluster ", "", rownames(coefficients))
  cluster_info_true <- cluster_info[cluster_info$cluster %in% coefficients_names, , drop = FALSE]
  
  ## -----------------------------
  ## Subset and order CpGs
  ## -----------------------------
  common_cpgs <- intersect(rownames(beta_matrix), rownames(cluster_info_true))
  
  if (length(common_cpgs) == 0) {
    stop("No overlapping of probes between samples and required CpGs.")
  }
  
  beta_sub <- beta_matrix[common_cpgs,
                          colSums(is.na(beta_matrix[common_cpgs, , drop = FALSE])) == 0,
                          drop = FALSE]
  
  ## Ensure identical ordering
  beta_sub <- beta_sub[rownames(cluster_info), , drop = FALSE]
  beta_sub[is.na(beta_sub)] <- 0
  rownames(beta_sub) <- rownames(cluster_info)
  
  ## -----------------------------
  ## Aggregate CpGs by cluster
  ## -----------------------------
  beta_clustered <- cbind(cluster = cluster_info$cluster, beta_sub)
  
  beta_clustered <- aggregate(
    . ~ cluster,
    data = beta_clustered,
    FUN = aggregate_fun
  )
  
  rownames(beta_clustered) <- beta_clustered$cluster
  beta_clustered$cluster <- NULL
  
  ## -----------------------------
  ## Predict age (back-transform)
  ## -----------------------------
  y_pred <- predict(model, t(beta_clustered))
  y_pred <- as.numeric(y_pred)
  
  # back-transform from log(age + 1)
  y_pred <- exp(y_pred) - 1
  
  names(y_pred) <- colnames(beta_sub)
  
  return(y_pred)
}

#------------------------------------------------------------------------------

# load essential libraries
library(data.table)
library(glmnet)

#load the test methylation beta matrix, rows = CpGs, columns = Samples
beta_test <- fread("./beta_test_matrix.csv")
beta_test <- as.data.frame(beta_test)
rownames(beta_test) <- beta_test[[1]]
beta_test[[1]] <- NULL

#beta_test structure: 
#              sample1  sample2  sample3  sample4  sample5
# cg00000957      0.85     0.84     0.84     0.83     0.84
# cg00002028      0.10     0.11     0.07     0.10     0.12
# cg00002719      0.07     0.08     0.06     0.07     0.04
# cg00002837      0.32     0.31     0.34     0.38     0.30
# cg00003287      0.17     0.18     0.19     0.14     0.15

#load the trained model, and the clustering information
model <- readRDS("./TFMethyl_model.rds")
cluster_info <- readRDS("./cluster_info.rds")

#the TFMethyl clock prediction function
#the predicted chronological-age (in years) 
predicted_age <- predict_age_from_beta(
  beta_matrix  = beta_test,
  model        = model,
  cluster_info = cluster_info)