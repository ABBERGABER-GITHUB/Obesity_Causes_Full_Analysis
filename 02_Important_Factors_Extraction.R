# ============================
# Important Factors Extraction 
# ============================

# 0) Paths
data_path  <- "D:/Data Analysis/Finished Projects/Obesity_Causes_Analysis/Cleaned_Final_Data.csv"
output_dir <- "D:/Data Analysis/Finished Projects/Obesity_Causes_Analysis/Important_Factors"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 1) Required packages
req <- c("dplyr","broom","randomForest","glmnet","ggplot2","scales","forcats","caret","ggrepel")
new <- req[!req %in% installed.packages()[,"Package"]]
if(length(new)) install.packages(new, repos = "https://cloud.r-project.org")
lapply(req, library, character.only = TRUE)

# 2) Read data & basic cleaning
data <- read.csv(data_path, stringsAsFactors = FALSE)
data$Country <- stringr::str_squish(as.character(data$Country))
data$Income_Level <- stringr::str_squish(as.character(data$Income_Level))
num_cols <- setdiff(names(data), c("Country","Income_Level"))
for(col in num_cols){
  data[[col]] <- suppressWarnings(as.numeric(as.character(data[[col]])))
}

# 3) Target variable
target <- "Avg_Adult_Obesity_BMI30"

# 4) Create output folder
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Create folder
imp_dir <- file.path(output_dir, "Important_Factors")
if(!dir.exists(imp_dir)) dir.create(imp_dir, recursive = TRUE)

# 0) Required libs (install if missing)
req <- c("dplyr","broom","randomForest","glmnet","ggplot2","scales","forcats","caret","ggrepel")
new <- req[!req %in% installed.packages()[,"Package"]]
if(length(new)) install.packages(new, repos = "https://cloud.r-project.org")
lapply(req, library, character.only = TRUE)
library(broom); library(ggplot2); library(dplyr); library(forcats); library(scales)

# 1) Correlations with target (pairwise.complete.obs)
num_vars <- data %>% dplyr::select(where(is.numeric)) %>% names()
num_vars <- setdiff(num_vars, target)  # exclude target from predictors list
corrs <- sapply(num_vars, function(v) cor(data[[v]], data[[target]], use = "pairwise.complete.obs"))
corr_df <- data.frame(variable = names(corrs), correlation_with_target = as.numeric(corrs))
corr_df <- corr_df %>% arrange(desc(abs(correlation_with_target)))
write.csv(corr_df, file.path(imp_dir, "correlation_with_target.csv"), row.names = FALSE)

# Plot correlations (top 15 by absolute corr)
top_corr <- corr_df %>% slice_max(order_by = abs(correlation_with_target), n = 15)
p_corr <- ggplot(top_corr, aes(x = reorder(variable, correlation_with_target), y = correlation_with_target)) +
  geom_col(aes(fill = correlation_with_target > 0)) +
  coord_flip() +
  labs(title = "Top 15 variables by absolute correlation with Adult Obesity",
       x = "", y = "Pearson correlation") +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "steelblue"), guide = FALSE) +
  theme_minimal(base_size = 12)
ggsave(file.path(imp_dir, "top15_correlation_with_target.png"), p_corr, width = 8, height = 6)

# 2) Standardized linear regression (standardize predictors -> get standardized coeff)
# Prepare dataframe with complete cases for chosen candidate predictors (use those with non-NA correlation)
candidates <- top_corr$variable  # top correlated variables (up to 15)
# ensure at least some predictors exist
candidates <- candidates[candidates %in% names(data)]
lm_df <- data %>% dplyr::select(all_of(c(target, candidates))) %>% drop_na()

if(nrow(lm_df) >= 10 && length(candidates) >= 1){
  # Standardize predictors and target
  std_df <- lm_df
  std_df[candidates] <- scale(std_df[candidates])
  std_df[[target]] <- scale(std_df[[target]])  # standardizing target permits direct comparison of betas
  
  # Fit linear model on standardized data
  lm_formula <- as.formula(paste(target, "~", paste(candidates, collapse = " + ")))
  lm_std <- lm(lm_formula, data = std_df)
  lm_tidy <- broom::tidy(lm_std) %>% filter(term != "(Intercept)")
  lm_tidy <- lm_tidy %>% arrange(desc(abs(estimate)))
  write.csv(lm_tidy, file.path(imp_dir, "standardized_lm_coefficients.csv"), row.names = FALSE)
  
  # Plot top coefficients
  top_coefs <- lm_tidy %>% slice_max(order_by = abs(estimate), n = 15)
  p_coef <- ggplot(top_coefs, aes(x = reorder(term, estimate), y = estimate, fill = estimate > 0)) +
    geom_col() + coord_flip() +
    labs(title = "Standardized coefficients (linear model)", x = "", y = "Standardized beta") +
    scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "steelblue"), guide = FALSE) +
    theme_minimal(base_size = 12)
  ggsave(file.path(imp_dir, "standardized_lm_top_coefs.png"), p_coef, width = 8, height = 6)
} else {
  message("not enough data to build a standard linear model (standardized).")
}

# 3) Random Forest variable importance
# Prepare RF data (complete cases)
rf_vars <- unique(c(target, candidates))
rf_df <- data %>% dplyr::select(any_of(rf_vars)) %>% drop_na()
if(nrow(rf_df) >= 30 && length(candidates) >= 2){
  set.seed(123)
  rf_fit <- randomForest::randomForest(as.formula(paste(target, "~ .")), data = rf_df, importance = TRUE, ntree = 500)
  rf_imp <- importance(rf_fit, type = 1) # %IncMSE
  rf_imp_df <- data.frame(variable = rownames(rf_imp), IncMSE = rf_imp[,1])
  rf_imp_df <- rf_imp_df %>% arrange(desc(IncMSE))
  write.csv(rf_imp_df, file.path(imp_dir, "rf_variable_importance_IncMSE.csv"), row.names = FALSE)
  
  # Plot top 15 RF importance
  p_rf <- ggplot(rf_imp_df %>% slice_max(order_by = IncMSE, n = 15), aes(x = reorder(variable, IncMSE), y = IncMSE)) +
    geom_col(fill = "tomato") + coord_flip() +
    labs(title = "Random Forest: variable importance (%IncMSE)", x = "", y = "%IncMSE") +
    theme_minimal(base_size = 12)
  ggsave(file.path(imp_dir, "rf_top15_importance.png"), p_rf, width = 8, height = 6)
  
  # Save RF model summary
  sink(file.path(imp_dir, "rf_model_summary.txt"))
  print(rf_fit)
  sink()
} else {
  message("بيانات غير كافية لتشغيل Random Forest (احتاج >=30 صفًا و >=2 متغيرات).")
}

# 4) Lasso (glmnet) selection
# Prepare matrix
glmnet_vars <- candidates  # use same candidate predictors
glmnet_df <- data %>% dplyr::select(any_of(c(target, glmnet_vars))) %>% drop_na()
if(nrow(glmnet_df) >= 30 && length(glmnet_vars) >= 2){
  x <- model.matrix(as.formula(paste("~", paste(glmnet_vars, collapse = "+"))), data = glmnet_df)[,-1]
  y <- glmnet_df[[target]]
  set.seed(123)
  cv_lasso <- cv.glmnet(x, y, alpha = 1, nfolds = 5)
  best_lam <- cv_lasso$lambda.min
  fit_lasso <- glmnet::glmnet(x, y, alpha = 1, lambda = best_lam)
  coefs <- as.matrix(coef(fit_lasso))
  coefs_df <- data.frame(variable = rownames(coefs), coefficient = as.numeric(coefs))
  coefs_df <- coefs_df %>% filter(variable != "(Intercept)") %>% arrange(desc(abs(coefficient)))
  write.csv(coefs_df, file.path(imp_dir, "lasso_coefficients_lambda_min.csv"), row.names = FALSE)
  
  # Non-zero coefficients (selected variables)
  nonzero <- coefs_df %>% filter(coefficient != 0)
  write.csv(nonzero, file.path(imp_dir, "lasso_nonzero_coefficients.csv"), row.names = FALSE)
  
  # Plot non-zero coefficients
  if(nrow(nonzero) > 0){
    p_lasso <- ggplot(nonzero, aes(x = reorder(variable, coefficient), y = coefficient, fill = coefficient > 0)) +
      geom_col() + coord_flip() +
      labs(title = paste("Lasso selected coefficients (lambda.min =", round(best_lam,6),")"), x = "", y = "Coefficient") +
      scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "steelblue"), guide = FALSE) +
      theme_minimal(base_size = 12)
    ggsave(file.path(imp_dir, "lasso_selected_coeffs.png"), p_lasso, width = 8, height = 6)
  }
  
  # Save CV plot
  png(file.path(imp_dir, "cv_lasso_plot.png"), width = 800, height = 600); plot(cv_lasso); dev.off()
} else {
  message("بيانات غير كافية لتشغيل Lasso (احتاج >=30 صفًا و >=2 متغيرات).")
}

# 5) Combine importance ranks into one table
# We'll combine correlation rank, standardized lm rank (if exists), RF rank, Lasso selection
rank_tab <- corr_df %>% rename(corr = correlation_with_target) %>% mutate(corr_rank = rank(-abs(corr)))
if(exists("lm_tidy")) {
  lm_rank <- lm_tidy %>% mutate(lm_rank = rank(-abs(estimate))) %>% dplyr::select(term, lm_rank) %>% rename(variable = term)
  rank_tab <- left_join(rank_tab, lm_rank, by = "variable")
}
if(exists("rf_imp_df")) {
  rf_rank <- rf_imp_df %>% mutate(rf_rank = rank(-IncMSE)) %>% dplyr::select(variable, rf_rank)
  rank_tab <- left_join(rank_tab, rf_rank, by = "variable")
}
if(exists("coefs_df")) {
  lasso_rank <- coefs_df %>% mutate(lasso_nonzero = (coefficient != 0)) %>% dplyr::select(variable, lasso_nonzero)
  rank_tab <- left_join(rank_tab, lasso_rank, by = "variable")
}
# Replace NAs with large rank or FALSE
rank_tab <- rank_tab %>% mutate(lm_rank = ifelse(is.na(lm_rank), Inf, lm_rank),
                                rf_rank = ifelse(is.na(rf_rank), Inf, rf_rank),
                                lasso_nonzero = ifelse(is.na(lasso_nonzero), FALSE, lasso_nonzero))
# Compute aggregate score (lower better): sum of ranks (ignore Inf by replacing with a big number)
rank_tab <- rank_tab %>% mutate(agg_score = corr_rank + ifelse(is.infinite(lm_rank), 100, lm_rank) + ifelse(is.infinite(rf_rank), 100, rf_rank))
rank_tab <- rank_tab %>% arrange(agg_score)
write.csv(rank_tab, file.path(imp_dir, "combined_importance_ranks.csv"), row.names = FALSE)

# 6) Simple summary file with clear recommendations
summary_lines <- c(
  "Important Factors Extraction - Summary",
  paste("Data rows used for analyses:", nrow(data)),
  "Files created in folder: Important_Factors (correlation, lm coeffs, RF importance, Lasso coeffs if available).",
  "",
  "How to read results:",
  "- correlation_with_target.csv: Pearson correlation (sign indicates direction). Check top abs() values.",
  "- standardized_lm_coefficients.csv: standardized betas; magnitude shows relative effect size (controlling other predictors).",
  "- rf_variable_importance_IncMSE.csv: Random forest importance (%IncMSE) — higher means more impact on prediction.",
  "- lasso_nonzero_coefficients.csv: variables Lasso kept (non-zero) — good for feature selection.",
  "- combined_importance_ranks.csv: aggregated rank across methods (lower = more consistently important)."
)
writeLines(summary_lines, file.path(imp_dir, "how_to_interpret_results.txt"))

cat("Important factors extraction finished. Check folder:", imp_dir, "\n")
print(list.files(imp_dir, recursive = TRUE))
