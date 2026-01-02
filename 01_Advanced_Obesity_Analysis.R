###############################################################################
# Obesity Advanced Analysis
# Author: Abeer Gaber
# All select() calls use dplyr::select() to avoid conflicts
###############################################################################

# ====== 0. Paths & options ======
data_path  <- "D:/Data Analysis/Finished Projects/Obesity_Causes_Analysis/Cleaned_Final_Data.csv"
output_dir <- "D:/Data Analysis/Finished Projects/Obesity_Causes_Analysis/Cleaned_Final_DataOutputs_Advanced"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
options(scipen = 999)

# ====== 1. Load packages ======
pkgs <- c("tidyverse","car","broom","glmnet","caret","GGally","corrplot","psych",
          "ggpubr","multcomp","Matrix","ggrepel")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(new)) install.packages(new, repos = "https://cloud.r-project.org")
lapply(pkgs, library, character.only = TRUE)
library(dplyr, warn.conflicts = FALSE)  # Prevent select() conflict

# ====== 2. Read data & basic cleaning ======
data <- read.csv(data_path, stringsAsFactors = FALSE)
data$Country <- stringr::str_squish(as.character(data$Country))
data$Income_Level <- stringr::str_squish(as.character(data$Income_Level))

num_cols <- setdiff(names(data), c("Country","Income_Level"))
for(col in num_cols){
  data[[col]] <- suppressWarnings(as.numeric(as.character(data[[col]])))
}

income_levels <- c("Low income","Lower-middle income","Upper-middle income","High income")
present <- unique(na.omit(data$Income_Level))
if(all(present %in% income_levels)) {
  data$Income_Level <- factor(data$Income_Level, levels = income_levels, ordered = TRUE)
} else data$Income_Level <- factor(data$Income_Level)

target <- "Avg_Adult_Obesity_BMI30"

# ====== 3. Stratified summaries & boxplots by Income_Level ======
if("Income_Level" %in% names(data)){
  strat_summary <- data %>%
    group_by(Income_Level) %>%
    summarise(across(where(is.numeric), list(mean = ~mean(., na.rm=TRUE),
                                             sd = ~sd(., na.rm=TRUE),
                                             n = ~sum(!is.na(.))),
                     .names = "{col}_{fn}"))
  write.csv(strat_summary, file.path(output_dir, "stratified_summary_by_income.csv"), row.names = FALSE)
}

if(all(c("Income_Level", target) %in% names(data))){
  p_box <- ggplot(data, aes(x = Income_Level, y = .data[[target]], fill = Income_Level)) +
    geom_boxplot(alpha = 0.9) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
    theme_minimal(base_size = 13) + theme(legend.position = "none") +
    labs(title = "Adult Obesity by Income Level", x = "Income Level", y = "Adult Obesity (BMI ≥30)")
  ggsave(file.path(output_dir, "boxplot_adult_obesity_by_income.png"), p_box, width = 8, height = 6)
}

# ====== 4. ANOVA + post-hoc ======
if("Income_Level" %in% names(data) && nlevels(data$Income_Level) > 1){
  anova_df <- data %>% dplyr::select(all_of(c(target, "Income_Level"))) %>% drop_na()
  anova_fit <- aov(as.formula(paste(target, "~ Income_Level")), data = anova_df)
  sink(file.path(output_dir, "anova_income_on_adult_obesity.txt"))
  cat("### ANOVA: Adult Obesity ~ Income_Level/n/n")
  print(summary(anova_fit))
  sink()
  
  tuk <- TukeyHSD(anova_fit)
  capture.output(tuk, file = file.path(output_dir, "tukeyHSD_income_adult_obesity.txt"))
  
  library(multcomp)
  mc <- glht(anova_fit, linfct = mcp(Income_Level = "Tukey"))
  cld <- cld(mc)
  sink(file.path(output_dir, "multcomp_cld_income.txt"))
  print(cld)
  sink()
}

# ====== 5. Multicollinearity check (VIF) ======
candidates <- c("Avg_SugarCons_PerCapita_2019_2022","Annual_FastFood_Sales_USD",
                "Animal_Protein_Intake","Avg_Fruit","Avg_Vege_cons",
                "Avg_Physical_Act2000_2030","Avg_Age")
candidates <- candidates[candidates %in% names(data)]

vif_df <- data %>% dplyr::select(all_of(c(target, candidates))) %>% drop_na()
if(nrow(vif_df) > 10 && length(candidates) >= 2){
  lm_full <- lm(as.formula(paste(target, "~", paste(candidates, collapse = " + "))), data = vif_df)
  vif_vals <- car::vif(lm_full)
  vif_tab <- data.frame(variable = names(vif_vals), VIF = as.numeric(vif_vals))
  write.csv(vif_tab, file.path(output_dir, "vif_values.csv"), row.names = FALSE)
} else {
  message("لا توجد بيانات كافية أو متغيرات قليلة جداً للفحص VIF.")
}

# ====== 6. Stepwise model selection (AIC) ======
step_vars <- unique(c(target, candidates, "Income_Level"))
step_df <- data %>% dplyr::select(any_of(step_vars)) %>% drop_na()
if(nrow(step_df) > 20){
  null_mod <- lm(as.formula(paste(target, "~ 1")), data = step_df)
  full_mod <- lm(as.formula(paste(target, "~", paste(setdiff(names(step_df), target), collapse = " + "))), data = step_df)
  step_mod <- step(null_mod, scope = list(upper = full_mod), direction = "both", trace = 0)
  sink(file.path(output_dir, "stepwise_selection_summary.txt"))
  cat("### Stepwise (AIC) selected model:/n")
  print(summary(step_mod))
  sink()
}

# ====== 7. Interaction effects ======
inter_vars <- c("Avg_SugarCons_PerCapita_2019_2022","Annual_FastFood_Sales_USD")
inter_vars <- inter_vars[inter_vars %in% names(data)]
if(length(inter_vars) == 2 && target %in% names(data)){
  inter_df <- data %>% dplyr::select(all_of(c(target, inter_vars, "Income_Level"))) %>% drop_na()
  formula_inter <- as.formula(paste(target, "~", inter_vars[1], "*", inter_vars[2], "+", "Income_Level"))
  inter_fit <- try(lm(formula_inter, data = inter_df), silent = TRUE)
  if(!inherits(inter_fit, "try-error")){
    sink(file.path(output_dir, "interaction_model_summary.txt"))
    cat("### Interaction model: sugar * fastfood + Income_Level/n")
    print(summary(inter_fit))
    sink()
    p_int <- ggplot(inter_df, aes_string(x = inter_vars[1], y = target, color = "Income_Level")) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", aes_string(group = "Income_Level"), se = FALSE) +
      labs(title = paste("Interaction:", inter_vars[1], "x", inter_vars[2], "by Income_Level"),
           x = inter_vars[1], y = target) +
      theme_minimal(base_size = 12)
    ggsave(file.path(output_dir, "interaction_sugar_income_plot.png"), p_int, width = 9, height = 6)
  }
}

# ====== 8. Regularized regression (Lasso & Ridge) ======
glmnet_vars <- c("Avg_SugarCons_PerCapita_2019_2022","Annual_FastFood_Sales_USD",
                 "Animal_Protein_Intake","Avg_Fruit","Avg_Vege_cons",
                 "Avg_Physical_Act2000_2030","Avg_Age")
glmnet_vars <- glmnet_vars[glmnet_vars %in% names(data)]
glmnet_df <- data %>% dplyr::select(all_of(c(target, glmnet_vars))) %>% drop_na()
if(nrow(glmnet_df) >= 30 && length(glmnet_vars) >= 2){
  x <- model.matrix(as.formula(paste("~", paste(glmnet_vars, collapse = "+"))), data = glmnet_df)[,-1]
  y <- glmnet_df[[target]]
  set.seed(123)
  cv_ridge <- cv.glmnet(x, y, alpha = 0, nfolds = 5)
  best_lambda_ridge <- cv_ridge$lambda.min
  fit_ridge <- glmnet(x, y, alpha = 0, lambda = best_lambda_ridge)
  cv_lasso <- cv.glmnet(x, y, alpha = 1, nfolds = 5)
  best_lambda_lasso <- cv_lasso$lambda.min
  fit_lasso <- glmnet(x, y, alpha = 1, lambda = best_lambda_lasso)
  write.csv(as.matrix(coef(fit_ridge)), file.path(output_dir, "glmnet_ridge_coefficients.csv"))
  write.csv(as.matrix(coef(fit_lasso)), file.path(output_dir, "glmnet_lasso_coefficients.csv"))
  png(file.path(output_dir, "cv_ridge_plot.png"), width = 800, height = 600); plot(cv_ridge); dev.off()
  png(file.path(output_dir, "cv_lasso_plot.png"), width = 800, height = 600); plot(cv_lasso); dev.off()
  set.seed(123)
  train_idx <- createDataPartition(y, p = 0.8, list = FALSE)
  x_train <- x[train_idx,]; y_train <- y[train_idx]
  x_test  <- x[-train_idx,]; y_test  <- y[-train_idx]
  pred_ridge <- predict(fit_ridge, newx = x_test)
  pred_lasso <- predict(fit_lasso, newx = x_test)
  perf <- data.frame(
    Model = c("Ridge","Lasso"),
    RMSE = c(sqrt(mean((y_test - pred_ridge)^2)), sqrt(mean((y_test - pred_lasso)^2))),
    MAE  = c(mean(abs(y_test - pred_ridge)), mean(abs(y_test - pred_lasso)))
  )
  write.csv(perf, file.path(output_dir, "glmnet_model_performance.csv"), row.names = FALSE)
}

# ====== 9. Diagnostics ======
model_for_diag <- if(exists("step_mod") && inherits(step_mod, "lm")) step_mod else if(exists("lm_full")) lm_full else NULL
if(!is.null(model_for_diag)){
  png(file.path(output_dir, "residuals_plot.png"), width = 900, height = 600)
  par(mfrow=c(2,2))
  plot(model_for_diag)
  dev.off()
  cd <- cooks.distance(model_for_diag)
  inf_df <- data.frame(Country = data$Country[as.numeric(names(cd))], cooksD = cd)
  inf_df <- inf_df %>% arrange(desc(cooksD)) %>% slice_head(n = 10)
  write.csv(inf_df, file.path(output_dir, "top_influential_points_cooksD.csv"), row.names = FALSE)
}

# ====== 10. Summary & Interpretation Hints ======
summary_lines <- c(
  paste("Data rows:", nrow(data), "Columns:", ncol(data)),
  if(exists("vif_tab")) "VIF computed.",
  if(exists("step_mod")) "Stepwise model computed.",
  if(exists("cv_lasso")) paste("Lasso λ.min:", best_lambda_lasso),
  if(exists("cv_ridge")) paste("Ridge λ.min:", best_lambda_ridge)
)
writeLines(summary_lines, file.path(output_dir, "advanced_analysis_summary.txt"))

interpret <- c(
  "Interpretation hints:",
  "- Check VIF.csv: values >5-10 indicate multicollinearity.",
  "- Significant ANOVA results (p<0.05) indicate differences between income levels.",
  "- Lasso coefficients = 0 indicate the variable is excluded (not influential).",
  "- If the sugar*fastfood interaction is significant, the relationship between sugar and obesity depends on fast food consumption.",
  "- Check Cook's D to identify highly influential countries.".
)
writeLines(interpret, file.path(output_dir, "advanced_interpretation_hints.txt"))

cat("✅ Advanced analysis completed successfully. Check your output folder:/n", output_dir, "/n")
