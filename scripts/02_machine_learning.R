###########################################
#  - Bachelor Project                     #
#  - miRNA biomarkers in Colon Cancer     #
#  - Caret RUS + 10-Fold CV               #
#  - Copyright: @Malak Wael               #
###########################################

library(caret)
library(randomForest)
library(e1071)
library(glmnet)
library(nnet)
library(pROC)
library(PRROC)
library(ggplot2)

# ============================================================
# STEP 1 — Load data (outputs of 01_preprocessing_differential_expression.R)
# ============================================================

dems.exp.raw <- read.csv(
  "DEMs_normalized_expression.csv",
  row.names = 1, check.names = FALSE
)

pheno.ml <- read.csv(
  "pheno_clean_table.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)

common.samples <- intersect(colnames(dems.exp.raw), pheno.ml$Sample_ID)
dems.exp.raw   <- dems.exp.raw[, common.samples]
pheno.ml       <- pheno.ml[match(common.samples, pheno.ml$Sample_ID), ]

exp_log    <- log2(dems.exp.raw + 1)
X          <- as.data.frame(t(exp_log))
X$label    <- factor(pheno.ml$Tissue_Type, levels = c("Normal", "Tumor"))
X_features <- X[, !colnames(X) %in% "label"]
y          <- X$label

cat("Dimensions:", nrow(X), "samples x", ncol(X_features), "miRNAs\n")
print(table(y))

# ============================================================
# STEP 2 — Helper functions
# ============================================================

evaluate_caret <- function(model, model_name) {
  preds <- model$pred

  cm  <- confusionMatrix(preds$pred, preds$obs, positive = "Normal")
  roc <- roc(preds$obs, preds$Normal,
             levels = c("Tumor", "Normal"), quiet = TRUE)

  labels_binary <- ifelse(preds$obs == "Normal", 1, 0)
  pr <- pr.curve(
    scores.class0 = preds$Normal[labels_binary == 1],
    scores.class1 = preds$Normal[labels_binary == 0],
    curve = TRUE
  )
  PR_AUC <- pr$auc.integral

  TP  <- as.numeric(cm$table["Normal", "Normal"])
  TN  <- as.numeric(cm$table["Tumor",  "Tumor"])
  FP  <- as.numeric(cm$table["Normal", "Tumor"])
  FN  <- as.numeric(cm$table["Tumor",  "Normal"])
  MCC <- (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

  cat("\n========================================\n")
  cat("Model:", model_name, "\n")
  cat("========================================\n")
  cat("Accuracy:   ", round(cm$overall["Accuracy"]*100, 2), "%\n")
  cat("ROC-AUC:    ", round(auc(roc)*100, 2), "%\n")
  cat("PR-AUC:     ", round(PR_AUC*100, 2), "%\n")
  cat("Sensitivity:", round(cm$byClass["Sensitivity"]*100, 2), "%\n")
  cat("Specificity:", round(cm$byClass["Specificity"]*100, 2), "%\n")
  cat("F1:         ", round(cm$byClass["F1"]*100, 2), "%\n")
  cat("MCC:        ", round(MCC, 4), "\n")
  print(cm$table)

  return(data.frame(
    Model       = model_name,
    Accuracy    = round(cm$overall["Accuracy"]*100, 2),
    ROC_AUC     = round(auc(roc)*100, 2),
    PR_AUC      = round(PR_AUC*100, 2),
    Sensitivity = round(cm$byClass["Sensitivity"]*100, 2),
    Specificity = round(cm$byClass["Specificity"]*100, 2),
    F1          = round(cm$byClass["F1"]*100, 2),
    MCC         = round(MCC, 4)
  ))
}

normalize_0_100 <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  if (max(x) == min(x)) return(rep(0, length(x)))
  return((x - min(x)) / (max(x) - min(x)) * 100)
}

plot_confusion_matrix <- function(cm, model_name, fill_color) {
  cm_df <- as.data.frame(cm$table)
  colnames(cm_df) <- c("Predicted", "Actual", "Count")

  cm_df$RowTotal <- ave(cm_df$Count, cm_df$Predicted, FUN = sum)
  cm_df$Pct <- round(cm_df$Count / cm_df$RowTotal * 100, 1)
  cm_df$Label <- paste0(cm_df$Count, "\n(", cm_df$Pct, "%)")

  p <- ggplot(cm_df, aes(x = Actual, y = Predicted, fill = Pct)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = Label), size = 5, fontface = "bold") +
    scale_fill_gradient(low = "white", high = fill_color, limits = c(0, 100)) +
    labs(title = paste("Confusion Matrix —", model_name),
         x = "Actual", y = "Predicted") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(paste0("confusion_matrix_", model_name, ".png"),
         p, width = 6, height = 5, dpi = 300)
}

plot_top15 <- function(df, score_col, title, filename, fill_color) {
  df_top <- head(df[order(-df[[score_col]]), ], 15)
  df_top <- df_top[order(df_top[[score_col]]), ]
  df_top$miRNA <- factor(df_top$miRNA, levels = df_top$miRNA)

  p <- ggplot(df_top, aes(x = miRNA, y = .data[[score_col]])) +
    geom_bar(stat = "identity", fill = fill_color) +
    coord_flip() +
    labs(title = title, x = NULL, y = "Score") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 9))

  ggsave(filename, p, width = 9, height = 7, dpi = 300)
}

# ============================================================
# STEP 3 — CV Control
# ============================================================

ctrl_10fold <- trainControl(
  method          = "cv",
  number          = 10,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  sampling        = "down"
)

# ============================================================
# STEP 4 — Train all 4 models
# ============================================================

set.seed(42)
rf_model <- caret::train(
  x = X_features, y = y, method = "rf",
  trControl = ctrl_10fold, metric = "ROC",
  ntree = 300, importance = TRUE
)
rf_results <- evaluate_caret(rf_model, "RF")

set.seed(42)
svm_model <- caret::train(
  x = X_features, y = y, method = "svmRadial",
  trControl = ctrl_10fold, metric = "ROC",
  preProcess = c("center", "scale")
)
svm_results <- evaluate_caret(svm_model, "SVM")

set.seed(42)
lasso_model <- caret::train(
  x = X_features, y = y, method = "glmnet",
  trControl = ctrl_10fold, metric = "ROC",
  preProcess = c("center", "scale"),
  tuneGrid = expand.grid(alpha  = 1,
                         lambda = seq(0.001, 0.1, length = 20))
)
lasso_results <- evaluate_caret(lasso_model, "LASSO")

set.seed(42)
nn_model <- caret::train(
  x = X_features, y = y, method = "nnet",
  trControl = ctrl_10fold, metric = "ROC",
  preProcess = c("center", "scale"),
  tuneGrid = expand.grid(size = 10, decay = 0.001),
  maxit = 200, trace = FALSE, MaxNWts = 5000
)
nn_results <- evaluate_caret(nn_model, "NN")

# ============================================================
# STEP 5 — Summary Table
# ============================================================

summary_table <- rbind(rf_results, svm_results, lasso_results, nn_results)
rownames(summary_table) <- NULL

cat("\n========================================\n")
cat("FINAL SUMMARY — Caret Pipeline\n")
cat("========================================\n")
print(summary_table)
write.csv(summary_table, "metrics_caret.csv", row.names = FALSE)

# ============================================================
# STEP 6 — Confusion matrix plots
# ============================================================

cm_rf    <- confusionMatrix(rf_model$pred$pred,    rf_model$pred$obs,    positive = "Normal")
cm_svm   <- confusionMatrix(svm_model$pred$pred,   svm_model$pred$obs,   positive = "Normal")
cm_lasso <- confusionMatrix(lasso_model$pred$pred, lasso_model$pred$obs, positive = "Normal")
cm_nn    <- confusionMatrix(nn_model$pred$pred,    nn_model$pred$obs,    positive = "Normal")

plot_confusion_matrix(cm_rf,    "RF",             "#00BFFF")
plot_confusion_matrix(cm_svm,   "SVM",            "#FF1493")
plot_confusion_matrix(cm_lasso, "LASSO",          "#32CD32")
plot_confusion_matrix(cm_nn,    "Neural_Network", "#A020F0")

# ============================================================
# STEP 7 — Combined ROC curve
# ============================================================

roc_rf    <- roc(rf_model$pred$obs,    rf_model$pred$Normal,    levels = c("Tumor", "Normal"), quiet = TRUE)
roc_svm   <- roc(svm_model$pred$obs,   svm_model$pred$Normal,   levels = c("Tumor", "Normal"), quiet = TRUE)
roc_lasso <- roc(lasso_model$pred$obs, lasso_model$pred$Normal, levels = c("Tumor", "Normal"), quiet = TRUE)
roc_nn    <- roc(nn_model$pred$obs,    nn_model$pred$Normal,    levels = c("Tumor", "Normal"), quiet = TRUE)

roc_df <- rbind(
  data.frame(
    FPR = 1 - roc_rf$specificities,
    TPR = roc_rf$sensitivities,
    Model = paste0("RF AUC=", round(auc(roc_rf) * 100, 2), "%")
  ),
  data.frame(
    FPR = 1 - roc_svm$specificities,
    TPR = roc_svm$sensitivities,
    Model = paste0("SVM AUC=", round(auc(roc_svm) * 100, 2), "%")
  ),
  data.frame(
    FPR = 1 - roc_lasso$specificities,
    TPR = roc_lasso$sensitivities,
    Model = paste0("LASSO AUC=", round(auc(roc_lasso) * 100, 2), "%")
  ),
  data.frame(
    FPR = 1 - roc_nn$specificities,
    TPR = roc_nn$sensitivities,
    Model = paste0("NN AUC=", round(auc(roc_nn) * 100, 2), "%")
  )
)

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 1.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("#00BFFF", "#FF1493", "#32CD32", "#A020F0")) +
  labs(title = "Combined ROC Curves — Four Models",
       x = "False Positive Rate",
       y = "True Positive Rate",
       color = "Model") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

ggsave("combined_ROC_4models.png", p_roc, width = 8, height = 6, dpi = 300)

# ============================================================
# STEP 8 — Biomarker importance from the trained caret models
# ============================================================

# RF: MeanDecreaseAccuracy
rf_imp <- importance(rf_model$finalModel)
rf_df  <- data.frame(miRNA = rownames(rf_imp),
                     RF = normalize_0_100(rf_imp[, "MeanDecreaseAccuracy"]))

# SVM: caret varImp
svm_imp <- varImp(svm_model)$importance
svm_df  <- data.frame(miRNA = rownames(svm_imp),
                      SVM = normalize_0_100(svm_imp$Normal))

# LASSO: absolute coefficients at the selected lambda
lasso_coef <- coef(lasso_model$finalModel, s = lasso_model$bestTune$lambda)
lasso_df   <- data.frame(miRNA = rownames(lasso_coef),
                         coef = as.numeric(lasso_coef))
lasso_df   <- lasso_df[lasso_df$miRNA != "(Intercept)", ]
lasso_df$LASSO <- normalize_0_100(abs(lasso_df$coef))

cat("\nmiRNAs kept by LASSO (non-zero coefficients):\n")
print(lasso_df[lasso_df$coef != 0, c("miRNA", "coef")])

# NN: caret varImp
nn_imp <- varImp(nn_model)$importance
nn_df  <- data.frame(miRNA = rownames(nn_imp),
                     NN = normalize_0_100(nn_imp[, 1]))

# ============================================================
# STEP 9 — Top 15 biomarker plots per model
# ============================================================

write.csv(head(rf_df[order(-rf_df$RF), ], 15),          "biomarkers_RF_top15.csv",    row.names = FALSE)
write.csv(head(svm_df[order(-svm_df$SVM), ], 15),       "biomarkers_SVM_top15.csv",   row.names = FALSE)
write.csv(head(lasso_df[order(-lasso_df$LASSO), ], 15), "biomarkers_LASSO_top15.csv", row.names = FALSE)
write.csv(head(nn_df[order(-nn_df$NN), ], 15),          "biomarkers_NN_top15.csv",    row.names = FALSE)

plot_top15(rf_df, "RF",
           "Random Forest — Top 15 miRNA Biomarkers",
           "barplot_RF_top15.png", "#00BFFF")

plot_top15(svm_df, "SVM",
           "SVM — Top 15 miRNA Biomarkers",
           "barplot_SVM_top15.png", "#FF1493")

plot_top15(lasso_df, "LASSO",
           "LASSO — Top 15 miRNA Biomarkers",
           "barplot_LASSO_top15.png", "#32CD32")

plot_top15(nn_df, "NN",
           "Neural Network — Top 15 miRNA Biomarkers",
           "barplot_NN_top15.png", "#A020F0")

# ============================================================
# STEP 10 — Combined score across all DEMs
# ============================================================

combined_scores <- Reduce(function(a, b) merge(a, b, by = "miRNA", all = TRUE),
                          list(rf_df,
                               svm_df,
                               lasso_df[, c("miRNA", "LASSO")],
                               nn_df))

combined_scores[is.na(combined_scores)] <- 0
combined_scores$Combined_Score <- rowMeans(combined_scores[, c("RF", "SVM", "LASSO", "NN")])
combined_scores <- combined_scores[order(-combined_scores$Combined_Score), ]
combined_scores$Rank <- 1:nrow(combined_scores)

top10_combined <- head(combined_scores, 10)

write.csv(combined_scores, "all_DEMs_combined_scores_4models.csv", row.names = FALSE)
write.csv(top10_combined,  "top10_combined_biomarkers_4models.csv", row.names = FALSE)

cat("\n========================================\n")
cat("TOP 10 COMBINED BIOMARKERS — 4 MODELS\n")
cat("========================================\n")
print(top10_combined)

# Combined top 10 plot
top10_plot <- top10_combined[order(top10_combined$Combined_Score), ]
top10_plot$miRNA <- factor(top10_plot$miRNA, levels = top10_plot$miRNA)

p_combined <- ggplot(top10_plot, aes(x = miRNA, y = Combined_Score)) +
  geom_bar(stat = "identity", fill = "#FF8C00") +
  coord_flip() +
  labs(title = "Top 10 Combined miRNA Biomarkers",
       subtitle = "Mean normalized score across RF, SVM, LASSO, and Neural Network",
       x = NULL,
       y = "Combined Score (0–100)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 9))

ggsave("barplot_top10_combined_4models.png", p_combined,
       width = 9, height = 7, dpi = 300)

# ============================================================
# STEP 11 — Bootstrap (RUS, 100 iterations, out-of-bag evaluation)
# ============================================================

n.samples <- nrow(X)

run.bootstrap.rus <- function(model.type, n.boot = 100) {

  cat("\nBootstrap:", model.type, "— RUS (100 iterations)\n")

  boot.acc  <- c(); boot.sens <- c()
  boot.spec <- c(); boot.f1   <- c()
  boot.auc  <- c(); boot.mcc  <- c()

  set.seed(42)

  for(i in 1:n.boot) {

    # Bootstrap sample with replacement
    boot.idx   <- sample(1:n.samples, n.samples, replace = TRUE)
    oob.idx    <- setdiff(1:n.samples, boot.idx)
    oob.labels <- y[oob.idx]

    if(sum(oob.labels == "Normal") == 0) next

    boot.feat <- X_features[boot.idx, ]
    boot.lab  <- y[boot.idx]
    oob.feat  <- X_features[oob.idx, ]

    # RUS — undersample Tumor to match Normal count
    n.norm <- sum(boot.lab == "Normal")
    n.tum  <- sum(boot.lab == "Tumor")
    if(n.norm == 0) next

    normal.idx <- which(boot.lab == "Normal")
    tumor.idx  <- which(boot.lab == "Tumor")
    tumor.keep <- sample(tumor.idx, min(n.norm, length(tumor.idx)))
    keep.idx   <- c(normal.idx, tumor.keep)

    boot.feat.rus <- boot.feat[keep.idx, ]
    boot.lab.rus  <- boot.lab[keep.idx]

    # Scale for SVM and NN
    feat.means <- colMeans(boot.feat.rus)
    feat.sds   <- apply(boot.feat.rus, 2, sd)
    feat.sds[feat.sds == 0] <- 1

    boot.feat.scaled <- as.data.frame(
      scale(boot.feat.rus,
            center = feat.means,
            scale  = feat.sds))
    oob.feat.scaled  <- as.data.frame(
      scale(oob.feat,
            center = feat.means,
            scale  = feat.sds))

    model <- tryCatch({

      if(model.type == "RF") {
        randomForest(x     = boot.feat.rus,
                     y     = boot.lab.rus,
                     ntree = 300,
                     mtry  = floor(sqrt(ncol(boot.feat.rus))))

      } else if(model.type == "SVM") {
        svm(x           = boot.feat.scaled,
            y           = boot.lab.rus,
            kernel      = "radial",
            cost        = 1,
            gamma       = 1/ncol(boot.feat.scaled),
            probability = TRUE)

      } else if(model.type == "LASSO") {
        bx     <- as.matrix(boot.feat.scaled)
        by     <- ifelse(boot.lab.rus == "Normal", 1, 0)
        cv.fit <- cv.glmnet(bx, by,
                            family       = "binomial",
                            alpha        = 1,
                            nfolds       = 5,
                            type.measure = "auc")
        glmnet(bx, by,
               family = "binomial",
               alpha  = 1,
               lambda = cv.fit$lambda.min)

      } else if(model.type == "NN") {
        nnet(x       = boot.feat.scaled,
             y       = class.ind(boot.lab.rus),
             size    = 10,
             decay   = 0.001,
             maxit   = 200,
             softmax = TRUE,
             trace   = FALSE,
             MaxNWts = 5000)
      }

    }, error = function(e) NULL)

    if(is.null(model)) next

    # Predict
    if(model.type == "RF") {
      pred.class <- predict(model, oob.feat)
      pred.prob  <- predict(model, oob.feat,
                            type = "prob")[, "Normal"]

    } else if(model.type == "SVM") {
      pred.class <- predict(model, oob.feat.scaled)
      pred.prob  <- attr(
        predict(model, oob.feat.scaled, probability = TRUE),
        "probabilities")[, "Normal"]

    } else if(model.type == "LASSO") {
      pred.prob  <- predict(model,
                            newx = as.matrix(oob.feat.scaled),
                            type = "response")[,1]
      pred.class <- factor(ifelse(pred.prob > 0.5,
                                  "Normal", "Tumor"),
                           levels = c("Normal", "Tumor"))

    } else if(model.type == "NN") {
      pred.raw   <- predict(model, oob.feat.scaled)
      pred.prob  <- pred.raw[, "Normal"]
      pred.class <- factor(
        colnames(pred.raw)[apply(pred.raw, 1, which.max)],
        levels = c("Normal", "Tumor"))
    }

    # Metrics
    cm <- confusionMatrix(pred.class, oob.labels,
                          positive = "Normal")

    tp    <- cm$table[1,1]; tn <- cm$table[2,2]
    fp    <- cm$table[1,2]; fn <- cm$table[2,1]
    denom <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    mcc   <- ifelse(denom == 0, NA,
                    (tp*tn - fp*fn) / denom)

    auc.val <- tryCatch(
      as.numeric(auc(roc(
        as.numeric(oob.labels == "Normal"),
        pred.prob, quiet = TRUE))),
      error = function(e) NA)

    boot.acc  <- c(boot.acc,  cm$overall["Accuracy"])
    boot.sens <- c(boot.sens, cm$byClass["Sensitivity"])
    boot.spec <- c(boot.spec, cm$byClass["Specificity"])
    boot.f1   <- c(boot.f1,   cm$byClass["F1"])
    boot.auc  <- c(boot.auc,  auc.val)
    boot.mcc  <- c(boot.mcc,  mcc)

    if(i %% 20 == 0) cat("  Iteration:", i, "\n")
  }

  data.frame(
    Model       = model.type,
    Accuracy    = round(mean(boot.acc,  na.rm=TRUE)*100, 2),
    Sensitivity = round(mean(boot.sens, na.rm=TRUE)*100, 2),
    Specificity = round(mean(boot.spec, na.rm=TRUE)*100, 2),
    F1          = round(mean(boot.f1,   na.rm=TRUE)*100, 2),
    AUC         = round(mean(boot.auc,  na.rm=TRUE)*100, 2),
    MCC         = round(mean(boot.mcc,  na.rm=TRUE), 4),
    Acc_SD      = round(sd(boot.acc,    na.rm=TRUE)*100, 2),
    Sens_SD     = round(sd(boot.sens,   na.rm=TRUE)*100, 2),
    Spec_SD     = round(sd(boot.spec,   na.rm=TRUE)*100, 2),
    F1_SD       = round(sd(boot.f1,     na.rm=TRUE)*100, 2),
    AUC_SD      = round(sd(boot.auc,    na.rm=TRUE)*100, 2),
    MCC_SD      = round(sd(boot.mcc,    na.rm=TRUE), 4)
  )
}

boot_results <- data.frame()
boot_results <- rbind(boot_results, run.bootstrap.rus("RF",    100))
boot_results <- rbind(boot_results, run.bootstrap.rus("SVM",   100))
boot_results <- rbind(boot_results, run.bootstrap.rus("LASSO", 100))
boot_results <- rbind(boot_results, run.bootstrap.rus("NN",    100))

cat("\n===========================================\n")
cat("BOOTSTRAP RESULTS — RF, SVM, LASSO, NN\n")
cat("RUS | Mean ± SD | 100 Iterations\n")
cat("===========================================\n")
cat(sprintf("%-8s %-12s %-14s %-14s %-12s %-12s %-10s\n",
            "Model", "Accuracy",
            "Sensitivity", "Specificity",
            "F1", "AUC", "MCC"))
cat(sprintf("%-8s %-12s %-14s %-14s %-12s %-12s %-10s\n",
            "-----", "--------",
            "-----------", "-----------",
            "--", "---", "---"))

for(r in 1:nrow(boot_results)) {
  cat(sprintf("%-8s %-12s %-14s %-14s %-12s %-12s %-10s\n",
              boot_results$Model[r],
              paste0(boot_results$Accuracy[r],    "±", boot_results$Acc_SD[r]),
              paste0(boot_results$Sensitivity[r], "±", boot_results$Sens_SD[r]),
              paste0(boot_results$Specificity[r], "±", boot_results$Spec_SD[r]),
              paste0(boot_results$F1[r],          "±", boot_results$F1_SD[r]),
              paste0(boot_results$AUC[r],         "±", boot_results$AUC_SD[r]),
              paste0(boot_results$MCC[r],         "±", boot_results$MCC_SD[r])))
}
cat("===========================================\n")

write.csv(boot_results, "Bootstrap_RUS_4Models.csv", row.names = FALSE)
cat("\nSaved: Bootstrap_RUS_4Models.csv\n")

# ============================================================
# STEP 12 — Permutation test (labels shuffled, 20 iterations)
# Re-implementation: the original permutation code was lost, so the
# AUC ranges will not match the thesis exactly.
# ============================================================

run.permutation <- function(method, n.perm = 20, ...) {

  cat("\nPermutation test:", method, "(", n.perm, "iterations )\n")
  perm.auc <- c()

  for (i in 1:n.perm) {
    set.seed(100 + i)
    y.perm <- sample(y)                     # shuffle the labels

    set.seed(42)
    fit <- caret::train(
      x = X_features, y = y.perm, method = method,
      trControl = ctrl_10fold, metric = "ROC", ...
    )

    # direction fixed so a worse-than-random model gets AUC < 50%
    roc.perm <- roc(fit$pred$obs, fit$pred$Normal,
                    levels = c("Tumor", "Normal"),
                    direction = "<", quiet = TRUE)
    perm.auc <- c(perm.auc, as.numeric(auc(roc.perm)) * 100)

    if (i %% 5 == 0) cat("  Iteration:", i, "\n")
  }

  return(perm.auc)
}

perm_rf <- run.permutation("rf", ntree = 300)

perm_svm <- run.permutation("svmRadial",
                            preProcess = c("center", "scale"))

perm_lasso <- run.permutation("glmnet",
                              preProcess = c("center", "scale"),
                              tuneGrid = expand.grid(alpha  = 1,
                                                     lambda = seq(0.001, 0.1, length = 20)))

perm_nn <- run.permutation("nnet",
                           preProcess = c("center", "scale"),
                           tuneGrid = expand.grid(size = 10, decay = 0.001),
                           maxit = 200, trace = FALSE, MaxNWts = 5000)

real_auc <- c(rf_results$ROC_AUC, svm_results$ROC_AUC,
              lasso_results$ROC_AUC, nn_results$ROC_AUC)
perm_list <- list(perm_rf, perm_svm, perm_lasso, perm_nn)

perm_table <- data.frame(
  Model     = c("RF", "SVM", "LASSO", "NN"),
  Real_AUC  = real_auc,
  Perm_Min  = round(sapply(perm_list, min), 2),
  Perm_Max  = round(sapply(perm_list, max), 2),
  Perm_Mean = round(sapply(perm_list, mean), 2),
  # share of permuted runs with AUC >= the real AUC
  p_value   = sapply(1:4, function(k) mean(perm_list[[k]] >= real_auc[k]))
)

cat("\n========================================\n")
cat("PERMUTATION TEST — 20 iterations\n")
cat("========================================\n")
print(perm_table)
write.csv(perm_table, "Permutation_test_4Models.csv", row.names = FALSE)

# ============================================================
# STEP 12 — Permutation test (20 iterations)
# Re-implementation: the original permutation code was lost, so the
# permuted AUC ranges will differ slightly from those in the thesis.
# Labels are shuffled and the same caret 10-fold CV + RUS setup is rerun.
# ============================================================

run.permutation <- function(model.type, n.perm = 20) {

  cat("\nPermutation test:", model.type, "(", n.perm, "iterations)\n")

  perm.auc <- c()

  for(i in 1:n.perm) {

    set.seed(42 + i)
    y.perm <- sample(y)   # shuffle labels, same class counts

    set.seed(42)
    model <- tryCatch({

      if(model.type == "RF") {
        caret::train(x = X_features, y = y.perm, method = "rf",
                     trControl = ctrl_10fold, metric = "ROC",
                     ntree = 300)

      } else if(model.type == "SVM") {
        caret::train(x = X_features, y = y.perm, method = "svmRadial",
                     trControl = ctrl_10fold, metric = "ROC",
                     preProcess = c("center", "scale"))

      } else if(model.type == "LASSO") {
        caret::train(x = X_features, y = y.perm, method = "glmnet",
                     trControl = ctrl_10fold, metric = "ROC",
                     preProcess = c("center", "scale"),
                     tuneGrid = expand.grid(alpha  = 1,
                                            lambda = seq(0.001, 0.1, length = 20)))

      } else if(model.type == "NN") {
        caret::train(x = X_features, y = y.perm, method = "nnet",
                     trControl = ctrl_10fold, metric = "ROC",
                     preProcess = c("center", "scale"),
                     tuneGrid = expand.grid(size = 10, decay = 0.001),
                     maxit = 200, trace = FALSE, MaxNWts = 5000)
      }

    }, error = function(e) NULL)

    if(is.null(model)) next

    preds   <- model$pred
    auc.val <- as.numeric(auc(roc(preds$obs, preds$Normal,
                                  levels = c("Tumor", "Normal"), quiet = TRUE)))
    perm.auc <- c(perm.auc, auc.val * 100)

    if(i %% 5 == 0) cat("  Iteration:", i, "\n")
  }

  return(perm.auc)
}

real.auc <- c(RF    = rf_results$ROC_AUC,
              SVM   = svm_results$ROC_AUC,
              LASSO = lasso_results$ROC_AUC,
              NN    = nn_results$ROC_AUC)

perm_results <- data.frame()

for(m in c("RF", "SVM", "LASSO", "NN")) {
  perm.auc <- run.permutation(m, 20)
  # p-value with the +1 correction, so it can never be exactly 0
  p.val <- (sum(perm.auc >= real.auc[m]) + 1) / (length(perm.auc) + 1)

  perm_results <- rbind(perm_results, data.frame(
    Model         = m,
    Real_AUC      = real.auc[m],
    Perm_AUC_min  = round(min(perm.auc), 2),
    Perm_AUC_max  = round(max(perm.auc), 2),
    Perm_AUC_mean = round(mean(perm.auc), 2),
    n_perm        = length(perm.auc),
    p_value       = round(p.val, 4)
  ))
}
rownames(perm_results) <- NULL

cat("\n===========================================\n")
cat("PERMUTATION TEST — 20 iterations\n")
cat("===========================================\n")
print(perm_results)

write.csv(perm_results, "Permutation_test_4Models.csv", row.names = FALSE)
cat("\nSaved: Permutation_test_4Models.csv\n")

cat("\n=== PIPELINE COMPLETE ===\n")
