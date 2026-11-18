# Machine Learning-Guided Discovery of miRNA Biomarkers in Colon Cancer

Bachelor's thesis project, German University in Cairo (2026)
**Author:** Malak Wael Talaat · **Supervisor:** Dr. Mohamed Hamed

## Overview

This project looks for microRNAs (miRNAs) that distinguish colon tumour tissue from adjacent normal tissue, as candidate biomarkers that could support simpler, less invasive colon cancer testing. It combines differential expression analysis with four machine learning classifiers and ranks miRNAs by how consistently the models rely on them.

**Data:** TCGA-COAD miRNA-seq, 1,881 miRNAs across 458 samples (450 tumour, 8 adjacent normal).

## Pipeline

| Step | Script | What it does |
|---|---|---|
| 1 | `01_preprocessing_differential_expression.R` | Loads raw GDC counts, QC plots, filtering (1,881 → 396 miRNAs), PCA-based outlier removal, DESeq2 (adjusted p < 0.01, \|log2FC\| > 2), heatmaps, volcano plot, plots of the TAM 2.0 enrichment results |
| 2 | `02_machine_learning.R` | Random Forest, SVM, LASSO and a neural network trained with caret, 10-fold cross-validation with random under-sampling of tumour samples inside each training fold; ROC and precision-recall evaluation; feature importance per model; consensus ranking; bootstrap (100 iterations) and label-permutation test (20 iterations) |

## Key results

**Differential expression:** 201 differentially expressed miRNAs (DEMs) between tumour and normal tissue.

**Classification** (10-fold CV with random under-sampling, Normal as the positive class):

| Model | Accuracy | ROC-AUC | Sensitivity | Specificity | F1 | MCC |
|---|---|---|---|---|---|---|
| Random Forest | 99.78% | 100% | 100% | 99.78% | 94.12% | 0.9417 |
| SVM | 99.78% | 100% | 100% | 99.78% | 94.12% | 0.9417 |
| LASSO | 99.12% | 99.80% | 87.5% | 99.33% | 77.78% | 0.7783 |
| Neural network | 99.78% | 100% | 100% | 99.78% | 94.12% | 0.9417 |

LASSO retained 4 of the 201 DEMs: hsa-mir-19a, hsa-mir-15b, hsa-mir-379 and hsa-mir-1468.

**Top 10 candidate biomarkers** (mean of 0–100 normalised importance across the four models):

| Rank | miRNA | Combined score |
|---|---|---|
| 1 | hsa-mir-19a | 77.1 |
| 2 | hsa-mir-15b | 75.0 |
| 3 | hsa-mir-379 | 73.5 |
| 4 | hsa-mir-3614 | 60.5 |
| 5 | hsa-mir-30b | 60.3 |
| 6 | hsa-mir-592 | 59.3 |
| 7 | hsa-mir-135b | 59.3 |
| 8 | hsa-mir-3150b | 58.7 |
| 9 | hsa-mir-1468 | 58.7 |
| 10 | hsa-mir-484 | 57.7 |

Five of these (hsa-mir-19a, hsa-mir-379, hsa-mir-30b, hsa-mir-1468, hsa-mir-3150b) were not found to be previously characterised in colorectal cancer in the thesis literature review, and are candidates for experimental validation.

## How to run

1. Download the data as described in [`data/README.md`](data/README.md).
2. Install the R packages:
   ```r
   install.packages(c("data.table", "readr", "ggplot2", "reshape2", "ggfortify",
                      "circlize", "rgl", "dplyr", "scales", "caret", "randomForest", "e1071",
                      "glmnet", "nnet", "pROC", "PRROC"))
   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   BiocManager::install(c("DESeq2", "ComplexHeatmap", "EnhancedVolcano"))
   ```
3. Open R (or RStudio) with the repository root as the working directory and run the scripts in order:
   ```r
   source("scripts/01_preprocessing_differential_expression.R")
   source("scripts/02_machine_learning.R")
   ```
   Tables and saved figures are written to the working directory. Some plots (density, PCA, heatmaps, volcano, the interactive 3D PCA) open in the plot window and can be exported from there. Script 2 reads the CSV files written by script 1.

The thesis analysis was run in R 4.2.0.

**Reproducibility:** rerunning script 2 in R 4.5.1 reproduces the classification metrics exactly, and the LASSO, SVM and neural network importance scores. Random Forest importance (permutation-based) varies slightly between R and package versions, so the consensus ranking can shift a little: in the rerun, the top 3 were identical and 9 of the top 10 were the same, with hsa-mir-484 moving from 10th to 11th.

## Notes and limitations

- **Very few normal samples.** TCGA-COAD contains only 8 adjacent normal miRNA samples, against 445 tumour samples after QC. With 10 folds, each test fold holds at most one normal sample, so the near-perfect metrics should be read with caution and need confirmation in an independent cohort.
- **Feature selection before cross-validation.** The 201 DEMs were selected with DESeq2 on all samples before cross-validation, so the CV estimates are optimistic. Nesting the differential expression step inside each fold would give an unbiased estimate.
- **Fold-change direction.** As in the thesis, log2 fold changes are computed as normal relative to tumour (`contrast = c("Tissue_Type", "Normal", "Tumor")`), so a positive log2FC means higher expression in normal tissue.
- **Permutation test.** The original permutation-test code was lost; step 12 of `02_machine_learning.R` is a re-implementation, so its AUC ranges differ slightly from those in the thesis (no permuted run reaches the real AUC in either).
- **Tissue, not blood.** The candidates were found in tumour tissue. Their value as non-invasive (e.g. blood-based) markers would need to be tested in circulating miRNA data.
