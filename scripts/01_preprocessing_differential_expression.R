###########################################
#  - Bachelor Project                     #
#  - miRNA biomarkers in Colon Cancer     #
#  - TCGA-COAD project - miRNAseq data    #
#  - Bioinformatics analysis              #
#  - Copyright: @Malak Wael               #
###########################################

#### load libraries ########
library(readr)
library(data.table)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(ggfortify)
library(rgl)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)
library(dplyr)
library(scales)

# Change these paths to where the data is on your computer (see data/README.md)
data.path    <- "data/raw"
pheno.path   <- "data/gdc_sample_sheet.tsv"
file_pattern <- "mirbase21.mirnas.quantification.txt"


# ============================================
# LOAD miRNA COUNTS
# ============================================

generate_count_mat <- function(path, pattern, readcount_col_index, total_number_ofcolumns) {
  files <- list.files(path, pattern, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)

  # Exclude files inside /logs/ folders (GDC download cache files)
  files <- files[!grepl("/logs/", files)]

  cat("Loading", length(files), "valid files...\n")

  mat <- as.data.frame(do.call(cbind, lapply(files, function(x) fread(x, stringsAsFactors = FALSE))))
  rownames(mat) <- mat[, 1]
  mat <- as.data.frame(mat[, seq(readcount_col_index, ncol(mat), total_number_ofcolumns)])
  # Each file sits in a folder named after its GDC File ID
  colnames(mat) <- basename(dirname(files))
  return(mat)
}

exp <- generate_count_mat(data.path, file_pattern, 2, 4)


# ============================================
# LOAD SAMPLE SHEET / PHENO TABLE
# ============================================

pheno <- read_delim(pheno.path, "\t", escape_double = FALSE, trim_ws = TRUE)
table(pheno$`Tissue Type`)
names(pheno) <- sub(" ", "_", names(pheno))   # replace " " with "_" in column names

# Match file IDs to pheno$File_ID and rename columns to Sample_ID
index.files <- match(colnames(exp), pheno$File_ID)
cat("Matched:", sum(!is.na(index.files)), "out of", ncol(exp), "\n")
colnames(exp) <- pheno$Sample_ID[index.files]
head(colnames(exp))  # should show TCGA-XX-XXXX-XXX

# Group indices using pheno matched to exp
pheno_matched <- pheno[match(colnames(exp), pheno$Sample_ID), ]
normal_idx <- pheno_matched$Tissue_Type == "Normal"
tumor_idx  <- pheno_matched$Tissue_Type == "Tumor"

cat("Normal:", sum(normal_idx, na.rm = TRUE), "\n")
cat("Tumor:",  sum(tumor_idx,  na.rm = TRUE), "\n")

save(exp, pheno, file = "MODA_miRNA_Seq.RDATA")
dim(exp)
sum(is.na(exp))
table(pheno$Tissue_Type)


# ============================================
# DENSITY PLOTS - Tumor alone and Normal alone
# ============================================

normal_exp <- log2(as.matrix(exp[, which(normal_idx)]) + 1)
tumor_exp  <- log2(as.matrix(exp[, which(tumor_idx)]) + 1)

normal_melted <- melt(normal_exp)
colnames(normal_melted) <- c("miRNA", "Sample", "Expression")
normal_melted$Group <- "Normal"

tumor_melted <- melt(tumor_exp)
colnames(tumor_melted) <- c("miRNA", "Sample", "Expression")
tumor_melted$Group <- "Tumor"

combined_data <- rbind(normal_melted, tumor_melted)

# Plot 1: Combined density plot (Normal vs Tumor)
ggplot(combined_data, aes(x = Expression, color = Group, fill = Group)) +
  geom_density(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Normal" = "steelblue", "Tumor" = "tomato")) +
  scale_fill_manual(values = c("Normal" = "steelblue", "Tumor" = "tomato")) +
  labs(title = "Density Plot - Normal vs Tumor Tissue",
       x = "log2(counts + 1)",
       y = "Density") +
  theme_bw() +
  theme(legend.position = "top")

# Plot 2: Normal tissue density plot (alone)
ggplot(normal_melted, aes(x = Expression)) +
  geom_density(color = "steelblue", fill = "steelblue", alpha = 0.5, size = 1) +
  labs(title = paste0("Density Plot - Nearby Normal Tissue (", ncol(normal_exp), " samples)"),
       x = "log2(counts + 1)",
       y = "Density") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Plot 3: Tumor tissue density plot (alone)
ggplot(tumor_melted, aes(x = Expression)) +
  geom_density(color = "tomato", fill = "tomato", alpha = 0.5, size = 1) +
  labs(title = paste0("Density Plot - Tumor Tissue (", ncol(tumor_exp), " samples)"),
       x = "log2(counts + 1)",
       y = "Density") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


# ============================================
# BOXPLOTS - 5 separate tumor boxplots + 1 normal boxplot
# Each saves as its own image file
# ============================================

tumor_samples  <- colnames(exp)[which(tumor_idx)]
normal_samples <- colnames(exp)[which(normal_idx)]

n_tumor  <- length(tumor_samples)
n_normal <- length(normal_samples)

# Divide tumor samples into 5 equal groups
samples_per_group <- floor(n_tumor / 5)
tumor_groups <- list()

for (i in 1:5) {
  start_idx <- (i - 1) * samples_per_group + 1
  end_idx   <- i * samples_per_group
  tumor_groups[[i]] <- tumor_samples[start_idx:end_idx]
}

for (i in 1:5) {
  tumor_group_exp <- log2(as.matrix(exp[, tumor_groups[[i]]]) + 1)

  plot_data <- melt(tumor_group_exp)
  colnames(plot_data) <- c("miRNA", "Sample", "Expression")

  p <- ggplot(plot_data, aes(x = Sample, y = Expression)) +
    geom_boxplot(fill = "tomato", alpha = 0.7, outlier.shape = NA) +
    labs(title = paste0("Tumor Group ", i, " (n=", ncol(tumor_group_exp), " samples)"),
         x = "Samples",
         y = "log2(counts + 1)") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.title = element_text(size = 12),
          panel.grid.minor = element_blank())

  ggsave(filename = paste0("Tumor_Group_", i, "_boxplot.png"),
         plot = p, width = 8, height = 6, dpi = 300)
  cat("Saved: Tumor_Group_", i, "_boxplot.png\n", sep = "")
}

normal_plot_data <- melt(normal_exp)
colnames(normal_plot_data) <- c("miRNA", "Sample", "Expression")

p_normal <- ggplot(normal_plot_data, aes(x = Sample, y = Expression)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.shape = NA) +
  labs(title = paste0("Normal Tissue (n=", ncol(normal_exp), " samples)"),
       x = "Samples",
       y = "log2(counts + 1)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.ticks.x = element_line(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank())

ggsave(filename = "Normal_Group_boxplot.png",
       plot = p_normal, width = 8, height = 6, dpi = 300)
cat("Saved: Normal_Group_boxplot.png\n")


# ============================================
# FILTERING
# ============================================

# Remove lowly expressed miRNAs: total read count across all samples <= 551
exp <- exp[which(rowSums(exp) > 551), ]
dim(exp)

# Remove the bottom 25% least variable miRNAs
variances    <- apply(exp, 1, var)
threshold    <- quantile(variances, 0.25)
exp_filtered <- exp[variances > threshold, ]

cat("miRNAs before:", nrow(exp), "\n")
cat("miRNAs after:", nrow(exp_filtered), "\n")


# ============================================
# PCA + OUTLIER REMOVAL
# ============================================

# Step 1 — Run PCA
exp_pca <- prcomp(t(exp_filtered), scale = TRUE)
summary(exp_pca)

# Step 2 — 2D plot
autoplot(exp_pca,
         data   = pheno_matched,
         colour = 'Tissue_Type',
         frame  = FALSE)

# Step 3 — 3D plot
color <- rep("tomato", nrow(pheno_matched))
color[which(pheno_matched$Tissue_Type == "Normal")] <- "steelblue"

plot3d(exp_pca$x[, 1:3],
       col = color, pch = 20, size = 12,
       axes = FALSE, xlab = "", ylab = "", zlab = "")
axes3d(edges = "bbox", box = TRUE, col = "gray50")
title3d(xlab = "PC1", ylab = "PC2", zlab = "PC3", line = 3)
grid3d(c("x", "y", "z"), col = "lightgray")
legend3d("topright",
         legend = c("Tumor", "Normal"),
         col = c("tomato", "steelblue"),
         pch = 20, cex = 0.6, inset = 0.05)

# Step 4 — Remove PC2 outliers (more than 3 SD from the mean)
outliers      <- which(abs(scale(exp_pca$x[, 2])) > 3)
outlier_names <- rownames(exp_pca$x)[outliers]
cat("Removing:", outlier_names, "\n")

exp_clean   <- exp_filtered[, !colnames(exp_filtered) %in% outlier_names]
pheno_clean <- pheno_matched[!pheno_matched$Sample_ID %in% outlier_names, ]
dim(exp_clean)

# Step 5 — Redo PCA after outlier removal
exp_pca2 <- prcomp(t(exp_clean), scale = TRUE)

autoplot(exp_pca2,
         data   = pheno_clean,
         colour = 'Tissue_Type',
         frame  = FALSE) +
  scale_color_manual(values = c("Normal" = "steelblue", "Tumor" = "tomato"))

color <- rep("tomato", nrow(pheno_clean))
color[which(pheno_clean$Tissue_Type == "Normal")] <- "steelblue"

plot3d(exp_pca2$x[, 1:3],
       col = color, pch = 20, size = 12,
       axes = FALSE, xlab = "", ylab = "", zlab = "")
axes3d(edges = "bbox", box = TRUE, col = "gray50")
title3d(xlab = "PC1", ylab = "PC2", zlab = "PC3", line = 3)
grid3d(c("x", "y", "z"), col = "lightgray")
legend3d("topright",
         legend = c("Tumor", "Normal"),
         col = c("tomato", "steelblue"),
         pch = 20, cex = 0.6, inset = 0.05)


# ============================================
# DIFFERENTIAL EXPRESSION ANALYSIS (DESeq2)
# ============================================

# Step 1 — Check conditions
table(pheno_clean$Tissue_Type)

# Step 2 — Prepare integer matrix
exp_deseq       <- round(as.matrix(exp_clean))
mode(exp_deseq) <- "integer"

# Step 3 — Match pheno to exp columns
pheno_deseq <- pheno_clean[match(colnames(exp_deseq), pheno_clean$Sample_ID), ]
pheno_deseq$Tissue_Type <- factor(pheno_deseq$Tissue_Type, levels = c("Normal", "Tumor"))

# Remove samples without matched metadata
keep <- !is.na(pheno_deseq$Sample_ID)
exp_deseq   <- exp_deseq[, keep]
pheno_deseq <- pheno_deseq[keep, ]
cat("Samples after removal:", ncol(exp_deseq), "\n")

write.csv(exp_clean,   file = "exp_clean_matrix.csv", row.names = TRUE)
write.csv(pheno_clean, file = "pheno_clean_table.csv", row.names = FALSE)

# Verify they match — MUST be TRUE before continuing
cat("Samples match:", all(colnames(exp_deseq) == pheno_deseq$Sample_ID), "\n")
cat("Normal:", sum(pheno_deseq$Tissue_Type == "Normal"), "\n")
cat("Tumor:",  sum(pheno_deseq$Tissue_Type == "Tumor"),  "\n")

# Step 4 — Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = exp_deseq,
                              colData   = pheno_deseq,
                              design    = ~ Tissue_Type)

# Step 5 — Run DESeq2
dds.run <- DESeq(dds)

# Step 6 — Get results
# log2FoldChange is Normal relative to Tumor (cond1 / cond2), as in the thesis.
cond1 <- "Normal"
cond2 <- "Tumor"
res <- results(dds.run, contrast = c("Tissue_Type", cond1, cond2))

# Step 7 — Remove NAs
res <- res[complete.cases(res), ]
summary(res)

res.df <- as.data.frame(res)

# Apply constraints: padj < 0.01 and |log2FC| > 2
res.dems <- res.df[res.df$padj < 0.01 & abs(res.df$log2FoldChange) > 2, ]
res.dems <- res.dems[order(res.dems$padj), ]
dems <- rownames(res.dems)

cat("Total DEmiRNAs:", length(dems), "\n")
cat("Upregulated (log2FC > 0):", sum(res.dems$log2FoldChange > 0), "\n")
cat("Downregulated (log2FC < 0):", sum(res.dems$log2FoldChange < 0), "\n")

write.table(dems, file = "dems.txt", quote = F, col.names = F, row.names = F)

# Normalized expression (log2(normalized counts + 1))
ntd      <- normTransform(dds.run)
exp.norm <- assay(ntd)
dems.exp <- exp.norm[dems, ]
dim(dems.exp)

# Top 50 dysregulated miRNAs
res.df_sorted <- res.df[order(res.df$padj, decreasing = FALSE), ]
dems.top50    <- rownames(res.df_sorted)[1:50]
res.top50_df  <- res.df[dems.top50, ]

upregulated   <- res.top50_df[res.top50_df$log2FoldChange > 0, ]
downregulated <- res.top50_df[res.top50_df$log2FoldChange < 0, ]
upregulated   <- upregulated[order(upregulated$log2FoldChange, decreasing = TRUE), ]
downregulated <- downregulated[order(downregulated$log2FoldChange, decreasing = FALSE), ]
res_top50_ordered <- rbind(upregulated, downregulated)

write.csv(res.df,            file = "all_miRNAs_results.csv")
write.csv(res.dems,          file = "significant_DEMs.csv")
write.table(dems,            file = "DEMs_names_list.txt", quote = F, col.names = F, row.names = F)
write.csv(dems.exp,          file = "DEMs_normalized_expression.csv")
write.csv(res_top50_ordered, file = "top50_miRNAs_ordered.csv")
write.table(dems.top50,      file = "top50_names.txt", quote = F, col.names = F, row.names = F)


# ============================================
# HEATMAPS
# ============================================

col_annotation <- HeatmapAnnotation(
  Group = pheno_deseq$Tissue_Type,
  col = list(Group = c("Normal" = "steelblue", "Tumor" = "tomato")),
  annotation_name_side = "right"
)
col_fun <- colorRamp2(c(-2, 0, 2), c("steelblue", "white", "tomato"))

# Heatmap 1: Top 50 DEMs
dems.top50.exp    <- exp.norm[dems.top50, ]
dems.top50.scaled <- t(scale(t(dems.top50.exp)))

Heatmap(dems.top50.scaled,
        name = "Z-score",
        col = col_fun,
        top_annotation = col_annotation,
        show_column_names = FALSE,
        show_row_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 5),
        column_title = "Top 50 Differentially Expressed miRNAs",
        heatmap_legend_param = list(title = "Z-score"))

# Heatmap 2: All DEMs
all_dems.scaled <- t(scale(t(dems.exp)))

Heatmap(all_dems.scaled,
        name = "Z-score",
        col = col_fun,
        top_annotation = col_annotation,
        show_column_names = FALSE,
        show_row_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = paste0("All DEMs (n=", nrow(dems.exp), ")"),
        heatmap_legend_param = list(title = "Z-score"))


# ============================================
# PCA - DEMs only (2D and 3D)
# ============================================

pca_dems <- prcomp(t(dems.exp), scale = TRUE)

autoplot(pca_dems,
         data   = pheno_deseq,
         colour = "Tissue_Type",
         frame  = FALSE) +
  scale_color_manual(values = c("Normal" = "steelblue", "Tumor" = "tomato")) +
  ggtitle("PCA — Differentially Expressed miRNAs (DEMs)")

colors_dems <- ifelse(pheno_deseq$Tissue_Type == "Normal", "steelblue", "tomato")

plot3d(pca_dems$x[, 1:3],
       col = colors_dems, pch = 19, size = 8,
       axes = FALSE, xlab = "", ylab = "", zlab = "")
axes3d(edges = "bbox", box = TRUE, col = "gray50")
title3d(xlab = "PC1", ylab = "PC2", zlab = "PC3", line = 3)
grid3d(c("x", "y", "z"), col = "lightgray")
legend3d("topright", legend = c("Tumor", "Normal"),
         col = c("tomato", "steelblue"), pch = 19, cex = 0.8)


# ============================================
# VOLCANO PLOT (EnhancedVolcano)
# ============================================

# Colours based on the DEM criteria (padj < 0.01 & |log2FC| > 2)
keyvals <- ifelse(
  res.df$padj < 0.01 & res.df$log2FoldChange > 2, 'red',          # upregulated
  ifelse(
    res.df$padj < 0.01 & res.df$log2FoldChange < -2, 'green',     # downregulated
    'grey'                                                         # not significant
  )
)
names(keyvals)[keyvals == 'red']   <- 'Upregulated'
names(keyvals)[keyvals == 'green'] <- 'Downregulated'
names(keyvals)[keyvals == 'grey']  <- 'Not significant'

selectLab <- rownames(res.df)[res.df$padj < 0.01 & abs(res.df$log2FoldChange) > 2]

EnhancedVolcano(
  res.df,
  lab = rownames(res.df),
  selectLab = selectLab,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Normal vs Tumor',
  pCutoff = 0.01,
  FCcutoff = 3,
  colCustom = keyvals,
  pointSize = 3.0,
  labSize = 5.0,
  colAlpha = 1,
  legendPosition = 'right',
  subtitle = NULL,
  caption = paste0("Total = ", nrow(res.df), " miRNAs"),
  legendLabSize = 12,
  legendIconSize = 4.0
)

all_mirnas <- rownames(res.df)
write.table(all_mirnas, file = "mirnas.txt", quote = F, col.names = F, row.names = F)


# ============================================
# REGULATION TABLES (Up/Down)
# ============================================

res.dems$regulation <- ifelse(res.dems$log2FoldChange > 0, "Upregulated", "Downregulated")
head(res.dems[, c("log2FoldChange", "padj", "regulation")])

res.top50_df$regulation <- ifelse(res.top50_df$log2FoldChange > 0, "Upregulated", "Downregulated")
print(res.top50_df[, c("log2FoldChange", "padj", "regulation")])

write.csv(res.dems,     file = "significant_DEMs_up_down.csv")
write.csv(res.df,       file = "res_all.csv")
write.csv(res.top50_df, file = "top50_miRNAs_up_down.csv")

top50_clean <- data.frame(
  miRNA = rownames(res.top50_df),
  log2FoldChange = round(res.top50_df$log2FoldChange, 2),
  padj = format(res.top50_df$padj, scientific = TRUE, digits = 2),
  Regulation = res.top50_df$regulation
)
print(top50_clean)

cat("Significant DEMs:", nrow(res.dems), "\n")
cat("  Upregulated:", sum(res.dems$regulation == "Upregulated"), "\n")
cat("  Downregulated:", sum(res.dems$regulation == "Downregulated"), "\n\n")
cat("Top 50 miRNAs:", nrow(res.top50_df), "\n")
cat("  Upregulated:", sum(res.top50_df$regulation == "Upregulated"), "\n")
cat("  Downregulated:", sum(res.top50_df$regulation == "Downregulated"), "\n")


# ============================================
# TAM 2.0 ENRICHMENT VISUALIZATION
# ============================================
# Enrichment was run on the DEMs with the TAM 2.0 web tool
# (http://www.lirmed.com/tam2/). Values below are copied from its output.

# FUNCTIONS DATA — Top 10 colon-relevant (FDR < 0.05)
functions_data <- data.frame(
  Term  = c(
    "Regulation of Stem Cell",
    "Cell Cycle",
    "Immune Response",
    "Epithelial-to-Mesenchymal Transition",
    "Inflammation",
    "Angiogenesis",
    "Cell Death",
    "Apoptosis",
    "Cell Differentiation",
    "Cell Proliferation"
  ),
  Count = c(36, 35, 44, 35, 55, 31, 40, 53, 27, 36),
  Fold  = c(1.43, 1.24, 1.38, 1.39, 1.38, 1.48, 1.38, 1.47, 1.42, 1.39),
  FDR   = c(2.10e-3, 2.90e-2, 2.13e-3, 4.43e-3, 5.81e-4, 2.03e-3, 2.96e-3, 7.33e-5, 7.68e-3, 4.11e-3)
)
functions_data$logFDR <- -log10(functions_data$FDR)
functions_data$Label  <- paste0(functions_data$Term, " (", functions_data$Count, ")")

# DISEASES DATA — Top 10 (FDR < 0.05)
diseases_data <- data.frame(
  Term  = c(
    "Carcinoma, Colon",
    "Carcinoma, Breast",
    "Carcinoma, Hepatocellular",
    "Carcinoma, Gastric",
    "Glioma",
    "Carcinoma, Lung",
    "Osteosarcoma",
    "Colon Neoplasms",
    "Adenocarcinoma, Colon",
    "Squamous Cell Carcinoma, Esophageal"
  ),
  Count = c(115, 125, 131, 105, 91, 95, 80, 42, 11, 63),
  Fold  = c(1.243, 1.288, 1.294, 1.242, 1.366, 1.468, 1.383, 1.532, 1.521, 1.414),
  FDR   = c(3.99e-5, 3.63e-7, 5.83e-8, 1.53e-4, 2.90e-6, 5.22e-9, 1.02e-5, 1.13e-4, 2.13e-2, 5.85e-5)
)
diseases_data$logFDR <- -log10(diseases_data$FDR)
diseases_data$Label  <- paste0(diseases_data$Term, " (", diseases_data$Count, ")")

# PLOT 1: Bar Plot — Diseases
p1 <- ggplot(diseases_data,
             aes(x = reorder(Label, Fold), y = Fold, fill = logFDR)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_gradient(low = "steelblue", high = "tomato", name = "-log10(FDR)") +
  labs(
    title    = "Top Enriched Diseases",
    subtitle = "TAM 2.0 Enrichment Analysis - Colon Cancer DEMs",
    x        = NULL,
    y        = "Fold Enrichment"
  ) +
  theme_bw() +
  theme(
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle      = element_text(hjust = 0.5),
    axis.text.y        = element_text(size = 10),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank()
  )

ggsave("TAM2_bar_diseases.png", p1, width = 10, height = 7, dpi = 300)
cat("Saved: TAM2_bar_diseases.png\n")

# PLOT 2: Bubble Plot — Functions
p2 <- ggplot(functions_data,
             aes(x = Fold, y = reorder(Label, Fold), size = Count, color = logFDR)) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(range = c(4, 14), name = "miRNA Count") +
  scale_color_gradient(low = "steelblue", high = "tomato", name = "-log10(FDR)") +
  labs(
    title    = "Bubble Plot - Enriched Functions",
    subtitle = "Size = Number of miRNAs, Color = Significance",
    x        = "Fold Enrichment",
    y        = NULL
  ) +
  theme_bw() +
  theme(
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle      = element_text(hjust = 0.5),
    axis.text.y        = element_text(size = 10),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey90")
  )

ggsave("TAM2_bubble_functions.png", p2, width = 10, height = 7, dpi = 300)
cat("Saved: TAM2_bubble_functions.png\n")
