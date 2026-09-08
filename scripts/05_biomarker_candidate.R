library(DESeq2)
library(dplyr)
library(ggplot2)
library(pROC)

# --- Load data ---
normalized_data <- readRDS("data/normalized_data.rds")
res_table_sig <- readRDS("data/res_table_sig.rds")
col_data <- readRDS("data/col_data.rds")

# --- Load data ---
normalized_data <- readRDS("data/normalized_data.rds")
res_table_sig <- readRDS("data/res_table_sig.rds")
col_data <- readRDS("data/col_data.rds")

# --- Select top candidate genes ---
# Exploratory candidate selection based on statistical significance
# (padj < 0.1, |log2FC| > 2, ranked by padj). This is NOT a validated 
# biomarker panel — see note at the end of this script.
top_n <- 10

candidates <- res_table_sig %>%
  arrange(padj) %>%
  head(top_n)

candidates[, c("gene_id", "gene_name", "log2FoldChange", "padj")]

# Create display label (gene_name if available, fallback to gene_id)
candidates$label <- ifelse(is.na(candidates$gene_name), candidates$gene_id, candidates$gene_name)
view(candidates)

# --- Prepare data for boxplots ---
# Extract normalized counts for candidate genes only
candidate_counts <- normalized_data[candidates$gene_id, ]
rownames(candidate_counts) <- candidates$label

# Reshape to long format for ggplot (one row per gene x sample)
candidate_long <- as.data.frame(candidate_counts) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(cols = -gene, names_to = "sample", values_to = "normalized_count")

# Add condition (Old/AD) from col_data
candidate_long$condition <- col_data$condition[match(candidate_long$sample, rownames(col_data))]

# --- Boxplots: expression per candidate gene, Old vs AD ---
ggplot(candidate_long, aes(x = condition, y = normalized_count, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  theme_bw() +
  theme(legend.position = "top") +
  labs(x = NULL, y = "Normalized counts")

# --- ROC/AUC per candidate gene ---
# Exploratory: with n=3/group, AUC estimates are unstable and should not
# be interpreted as validated biomarker performance.

roc_results <- lapply(unique(candidate_long$gene), function(g) {
  gene_data <- candidate_long %>% filter(gene == g)
  roc_obj <- roc(gene_data$condition, gene_data$normalized_count, quiet = TRUE)
  data.frame(gene = g, auc = as.numeric(auc(roc_obj)))
})

roc_summary <- do.call(rbind, roc_results)
roc_summary <- roc_summary[order(-roc_summary$auc), ]
roc_summary

# NOTE ON AUC = 1: these genes were selected specifically for having the 
# lowest padj in this same 6-sample dataset — meaning perfect separation 
# was guaranteed by construction, not evidence of true predictive power. 
# This is a circular evaluation (selecting on padj, then evaluating on 
# the same samples) and AUC values here should NOT be interpreted as 
# validated biomarker performance. A proper evaluation would require an 
# independent validation cohort not used in gene selection.