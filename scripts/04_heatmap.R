library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(dplyr)

# --- Load significant genes and normalized counts ---
normalized_data <- readRDS("data/normalized_data.rds")
res_table_sig <- readRDS("data/res_table_sig.rds")

# --- Subset normalized counts to significant genes only ---
data_hm <- normalized_data[rownames(normalized_data) %in% res_table_sig$gene_id, ]

# --- Replace gene IDs with gene symbols for readability ---
gene_name_map <- setNames(res_table_sig$gene_name, res_table_sig$gene_id)
rownames(data_hm) <- gene_name_map[rownames(data_hm)]

data_hm <- as.matrix(data_hm)

# --- Heatmap ---
heat_color <- brewer.pal(7, "RdYlGn")

pheatmap::pheatmap(data_hm,
                   color = heat_color,
                   breaks = seq(-2, 2, length.out = 7),
                   border_color = "grey25",
                   cellwidth = 15,
                   cluster_rows = TRUE,
                   show_rownames = TRUE,
                   cutree_cols = 2,
                   scale = "row")
