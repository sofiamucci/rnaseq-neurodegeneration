library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)

# Load significant genes
res_table_sig <- readRDS("data/res_table_sig.rds")

# Strip Ensembl version suffix (same fix as 03_volcano_plot.R)
res_table_sig$gene_id_clean <- sub("\\..*", "", res_table_sig$gene_id)

# GO over-representation analysis (ORA)
go_enrichment <- enrichGO(gene = res_table_sig$gene_id_clean,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENSEMBL",
                          ont = "ALL",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.1,
                          readable = TRUE)

as.data.frame(go_enrichment)

# Dotplot of top enriched categories
dotplot(go_enrichment, showCategory = 15) + 
  ggtitle("GO enrichment: significant DE genes (Old vs AD)")
