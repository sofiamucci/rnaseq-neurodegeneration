library(dplyr)
library(httr)
library(jsonlite)

# Load previous results
res_table_sig <- readRDS("data/res_table_sig.rds")

# Load candidate genes
candidates <- res_table_sig %>%
  arrange(padj) %>%
  head(10)

candidates$label <- ifelse(is.na(candidates$gene_name), candidates$gene_id, candidates$gene_name)

# --- Query GWAS Catalog for each candidate gene ---

query_gwas_catalog <- function(gene_symbol) {
  if (is.na(gene_symbol)) return(NULL)
  
  url <- paste0(
    "https://www.ebi.ac.uk/gwas/api/search/downloads?q=ensemblMappedGenes:",
    gene_symbol,
    "&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=",
    "&genotypingfilter%5B%5D=&traitfilter%5B%5D=&dateaddedfilter=&facet=association&efo=true"
  )
  
  res <- tryCatch(GET(url), error = function(e) NULL)
  if (is.null(res) || status_code(res) != 200) return(NULL)
  
  txt <- content(res, as = "text", encoding = "UTF-8")
  if (nchar(trimws(txt)) == 0) return(NULL)
  
  tryCatch(read_tsv(txt, show_col_types = FALSE), error = function(e) NULL)
}

gwas_results <- list()
for (i in seq_len(nrow(candidates))) {
  g <- candidates$gene_name[i]
  key <- if (is.na(g)) candidates$gene_id[i] else g
  gwas_results[key] <- list(query_gwas_catalog(g))  # [key] + list(), not [[key]]
  Sys.sleep(0.5)
}

summary_table <- lapply(names(gwas_results), function(g) {
  df <- gwas_results[[g]]
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(gene = g, n_associations = 0, alzheimer_related = FALSE))
  }
  n <- nrow(df)
  alz <- any(grepl("alzheimer", df$`DISEASE/TRAIT`, ignore.case = TRUE))
  data.frame(gene = g, n_associations = n, alzheimer_related = alz)
})

summary_table <- do.call(rbind, summary_table)
summary_table

# "alzheimer_related" = exact match to Alzheimer's trait in GWAS Catalog.
# 3 genes w/o gene_name skipped (search needs symbol, not Ensembl ID).