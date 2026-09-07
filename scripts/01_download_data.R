# =============================================================
# 01_download_data.R
# Download and subset RNA-seq counts from GSE104704 (GEO/SRA)
# via recount3. Subset: 3 "Old" + 3 "AD" donors, pairwise 
# matched by age, PMI, and RIN (see METHODS.md for rationale).
# =============================================================

library(recount3)
library(dplyr)

# 1. Retrieve the full study (all ~31 samples)
human_projects <- available_projects()
proj_info <- subset(human_projects, project == "SRP119561")
rse <- create_rse(proj_info)

# 2. Our subset: 3 Old + 3 AD (see data/sample_metadata.csv for full details)
srr_subset <- c(
  "SRR6145425", "SRR6145429", "SRR6145431",  # Old
  "SRR6145433", "SRR6145436", "SRR6145443"   # AD
)

# recount3 stores the SRR accession in colData$external_id
rse_subset <- rse[, colData(rse)$external_id %in% srr_subset]

# 3. Sanity check: confirm we got exactly our 6 samples
colData(rse_subset)$external_id
dim(rse_subset)  # expect: genes x 6

# 4. Extract raw counts matrix and save
counts <- assay(rse_subset, "raw_counts")
dim(counts)       # sanity check: same 6 columns
saveRDS(rse_subset, "data/rse_subset.rds")
file.exists("data/rse_subset.rds")  # sanity check: saved correctly

# 5. Extract raw counts matrix and export as CSV (if prefers CSV)
counts <- assay(rse_subset, "raw_counts")
dim(counts)  # sanity check: genes x 6

write.csv(counts, "data/counts.csv", row.names = TRUE)
file.exists("data/counts.csv")  # sanity check: saved correctly
