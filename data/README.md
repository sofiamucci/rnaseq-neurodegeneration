# Data

This directory does not contain raw sequencing data. Raw counts are 
downloaded programmatically via the `scripts/01_download_data.R` script, 
which pulls data from GEO/SRA using the `recount3` R package.

## Source
- **GEO accession**: GSE104704
- **SRA study**: SRP119561
- **Samples used**: 3 "Old" (cognitively normal) + 3 "AD" (Alzheimer's disease) 
  donors, pairwise matched by age, PMI, and RIN (see `sample_metadata.csv` 
  and [METHODS.md](../METHODS.md) for full selection rationale).

## Reproducing the download
Run `scripts/01_download_data.R`. This will fetch counts via `recount3` 
and save the subset as `data/rse_subset.rds`.