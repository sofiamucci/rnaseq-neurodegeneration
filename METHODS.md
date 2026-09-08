# Sample Selection Rationale

To adapt this pipeline to the GSE104704 dataset (Old vs. AD groups, lateral 
temporal lobe RNA-seq), a subset of 3 biological replicates per condition 
was selected to mirror the original study design.

Samples were restricted to the "Old" (cognitively normal, aged) and "AD" 
(Alzheimer's disease) groups, excluding the "Young" group present in the 
original dataset. This choice was made to isolate disease-associated 
transcriptional changes from those attributable to normal aging, which 
would otherwise confound the comparison given the substantial age gap 
between the young and aged cohorts.

Within the Old and AD groups, donors were selected and pairwise matched on 
age at death, post-mortem interval (PMI), and RNA integrity number (RIN) 
to minimize technical and biological variability unrelated to disease 
status. Only male donors were included, as each group contained a single 
female donor; given the small sample size (n=3 per group), including a 
single unreplicated sex category would introduce unmodeled variance rather 
than allow sex to be properly controlled for as a covariate. We acknowledge 
this as a limitation, as sex-specific differences in Alzheimer's disease 
neuropathology are an active area of research.

| Pair | Old (Control) | Age | PMI (hr) | RIN | AD | Age | PMI (hr) | RIN |
|---|---|---|---|---|---|---|---|---|
| 1 | GSM2806224 | 62 | 16 | 7.9 | GSM2806232 | 64 | 16.5 | 7.8 |
| 2 | GSM2806228 | 61 | 6 | 8.0 | GSM2806242 | 64 | 5 | 7.7 |
| 3 | GSM2806230 | 72 | 17 | 7.1 | GSM2806235 | 74 | 18 | 6.3 |

## Notes on Small Sample Size

Several patterns observed in this analysis are consequences of the small 
sample size (n=3 per group) rather than pipeline errors, and are documented 
here for transparency:

- **PCA**: PC1 and PC2 together explain >50% of variance, but do not show 
  clean separation between Old and AD groups. This suggests other sources 
  of variability (individual donor differences, technical variation) may 
  dominate over disease status in this small, matched subset.

- **Differential expression outliers**: Cook's distance-based outlier 
  detection flagged a higher-than-typical percentage of genes (~8%) as 
  outliers. With n=3 per group, this method is more sensitive to individual 
  samples, which is expected behavior rather than an indication of poor 
  data quality.

- **Volcano plot gap**: after applying `lfcShrink()` (apeglm), the 
  up-regulated genes show a visible gap in log2FoldChange values (~5 to 
  ~8.5) rather than a smooth continuum. This is a known behavior of 
  shrinkage estimators under small sample sizes: genes with strong evidence 
  retain high fold-change estimates, while genes with moderate evidence are 
  shrunk more aggressively toward zero, producing a bimodal-like 
  distribution rather than a continuous one.
  
- **Effect of log2FoldChange shrinkage**: applying `lfcShrink()` (apeglm) 
  substantially reduced the number of genes passing the significance 
  threshold (padj < 0.1, |log2FC| > 2) — from 324 down / 448 up 
  (unshrunk) to 5 down / 58 up (shrunk). This confirms that most of the 
  original "significant" genes had inflated fold-change estimates driven 
  by low counts rather than genuine differential expression, underscoring 
  the importance of shrinkage correction with small sample sizes.

## Functional Interpretation

GO over-representation analysis (63 significant genes) identified enrichment 
in ion transport (SLC5A5, SLC13A4, ATP4A) and chemokine/growth factor 
signaling (CCL26, CCL19, VGF) — both pathways with documented relevance to 
neurodegeneration. Cross-referencing candidates against GWAS Catalog found 
one gene (CCL26) with prior Alzheimer's disease trait associations; the 
remaining candidates had GWAS associations with other conditions or no 
GWAS Catalog entry, which does not rule out biological relevance given the 
limited power of this small subset.

