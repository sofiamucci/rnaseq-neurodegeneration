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