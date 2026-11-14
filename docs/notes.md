# Project notes: regional expression of brain-disorder risk genes

## Question
Do GWAS risk genes for brain disorders show elevated expression in particular brain regions (Allen Human Brain Atlas)?

## Environment (Colab)
- Python 3.13, abagen 0.1.3, pandas 2.2.3, numpy 2.1.3
- abagen 0.1.3 needed two compatibility fixes on this stack:
  1. Patch: removed `inplace=False` from `set_axis(...)` calls in abagen/probes_.py (pandas 2 removed that argument).
  2. Shim: re-created `DataFrame.append` using `pd.concat` (pandas 2 removed it). Held in memory only, so re-run it after any restart.
- Do not install nilearn/neuromaps unpinned: it upgraded pandas to 3.x and broke abagen.

## Brain data
- Atlas: Desikan-Killiany, 83 regions (68 cortex, 15 subcortex/brainstem).
- `get_expression_data(..., lr_mirror='bidirectional')` because only 2 of 6 donors have right-hemisphere samples.
- Result: 83 regions x 15,633 genes, no empty regions. Saved as expression_dk.csv / expression_dk_named.csv.
- Primary analysis set: left hemisphere + brainstem = 42 regions (the right hemisphere is a mirrored copy, L-vs-R correlation 0.98, so it would double-count).
- Atlas cannot resolve substantia nigra, so Parkinson's cannot be tested at its classic site.

## GWAS data
- File: gwas-catalog-associations_ontology-annotated-full.zip from the GWAS Catalog FTP (releases/latest). Date downloaded: FILL IN.
- Inner file: gwas-catalog-download-associations-alt-full.tsv, 1,205,644 rows. Slim copy: gwas_slim.csv.gz.
- Gene list rules (fixed before looking at results):
  - p < 5e-8
  - rows whose MAPPED_TRAIT exactly equals the disease label (no multi-trait or polygenic-score rows)
  - MAPPED_GENE split on " - ", commas, semicolons, " x " (bare hyphens kept, e.g. HLA-DRB1)
  - keep only genes present in the brain table
  - 20-500 usable genes: more than 500 -> keep the 500 with strongest p; fewer than 20 -> excluded
  - "clean" version also drops MHC (chr6 25-34 Mb) and APOE (chr19 44.4-45.4 Mb). Coordinates assume GRCh38 (not verified).

## Statistics
- Score = mean expression of the gene set per region.
- Null = 10,000 random gene sets of the same size (seed 0), drawn excluding the real genes.
- p one-sided = (count + 1) / (10,000 + 1); BH FDR across the 42 regions.
- Matched nulls: random genes drawn from bins of equal PC1 loading (20 bins), or PC1 x PC2 (6 x 6 bins).
- A q printed as 0.000 is the permutation floor. Report z-scores.

## Pre-declared predictions (written before results; from general knowledge, NOT yet verified in the literature)
- Alzheimer's: entorhinal, hippocampus
- ALS: precentral
- Epilepsy: hippocampus, entorhinal, parahippocampal, temporal pole
- Schizophrenia: hippocampus, superior/rostral/caudal middle frontal
- Everything else: exploratory. Crohn's and coronary artery disease are negative controls.

## Results so far
Positive control (striatal markers, 9 of 10 usable): caudate z=7.34, putamen z=5.55, accumbens z=2.95 (q = 0.0021, 0.0021, 0.018). Other regions near or below zero.

Alzheimer's, original analysis (free-text trait names):
- strict 166 genes: best q 0.78 (thalamus); no-APOE 156 genes: best q 0.97
- broad list 756 usable genes (1,299 raw, above the 500 guideline): best q 0.17
- entorhinal+hippocampus: strict p=0.069, no-APOE p=0.165, broad p=0.041 (not corrected across lists; not claimed as significant)

Huntington's: GWAS file has 22 rows, 3 pass 5e-8, only FAN1 usable. Too small to test (underpowered). HTT is not enriched in striatum (caudate/putamen rank 37/38 of 42). Descriptive only.

Grid (regions with q<0.05 out of 42): standard null / PC1-matched (primary) / PC1-matched no MHC+APOE / PC1xPC2-matched
- schizophrenia 32 / 22 / 26 / 29
- bipolar 25 / 19 / 20 / 27
- major depression 27 / 14 / 20 / 15
- ADHD 10 / 0 / 0 / 1 (not robust)
- Alzheimer's, Parkinson's, ALS, MS, epilepsy, migraine, Crohn's, CAD: 0 in every version
- Autism: excluded (19 usable genes; 18 without MHC/APOE)
- Grid Alzheimer's list has 202 usable genes (mapped-trait definition) vs 166 in the original analysis (free-text definition).

Observations:
- Signal for the three robust psychiatric diseases is broad (14-29 regions), not focal. Top region changes between versions, so do not interpret it.
- PC1 explains 43.6% of regional variation (PC2 15.1%, PC3 8.7%); psychiatric standard-null profiles correlate -0.74 to -0.91 with PC1.
- Gene-list overlap between diseases is low (max Jaccard 0.12), so shared genes do not explain shared profiles.
- Clustering (primary, PC1-matched z): no clean psychiatric/neurodegenerative separation (e.g. bipolar pairs with Crohn's, schizophrenia with migraine). Null-disease profiles are mostly noise.
- Supplementary: standard-null psychiatric profiles correlate 0.79-0.92, but this is driven by the shared cortical gradient.
- 2-axis matching gave MORE significant regions than 1-axis for schizophrenia and bipolar; unexplained, possibly a narrower null. Do not read as stronger evidence.

## Analysis rule
All analyses run are reported. Primary = PC1-matched grid. Other versions are robustness checks. No further variants to be added after seeing results.

## Limitations
- 6 donors, bulk tissue (cell types averaged together)
- Probe filtering dropped real immune genes (CR1, MS4A6A, TREML2, LILRB2, RHOH...), so microglial risk is under-represented
- Nearest-gene mapping is not causal-gene assignment
- Matching on 1-2 axes does not rule out all confounds (e.g. cell-type composition)
- 11 diseases tested: disease-level multiple testing not formally corrected
- Atlas coarse; no substantia nigra
- Spatial autocorrelation nulls and a second parcellation (Schaefer) NOT done: future work

## Not yet verified (check before citing)
- Fulcher et al. 2021 paper on random-gene nulls
- Psychiatric GWAS genes enriched in neuron-expressed genes
- Predicted regions listed above
- Biological interpretation of HTT result
