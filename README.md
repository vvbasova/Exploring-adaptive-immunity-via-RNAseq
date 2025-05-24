# Exploring tissue specific adaptive immunity features via RNA sequencing data

---

**Students:** Victoria Basova, Ryazantsev Dmitrii

**Supervisors:** Daniil Luppov, Mikhail Shugai

## Table of context

- [Project description](#project-description)
  - [Project goal and objectives](#project-goal-and-objectives)
  - [Data](#data)
  - [Results](#results)
  - [References](#references)
- [Project navigation](#project-navigation)
---

# Project description

---

This study aimed to explore immune receptor repertoires derived from bulk RNA-seq data 
across various human tissues, both healthy and cancerous. While repertoire profiling is traditionally conducted using more targeted methods such as amplicon-based sequencing or single-cell RNA-seq, our goal was to evaluate the utility of bulk RNA-seq for immunogenomic analyses using large publicly available datasets.



## Project goal and objectives

---

**Project goal**: Study T-cell and B-cell interactions in different human tissues in healthy and pathological conditions

**Objectives**: 
- Process publicly available RNA-seq data and extract TCR/BCR sequences
- Analyze antigen receptor composition based on Variable gene usage, isotype and other features
- Analyze correlation patterns between deconvolution marker gene expression and repertoire composition 
- Explore repertoire differences across healthy, abnormal and immune-privileged tissues.

## Data

---

Healthy tissue RNA-seq data were obtained via manual search from the European Nucleotide Archive (ENA) [1]. Accession numbers for all included samples are available in the `csv` tables in `MIXCR_healthy_tissues` directory.

Cancer-related immune repertoires were extracted  TCGA (The Cancer Genome Atlas) bulk RNA-seq data [2], which includes more than ten thousand tumor samples across 33 cancer types. 

To assess the antigen specificity of T-cell receptors, we utilized the VDJdb database [3], which contains curated TCR sequences with known epitope specificities. 

## Results

---


### γδ T Cells in Solid Tumors

We investigated the usage of T cell receptor (TCR) δ-chain variable genes (TRDV) to characterize the distribution of γδ T cells—an unconventional T cell subset known for their rapid immune response and ability to recognize antigens in an MHC-independent manner.

In healthy tissues, TRDV1 expression is typically dominant, while TRDV2 is enriched primarily in blood. Interestingly, in our data, the three cancer types with the highest relative TRDV2 usage were glioblastoma, leukemia, and glioma. These findings are consistent with previous studies reporting preferential infiltration of γδ T cells, especially the Vγ9Jγ2–Vδ2 subtype, into glioblastoma multiforme (GBM), potentially contributing to anti-tumor immune responses [4], [5].

![trdv_hist](./plots/TRDV_hist.png)
### Antigen Specificity Using TCR–Epitope Mapping

Using the VDJdb database of annotated TCR CDR3 sequences, we assessed antigen specificity by mapping CDR3 sequences to known human epitopes. We quantified enrichment using Fisher's exact test and visualized the results as heatmaps of log-transformed odds ratios across cancer types.

We observed significant enrichment of TCRs specific to MLANA, a melanoma-associated antigen, within melanoma samples, and to BST-2 (tetherin), an immunotherapy target in several cancer types. These results demonstrate that immune repertoires derived from bulk RNA-seq retain antigen-specific features and may hold promise for applications in cancer immunotherapy.

![fisher](./plots/fisher_plot.png)

### Gene–Gene Correlation and Survival Analysis

We further analyzed co-usage and interaction between receptor genes across different cancer types. Using correlation analysis and Kaplan–Meier survival modeling, we identified gene pairs whose expression levels were statistically correlated within specific cancers (for example, IGHV2-10 and TRBV12-3 in GBM project). In some cases, receptor usage was also linked to patient survival (for example,high TRBV4-1 abundance in GBM is associated with a significantly increased risk of death: patients with high expression show more than twice the hazard of death (HR = 2.7, p < 0.0001) compared to those with low expression.).

![kaplan](./plots/gbm_kaplan_meier.png)


![cluster](./plots/TCGA-GBM_clustermap.png)

### MAIT Cells and Their Association with Renal Cancer
Mucosal-associated invariant T (MAIT) cells are a unique T cell subset characterized by semi-invariant TCRs that recognize microbial-derived metabolites. Due to their highly conserved TCR usage, MAIT cells are readily detectable in RNA-seq–based repertoires.

Our analysis revealed a notable enrichment of MAIT cells in papillary renal cell carcinoma (KIRP) compared to other cancer types. This finding is supported by previous studies indicating elevated expression of CD161 (KLRB1)—a known MAIT cell marker and immunoregulatory molecule—in renal tumors [6]. CD161 is associated with immune infiltration, favorable immunotherapy response, and is expressed by T cells, macrophages, and tumor cells, but predominantly by MAIT cells within the T cell compartment.

Moreover, recent literature suggests a dual role for MAIT cells in cancer: while they can exhibit cytotoxic anti-tumor activity, they may also become functionally exhausted and secrete pro-tumorigenic cytokines in certain contexts [7]. The enrichment of MAIT cells in KIRP is a potentially valuable observation that warrants further investigation, particularly in the context of immunotherapy response and prognosis in kidney cancer.

![mait](./plots/mait.png)

## References

---
1. European Nucleotide Archive (ENA) [Electronic resource]. – Available at: https://www.ebi.ac.uk/ena (accessed: 24 May 2025).

2. National Cancer Institute. The Cancer Genome Atlas (TCGA) [Electronic resource]. – Available at: https://www.cancer.gov/tcga (accessed: 24 May 2025).

3. Shugay M., Bagaev D.V. et al. VDJdb: a curated database of T-cell receptor sequences with known antigen specificity // Nucleic Acids Research. – 2018. – Vol. 46, Issue D1. – P. D419–D427. – DOI: 10.1093/nar/gkx760. – Available at: https://vdjdb.cdr3.net (accessed: 24 May 2025).

4. Kang I., Kim Y., Lee H.K.
γδ T cells as a potential therapeutic agent for glioblastoma // Frontiers in Immunology. – 2023. – Vol. 14. – Article ID: 1273986. – DOI: 10.3389/fimmu.2023.1273986. – Available at: https://pubmed.ncbi.nlm.nih.gov/37928546/ (accessed: 24 May 2025).

5.  Lee M., Park C., Woo J., et al.
Preferential Infiltration of Unique Vγ9Jγ2-Vδ2 T Cells Into Glioblastoma Multiforme // Frontiers in Immunology. – 2019. – Vol. 10. – Article ID: 555. – DOI: 10.3389/fimmu.2019.00555. – Available at: https://pubmed.ncbi.nlm.nih.gov/30967876/ (accessed: 24 May 2025).

6. Li H., Zhou K., Wang K., et al.
A pan-cancer and single-cell sequencing analysis of CD161, a promising onco-immunological biomarker in tumor microenvironment and immunotherapy // Frontiers in Immunology. – 2022. – Vol. 13. – Article ID: 1040289. – DOI: 10.3389/fimmu.2022.1040289. – Available at: https://pubmed.ncbi.nlm.nih.gov/36660546/ (accessed: 24 May 2025).

7. Yigit M., Basoglu O.F., Unutmaz D.
Mucosal-associated invariant T cells in cancer: dual roles, complex interactions and therapeutic potential // Frontiers in Immunology. – 2024. – Vol. 15. – Article ID: 1369236. – DOI: 10.3389/fimmu.2024.1369236. – Available at: https://pubmed.ncbi.nlm.nih.gov/38545100/ (accessed: 24 May 2025).

   

# Project navigation

---

```angular2html
Exploring-adaptive-immunity-via-RNAseq/
│
├── MIXCR_healthy_tissues/           # Calculated CSV files with MiXCR clone data for healthy tissues
│   ├── brain_clones.csv
│   ├── colon_clones.csv
│   ├── eye_clones.csv
│   └── skin_clones.csv
│
├── notebooks/                       # Jupyter Notebooks with main analysis
│   ├── 1_airr_data_statistical_analysis.ipynb    # Statistical analysis of TCGA repertoire data
│   ├── 2_vdjdb_analysis.ipynb                    # Analysis using VDJdb epitope mapping
│   └── 3_healthy_tissues.ipynb                   # Analysis of healthy tissues repertoire extracted via MIXCR
│
├── plots/                           # Generated plots and figures
│
├── scripts/                         # Python and R scripts used in the analysis
│   ├── __init__.py                  
│   ├── correlation_functions.py     # Functions for computing and plotting correlation matrices
│   ├── regression_functions.py      # Regression functions for gene expression vs cancer types
│   ├── survival_plots.R             # R script for survival analysis visualizations
│   └── vdjdb_functions.py           # Enrichment and overlap analysis with VDJdb epitopes
│
├── R.requirements.txt               # List of R packages needed for R scripts
├── requirements.txt                 # Python dependencies
└── README.md                        # Project overview 

```





