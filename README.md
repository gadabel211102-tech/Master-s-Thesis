# GSDMB Variant Analysis Pipeline

A fully integrated bioinformatics and statistical pipeline for targeted amplicon sequencing of the **GSDMB locus (17q21)** in breast and endometrial cancer cohorts.

Developed as part of a Master's thesis (TFM), this repository implements an end-to-end workflow spanning sequencing quality control, variant annotation, haplotype analysis, clinical harmonisation, and SNP–phenotype association modelling.

---

## Overview

This pipeline performs:

- BAM integrity verification and indexing  
- Coverage quality control and technical auditing  
- Variant calling and normalisation  
- Functional annotation (Ensembl VEP + pathogenicity predictors)  
- COSMIC contextualisation  
- Variant statistics and enrichment testing  
- Haplotype phasing and association testing  
- Clinical data integration and harmonisation  
- SNP–clinical association modelling  
- Survival analysis (Kaplan–Meier + Cox proportional hazards)  
- Interactive visualisation outputs  

The workflow integrates sequencing QC, genomic annotation, and multivariable statistical modelling into a reproducible analytical framework.

---

## Repository Structure

```
.
├── 01_check_and_index.sh
├── 02_dna_qc.sh
├── 02b_more-qc.sh
├── 02c_zerocovinfo.py
├── 03_qc_visualisation.py
├── 04_technical_audit.py
├── 05_variant_calling.sh
├── 06_annotation.sh
├── 07_merge-annotations.py
├── 07b_variant-qc.py
├── 07c_cosmic_integration.py
├── 08_mapping.py
├── 09_gsdmb-only.py
├── 10_variant-stats.py
├── 11_SNPs.py
├── 12_stats-enrichment.py
├── 13_permutation_testing.py
├── 14_interactive_dashboard.py
├── 15_haplotypes.sh
├── 15_haplo.stats.r
├── 16_excel_making.py
├── 16b_excel_harmonisation.py
└── 17_snp-association.py
```

Scripts are intended to be executed sequentially (`01` → `17`).

---

# Computational Environments

The pipeline uses **three isolated environments** due to distinct toolchains required for bioinformatics processing, annotation, and statistical modelling.

---

## 1. Bioinformatics Environment (Micromamba)

Used for shell-based genomic processing:

- `01_check_and_index.sh`
- `02_dna_qc.sh`
- `02b_more-qc.sh`
- `05_variant_calling.sh`
- `15_haplotypes.sh`

### Installation

```bash
micromamba create -n bam-steps \
  samtools bcftools mosdepth \
  -c bioconda -c conda-forge

micromamba activate bam-steps
```

### Tools

- samtools — BAM integrity checking and indexing  
- bcftools — VCF normalisation, filtering, merging  
- mosdepth — Per-base and per-region coverage metrics  

---

## 2. Python Analysis Environment

Used for all Python-based processing and modelling:

- `02c_zerocovinfo.py`
- `03_qc_visualisation.py`
- `04_technical_audit.py`
- `07`–`14`
- `16_excel_making.py`
- `16b_excel_harmonisation.py`
- `17_snp-association.py`

### Installation

```bash
python3 -m venv tfm_env
source tfm_env/bin/activate
pip install -r requirements.txt
```

### Core Dependencies

- pandas  
- numpy  
- scipy  
- statsmodels  
- seaborn  
- matplotlib  
- lifelines  
- openpyxl  

This environment handles data integration, harmonisation, regression modelling, survival analysis, and visualisation.

---

## 3. VEP Annotation Environment

Used exclusively for:

- `06_annotation.sh`

### Activation

```bash
conda activate vep_env
```

### Required Resources

- Ensembl VEP (GRCh38 cache)  
- Reference FASTA (GRCh38)  
- CADD scores  
- dbNSFP  
- COSMIC annotation data  

> Note: The annotation script automatically unsets `PERL5LIB` and `PERL_LOCAL_LIB_ROOT` to avoid WSL Perl conflicts.

---

# Pipeline Stages

---

## Stage 1 — Sequencing QC

```bash
micromamba activate bam-steps
./01_check_and_index.sh
./02_dna_qc.sh
./02b_more-qc.sh

source tfm_env/bin/activate
python 02c_zerocovinfo.py
python 03_qc_visualisation.py
python 04_technical_audit.py
```

Outputs:
- Coverage metrics  
- Zero-coverage region reports  
- PCA and QC plots  
- Technical audit summary  

---

## Stage 2 — Variant Calling & Annotation

```bash
micromamba activate bam-steps
./05_variant_calling.sh

conda activate vep_env
./06_annotation.sh
```

Functional predictors applied:
- gnomAD allele frequencies  
- CADD  
- SIFT  
- PolyPhen-2  
- REVEL  
- DANN  
- COSMIC  

---

## Stage 3 — Variant Processing & Enrichment

```bash
source tfm_env/bin/activate
python 07_merge-annotations.py
python 07b_variant-qc.py
python 07c_cosmic_integration.py
python 08_mapping.py
python 09_gsdmb-only.py
python 10_variant-stats.py
python 11_SNPs.py
python 12_stats-enrichment.py
python 13_permutation_testing.py
python 14_interactive_dashboard.py
```

Produces:
- Filtered SNP lists  
- Population frequency filtering  
- Enrichment testing  
- Permutation validation  
- Interactive dashboard  

---

## Stage 4 — Haplotype Analysis

```bash
micromamba activate bam-steps
./15_haplotypes.sh

Rscript 15_haplo.stats.r
```

- BEAGLE phasing  
- Haplotype frequency estimation  
- Association modelling  

---

# Clinical Integration & Harmonisation

## 16_excel_making.py

Builds a unified SNP master dataset by:
- Merging SNP metadata  
- Joining breast and endometrial clinical data  
- Applying HER2 case-level germline logic  

Output:

```
MASTER_SNP_plus_clinical__MERGED.xlsx
```

---

## 16b_excel_harmonisation.py

Creates the harmonised analytical layer:

- Preserves all original columns (lossless merge)  
- Adds canonical `canon__*` variables  
- Standardises:
  - Yes/No encodings  
  - FIGO stage  
  - Grade  
  - ER/PR status  
  - HER2 copy number  
  - MSI status  

Output:

```
MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx
```

---

# SNP–Clinical Association Modelling

## 17_snp-association.py

Performs:

### Pre-processing
- DNA-only filtering  
- gnomAD NFE AF > 1%  
- Manifest-based sequencing confirmation  
- Technical replicate removal  

### Statistical Analyses
- Tumour vs Healthy — Fisher + Odds Ratio  
- SNP × Clinical variables:
  - Continuous → Mann–Whitney U  
  - Binary → Fisher + age-adjusted logistic regression  
  - Nominal → Chi-square  
- Genotype-dose trend testing  
- Survival analysis (Kaplan–Meier + Cox PH)  

### Multiple Testing
Benjamini–Hochberg False Discovery Rate (default: FDR < 0.10)

### Outputs
- Multi-sheet Excel results  
- Volcano plots  
- Clinical heatmaps  
- Kaplan–Meier curves  
- Sample manifest audit  

---

# Statistical Framework

- Fisher’s exact test  
- Mann–Whitney U  
- Chi-square  
- Logistic regression (age-adjusted)  
- Cox proportional hazards  
- Benjamini–Hochberg FDR correction  

---

# System Notes

- Developed on Ubuntu 24 (WSL2)  
- Python 3.11  
- Ion Torrent amplicon sequencing  
- Assumed directory structure:

```
<root>/<cohort>/<tissue>/<sample>/
```

---

# Author

Research conducted as part of a Master's Thesis, February 2026.
