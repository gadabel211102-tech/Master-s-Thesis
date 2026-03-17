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
- Variant statistics and enrichment testing  
- Haplotype phasing and association testing  
- Clinical data integration and harmonisation  
- SNP–clinical association modelling  
- Haplotype–clinical association modelling  
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
├── 17_snp-association.py
└── 18_haplotype-association.py
```

Scripts are intended to be executed sequentially (`01` → `18`).

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
- `18_haplotype-association.py`

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
./01_check_and_index.sh -i "$HOME/tfm/endometrium/normal/dna" -t 8
./02_dna_qc.sh \
  -i "/home/gadeaalonsoj/tfm/breast/tumour/dna" \
  -o "/home/gadeaalonsoj/tfm/breast/tumour/dna_qc" \
  -r "/home/gadeaalonsoj/ref_alt/hg38_alt.fa" \
  -b "/home/gadeaalonsoj/tfm/dna_bed/IAD255368_167_Submitted.bed" \
  -t 8
./02b_more-qc.sh \
  -q /home/gadeaalonsoj/tfm/breast/tumour/dna_qc \
  -m /home/gadeaalonsoj/tfm/manifests/breast-tumour-pass_manifest.txt \
  -z /home/gadeaalonsoj/tfm/breast/tumour/zero_cov \

source tfm_env/bin/activate
python 02c_zerocovinfo.py --zero /home/gadeaalonsoj/tfm/endometrium/tumour/zero_coverage/zero_cov_tumour.tsv --qc /home/gadeaalonsoj/tfm/endometrium/tumour/dna_qc/qc_summary.tsv
python 03_qc_visualisation.py
python 04_technical_audit.py both
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
bash 05_variant_calling.sh \
  -m ~/tfm/manifests/pass_breast_tumour.txt \
  -r ~/ref_alt/hg38_canonical.fa \
  -b ~/tfm/dna_bed/IAD255368_167_Submitted.bed \
  -o ~/tfm/breast/tumour/dna_calls \
  -t 8 \

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

Performs SNP-level association analysis using individual variant genotypes as the
exposure variable. Mirrors the statistical framework of script 18 but operates at
single-variant resolution rather than haplotype level.

### Pre-processing
- DNA-only filtering  
- gnomAD NFE AF > 1%  
- Manifest-based sequencing confirmation  
- Technical replicate removal  

### Statistical Analyses
- Tumour vs Healthy — Fisher's exact test + Odds Ratio  
- SNP × Clinical variables:
  - Continuous → Mann–Whitney U  
  - Binary → Fisher's exact + age- and BMI-adjusted logistic regression  
  - Nominal → Chi-square  
- Genotype-dose trend testing (heterozygous / homozygous contrasts)  
- Survival analysis (Kaplan–Meier + age-adjusted Cox PH) for OS and PFS  
- Cancer risk — case-control logistic regression (tumour vs healthy)  

### Multiple Testing
Benjamini–Hochberg False Discovery Rate (default: FDR < 0.10)

### Outputs
- `17_SNP_Clinical_Association_Results.xlsx` — multi-sheet results workbook  
- Volcano plots (raw p and log₂ OR)  
- Clinical association heatmaps (raw p and FDR)  
- Forest plots (binary outcomes, cancer risk)  
- Kaplan–Meier curves (PDF)  
- Summary panel figure  
- Sample manifest audit  

### Run

```bash
python 17_snp-association.py \
  --gsdmb   /path/to/GSDMB_Annotated_Report_Fixed.xlsx \
  --master  /path/to/MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx \
  --out_dir /path/to/output/
```

---

# Haplotype–Clinical Association Modelling

## 18_haplotype-association.py

Mirrors script 17 but uses **BEAGLE-phased haplotype carrier status** as the
exposure variable instead of individual SNP genotypes. Haplotype strings are
constructed from the established SNP positions (gnomAD NFE AF > 1%) identified
in script 11; each sample carries two haplotype strings (one per chromosome)
and is classified as a carrier of a given haplotype if it appears on at least
one chromosome (dosage ≥ 1). Only haplotypes with a global frequency ≥ 2% are
tested.

### Data Sources
- `phased_genotypes.tsv` — BEAGLE output from `15_haplotypes.sh`  
- `MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx` — harmonised clinical master  
- `GSDMB_Annotated_Report_Fixed.xlsx` — source of established SNP positions  

### Statistical Analyses
- Tumour vs Healthy — Fisher's exact test per haplotype  
- Haplotype × Clinical variables (per cohort):
  - Continuous → Mann–Whitney U  
  - Binary → Fisher's exact + age- and BMI-adjusted logistic regression  
  - Nominal → Chi-square  
- Haplotype-dose trend testing (dosage 0 / 1 / 2) for nominally significant hits  
- Survival analysis (Kaplan–Meier + age-adjusted Cox PH) for OS and PFS  
- Cancer risk — case-control logistic regression (tumour vs healthy)  

### Multiple Testing
Benjamini–Hochberg False Discovery Rate (default: FDR < 0.10)

### Outputs
- `18_Haplo_Clinical_Association_Results.xlsx` — multi-sheet results workbook:
  - `haplotype_frequencies`
  - `tumour_vs_healthy`
  - `breast_clinical_assoc`
  - `endo_clinical_assoc`
  - `haplotype_dose`
  - `survival_cox`
  - `cancer_risk`
  - `summary_significant`
  - `sample_manifest`
- `18_Haplo_Volcano_raw_p.png` — tumour vs healthy volcano plot  
- `18_Haplo_Heatmap_Breast.png` — dual-panel heatmap (−log₁₀p + effect size), breast cohort  
- `18_Haplo_Heatmap_Endometrial.png` — dual-panel heatmap, endometrial cohort (raw p)  
- `18_Haplo_Heatmap_Endometrial_FDR.png` — dual-panel heatmap, endometrial cohort (FDR)  
- `18_Haplo_Forest_Breast.png` — forest plot, breast binary outcomes  
- `18_Haplo_Forest_Endometrial.png` — forest plot, endometrial binary outcomes  
- `18_Haplo_Forest_CancerRisk.png` — cancer risk forest plot with 95% CI whiskers  
- `18_KM_Curves_All_Cohorts.pdf` — Kaplan–Meier curves for all haplotypes and endpoints  

### Run

```bash
python 18_haplotype-association.py \
  --phased   /path/to/phased_genotypes.tsv \
  --master   /path/to/MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx \
  --annot    /path/to/GSDMB_Annotated_Report_Fixed.xlsx \
  --out_dir  /path/to/output/ \
  --min_freq 0.02
```

---

# Statistical Framework

- Fisher's exact test  
- Mann–Whitney U  
- Chi-square  
- Logistic regression (age-adjusted and age+BMI-adjusted)  
- Cox proportional hazards (age-adjusted)  
- Kaplan–Meier survival estimation  
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

Research conducted as part of a Master's Thesis
