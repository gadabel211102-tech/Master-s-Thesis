# GSDMB Variant Analysis Pipeline

A bioinformatics pipeline for targeted amplicon sequencing analysis of the **GSDMB locus** on chromosome 17q21, applied to cancer cohorts (breast and endometrial). The pipeline covers everything from raw BAM quality control through to haplotype association testing, with interactive visualisation outputs.

---

## Overview

This pipeline was developed as part of a master's thesis (*TFM*) and processes Ion Torrent amplicon sequencing data across multiple tumour cohorts. It performs:

- BAM file integrity checking and indexing
- Per-sample and per-amplicon quality control (coverage, uniformity, GC bias)
- Variant calling, normalisation, and filtering
- Functional annotation via Ensembl VEP
- COSMIC Cancer Gene Census contextualisation of detected germline variants
- Downstream variant statistics, enrichment, and permutation testing
- Clonality and tumour heterogeneity analysis
- Haplotype phasing and association statistics
- Interactive HTML dashboard for result exploration

---

## Repository Structure

```
.
├── 01_check_and_index.sh        # BAM integrity check and indexing
├── 02_dna_qc.sh                 # Primary DNA quality control (mosdepth)
├── 02b_more-qc.sh               # Extended QC metrics
├── 02c_zerocovinfo.py           # Zero-coverage region reporting
├── 03_qc_visualisation.py       # QC visualisation (violin, PCA, heatmaps, landscapes)
├── 04_technical_audit.py        # Technical audit report generation
├── 05_variant_calling.sh        # Variant calling, normalisation, and filtering
├── 06_annotation.sh             # VEP functional annotation
├── 07_merge-annotations.py      # Merge VEP annotations into master table
├── 07b_variant-qc.py            # Post-annotation variant QC
├── 07c_cosmic_integration.py     # COSMIC Cancer Gene Census contextualisation (germline)
├── 08_mapping.py                # Variant-to-gene mapping
├── 09_gsdmb-only.py             # GSDMB locus-specific extraction
├── 10_variant-stats.py          # Descriptive variant statistics
├── 11_SNPs.py                   # SNP filtering and established variant identification
├── 12_stats-enrichment.py       # Statistical enrichment analysis
├── 13_permutation_testing.py    # Permutation-based significance testing
├── 15_clonality_analysis.py     # Tumour heterogeneity and clonality assessment
├── 16_interactive_dashboard.py  # Interactive Plotly HTML dashboard
├── 17_haplotypes.sh             # BEAGLE haplotype phasing
└── 17_haplo_stats.r             # Haplotype-phenotype association (haplo.stats)
```

---

## Environments

This pipeline uses **three separate environments** due to the different tooling requirements of each stage.

### 1. Bioinformatics Environment (Micromamba)

Used for all shell-based steps (scripts `01`–`06`, `17_haplotypes.sh`).

```bash
micromamba create -n bam-steps samtools bcftools mosdepth -c bioconda -c conda-forge
micromamba activate bam-steps
```

**Tools included:**

| Tool | Purpose |
|---|---|
| `samtools` | BAM integrity checking, indexing, pileup |
| `bcftools` | VCF normalisation, filtering, merging |
| `mosdepth` | Fast per-base and per-region coverage calculation |

---

### 2. Visualisation Environment (Python venv)

Used for all Python analysis and visualisation scripts (`02c`, `03`–`16`).

```bash
python3 -m venv tfm_env
source tfm_env/bin/activate
pip install -r requirements.txt
```

See [`requirements.txt`](requirements.txt) for the full package list.

---

### 3. VEP Annotation Environment (Conda)

Used exclusively for script `06_annotation.sh`. Must be activated before running annotation.

```bash
conda activate vep_env
```

**Requires a working VEP installation** with the following plugins and cache files:

| Resource | Description |
|---|---|
| VEP cache | `homo_sapiens`, GRCh38, version 112 |
| Reference FASTA | `Homo_sapiens.GRCh38.dna.toplevel.fa.gz` |
| CADD | `CADD_GRCh38_whole_genome_SNVs.tsv.gz` |
| dbNSFP | `dbNSFP5.3.1a_grch38.gz` |

> ⚠️ **WSL users:** The annotation script unsets `PERL5LIB` and `PERL_LOCAL_LIB_ROOT` at runtime to avoid Perl environment conflicts common in WSL. This is handled automatically within `06_annotation.sh`.

---

## Pipeline Walkthrough

### Stage 1: QC

```bash
micromamba activate bam-steps

# Check BAM integrity and index
./01_check_and_index.sh -i /path/to/bams -t 4

# Run mosdepth coverage QC
./02_dna_qc.sh

# Extended QC metrics
./02b_more-qc.sh

# Zero-coverage region reporting
source tfm_env/bin/activate
python 02c_zerocovinfo.py
```

### Stage 2: Visualisation and Audit

```bash
source tfm_env/bin/activate
python 03_qc_visualisation.py   # Violin plots, PCA, coverage heatmaps
python 04_technical_audit.py    # Technical audit report
```

### Stage 3: Variant Calling

```bash
micromamba activate bam-steps
./05_variant_calling.sh \
  -m manifest.txt \
  -r ref/GRCh38.fa \
  -b targets.bed \
  -o output/
```

### Stage 4: Annotation

```bash
conda activate vep_env
./06_annotation.sh
```

### Stage 5: Downstream Analysis

```bash
source tfm_env/bin/activate
python 07_merge-annotations.py
python 07b_variant-qc.py
python 07c_cosmic_integration.py   # Biological context before SNP analysis
python 08_mapping.py
python 09_gsdmb-only.py
python 10_variant-stats.py
python 11_SNPs.py
python 12_stats-enrichment.py
python 13_permutation_testing.py
python 15_clonality_analysis.py
python 16_interactive_dashboard.py
```

### Stage 6: Haplotype Analysis

```bash
micromamba activate bam-steps
./17_haplotypes.sh   # BEAGLE phasing

Rscript 17_haplo_stats.r  # Association testing
```

---

## Annotations Applied

The VEP annotation step (`06_annotation.sh`) applies the following functional predictors:

| Annotation | Description |
|---|---|
| **CADD** | Combined Annotation Dependent Depletion – deleteriousness score |
| **SIFT** | Sorting Intolerant From Tolerant – functional impact prediction |
| **PolyPhen-2** | Polymorphism Phenotyping v2 – pathogenicity prediction |
| **REVEL** | Rare Exome Variant Ensemble Learner |
| **DANN** | Deep Neural Network pathogenicity score |
| **gnomAD** | Population allele frequencies |
| **COSMIC** | Catalogue of Somatic Mutations in Cancer |

---

## Configuration

Most scripts define configuration variables at the top of the file (typically under a `# --- CONFIGURATION ---` or `PATHS` block). Before running, review and update:

- `ROOT_DIR` / `BASE_PATH` – root directory for cohort data
- `REF_FA` – path to the reference FASTA
- `CACHE_DIR` – VEP cache directory
- `BED` – targeted amplicon BED file

---

## R Dependencies

Script `17_haplo_stats.r` requires the following R packages:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "openxlsx"))
install.packages("haplo.stats")
```

---

## Notes

- All scripts are written to be run in order (numbered `01` → `17`). N3`.
- Scripts assume a directory structure of `<root>/<cohort>/<tissue_type>/<sample>/`.
- The pipeline was developed and tested on **Ubuntu 24 (WSL2)** with Python 3.11.
- Ion Torrent variant calling in `05_variant_calling.sh` wraps `tvc` (Torrent Variant Caller); ensure this is available and in `PATH` if using a non-standard installation.

---

## Author

Research conducted as part of a Master's Thesis (TFM), February 2026.
