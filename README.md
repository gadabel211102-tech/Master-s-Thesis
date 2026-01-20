# GSDMB Targeted Sequencing Analysis Pipeline

##  Project Overview

This research investigates the clinical significance of **Gasdermin-B (GSDMB)** variants (SNPs and isoforms) in the context of breast and endometrial cancers.

The repository provides a modular bioinformatics pipeline for Quality Control (QC), coordinate harmonisation across different reference headers (e.g., `17` vs `chr17`), and high-resolution cohort-wide visualisation. It is specifically designed to ensure data interoperability across diverse clinical datasets.

---

##  Pipeline Structure & Logic

### 1. Data Validation & Indexing (`01_check_and_index.sh`)

* **Purpose:** Ensures raw data integrity before long-running steps.
* **Logic:** Uses `samtools quickcheck` to find corrupted BAMs and `samtools index` to generate missing `.bai` files.

### 2. Harmonised DNA QC (`02_dna_qc.sh`)

* **Purpose:** The primary technical gatekeeper.
* **Coordinate Harmonisation:** Automatically detects if the reference uses `17` or `chr17` and normalises the output.
* **Tools:** Uses `mosdepth` for rapid depth calculation.
* **Sorting:** Automatically segregates samples into `pass_bams` and `fail_bams` based on the thresholds below.

### 3. Cohort Analysis & Visualisation (`03_qc_analysis.py`)

* **Purpose:** Consolidates all samples into a single research dashboard.
* **Key Figures:** Generates **13 figures**, including PCA of technical variance, log10 coverage heatmaps (red-highlighted failing samples), and systemic failure reports using genomic coordinates.
* **Standards:** Uses British spellings (*normalise*, *tumour*, *colour*) for all clinical labeling.

### 4. Variant Calling (`04_variant_calling.sh`)

* **Purpose:** Clinical-grade mutation identification.
* **Engine:** Ion Torrent Variant Caller (TVC).
* **Post-Processing:** Uses `bcftools` for normalisation (splitting multi-allelic sites) and hard-filtering variants based on Gema's validated thresholds.

### 5. Final Technical Audit (`05_technical_audit.py`)

* **Purpose:** A "Last-Step QC" to identify systemic failures and annotate the final matrices.
* **Logic:** Cross-references the `Full_Panel.csv` matrices with the original BED file to map coordinates to **RS IDs** and **Gene Names**.
* **Audit:** Automatically flags the **Top 10 Worst Amplicons** for each cohort (based on median depth) to identify potential technical dropouts or "blind spots" in the panel.
* **Output:** Consolidates all annotations and audit results into a multi-sheet Excel workbook (`GSDMB_Annotated_Report_Fixed.xlsx`).

---

##  Validated Thresholds

### DNA Quality Control (QC)

| Metric | Threshold |
| --- | --- |
| **Total Reads** | ≥ 120,000 |
| **Mapped Reads** | ≥ 90.0% |
| **On-Target Rate** | ≥ 75.0% |
| **Mean Coverage** | ≥ 200× |
| **Uniformity (100x)** | ≥ 95.0% |
| **Amplicon Floor** | ≥ 50× |

### Variant Calling (TVC)

| Parameter | Final Clinical Threshold |
| --- | --- |
| **Depth per Variant (DP)** | **≥ 50x** |
| **Quality (QUAL)** | ≥ 20 |
| **TVC Filter** | Only `FILTER=PASS` |
| **Homopolymers** | Exclude if > 8 bases |
| **SNP vs Indels** | Unified filtering approach |

---

##  Installation & Usage

### Environment Setup

```bash
# 1. Bioinformatics (Micromamba)
micromamba create -n bam-steps samtools bcftools mosdepth -c bioconda -c conda-forge

# 2. Visualisation (Python venv)
python3 -m venv tfm_env
source tfm_env/bin/activate
pip install pandas seaborn matplotlib scikit-learn openpyxl

```

### Full Execution Sequence

```bash
# Step 1: Check Index
./01_check_and_index.sh -i ./raw_data/bams

# Step 2: Run QC
./02_dna_qc.sh -i ./raw_data/bams -o ./results -r hg38.fa -b GSDMB_targets.bed

# Step 3: Analysis & Plots
python3 03_qc_analysis.py

# Step 4: Call Variants
./04_variant_calling.sh -m ./results/qc_pass_manifest.txt -r hg38.fa -b GSDMB_targets.bed -o ./variants

# Step 5: Final Technical Audit
python3 05_technical_audit.py
```



