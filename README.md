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

### 2.1 Manifest Generation (´02.1_manifest.sh´)
* **Purpose:** Automates the bridge between quality control and downstream analysis by filtering validated samples.
* **Logic:** Parses the qc_summary.tsv generated in the previous step to identify samples with a PASS status. It cross-references these names against the physical files in the pass_bams directory to ensure data integrity.
* **Validation:** Employs a robust awk-based header detection to handle column mapping dynamically, ensuring the script remains resilient even if the summary file structure changes.
* **Output:** Generates a flat-file manifest containing absolute paths to high-quality BAM files, which serves as the definitive input list for cohort-level analysis or variant calling.

### 2.2 Calculating how many bases had 0 coverage per amplicon (´02.2_zero_base-coverage.sh´)
* **Purpose:** Extends the DNA‑QC workflow by measuring fine‑grained coverage completeness across all amplicons. Whereas mean depth and uniformity reflect overall performance, this stage explicitly quantifies how many bases within each target region received zero sequencing coverage per sample, per cohort.
* **Logic:** Iterates over the cleaned BED coordinates and re‑profiles each amplicon using temporary per‑base mosdepth output. This enables precise counting of uncovered bases, overcoming the limitation of region‑summary depth files. For each sample, the script constructs a coverage vector and computes the number of positions with depth 0 inside each amplicon.
* **Robustness:** Works seamlessly on both PASS and FAIL samples, ensuring that low‑quality sequencing runs are accurately represented. Its design avoids rewriting any core QC scripts and leaves all existing outputs untouched.
* **Output:** Produces a tabular report for each cohort where rows represent amplicons and columns represent samples. Each table entry indicates the number of uncovered bases for that amplicon in that sample, providing a granular measure of panel completeness suitable for downstream reporting or troubleshooting.

### 3. Cohort Analysis & Visualisation (`03_qc_analysis.py`)

* **Purpose:** Consolidates all samples into a single research dashboard.
* **Key Figures:** Generates **13 figures**, including PCA of technical variance, log10 coverage heatmaps (red-highlighted failing samples), and systemic failure reports using genomic coordinates.
* **Standards:** Uses British spellings (*normalise*, *tumour*, *colour*) for all clinical labeling.

### 4. Final Technical Audit (`04_technical_audit.py`)

* **Purpose:** A "Last-Step QC" to identify systemic failures and annotate the final matrices.
* **Logic:** Cross-references the `Full_Panel.csv` matrices with the original BED file to map coordinates to **RS IDs** and **Gene Names**.
* **Audit:** Automatically flags the **Top 10 Worst Amplicons** for each cohort (based on median depth) to identify potential technical dropouts or "blind spots" in the panel.
* **Output:** Consolidates all annotations and audit results into a multi-sheet Excel workbook (`GSDMB_Annotated_Report_Fixed.xlsx`).

### 5. Variant Calling (`05_variant_calling.sh`)

* **Purpose:** Clinical-grade mutation identification.
* **Engine:** Ion Torrent Variant Caller (TVC).
* **Post-Processing:** Uses `bcftools` for normalisation (splitting multi-allelic sites) and hard-filtering variants based on validated thresholds.

### 6. Variant Annotation ([`06_annotation.sh`](06_annotation.sh))

* **Purpose:** Adds biological and functional context to identified variants.
* **Engine:** **Ensembl Variant Effect Predictor (VEP)** using the GRCh38 assembly.
* **Key Features:** * Configured for **offline high-speed processing** using local cache directories.
    * Includes environment-specific "Secret Sauce" fixes for Perl library paths to ensure stability in WSL environments.
* **Data Extraction:** Automatically parses complex VEP headers to extract Chromosome, Position, Gene Symbol, Consequence, and Protein Change (HGVSp).

### 7. Clinical Report Integration ([`07_merge-annotations.py`](07_merge-annotations.py))

* **Purpose:** Consolidates technical QC data and biological variants into a single, researcher-friendly Excel report.
* **Core Logic:**
    * **Impact Prioritisation:** Automatically ranks variants by biological severity (High > Moderate > Low > Modifier).
    * **Robust Parsing:** Uses case-insensitive column matching and dynamic VCF header extraction to handle diverse annotation outputs.
* **Final Deliverable:** Appends a `Biological_Annotations` sheet to the master project report (`GSDMB_Annotated_Report_Fixed.xlsx`).
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
micromamba activate bam-steps

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

# Step 2.1: Make PASS manifests
 ./02.5_manifest.sh /path/to/cohort/dna_qc /path/to/output_manifest.txt

# Step 2.2: Zero base coverage
02.2_zero_base_coverage.sh <ROOT_DIR> <OUT_DIR>

# Step 3: Analysis & Plots
python3 03_qc_analysis.py

# Step 4: Final Technical Audit
python3 04_technical_audit.py

# Step 5: Call Variants
./05_variant_calling.sh -m ./results/qc_pass_manifest.txt -r hg38.fa -b GSDMB_targets.bed -o ./variants

# Step 6: Annotate Variants
./06_annotation.sh

# Step 7: Final Report Consolidation
python3 07_merge-annotations.py


```



