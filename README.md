# GSDMB Targeted Sequencing Analysis Pipeline

This repository contains a modular bioinformatics pipeline for the Quality Control (QC), cross-cohort visualisation, and variant calling of targeted DNA sequencing data. The pipeline is specifically optimised for the **GSDMB panel** and handles multiple tissue types (Breast and Endometrium) by harmonising genomic coordinates across different BAM reference headers.

## 🧬 Pipeline Overview

The workflow is divided into four main stages, moving from raw alignment files to finalised clinical-grade visualisations and variant calls.

### 1. Pre-processing (`01_check_and_index.sh`)

Performs an initial integrity check on raw BAM files using `samtools quickcheck` and ensures all files have a valid `.bai` index for downstream processing.

### 2. Harmonised DNA QC (`02_dna_qc.sh`)

This script executes the primary technical audit of the sequencing run.

* **Coordinate Harmonisation**: Automatically detects and forces the `chr17` prefix for all amplicons, ensuring samples from different pipelines (e.g., using `17` vs `chr17`) can be compared.
* **Metrics Tracked**: Total reads, mapping percentage, on-target percentage, mean coverage, and uniformity ().
* **Automated Sorting**: Segregates samples into `pass_bams` and `fail_bams` based on validated thresholds.

### 3. Cohort Analysis & Visualisation (`03_qc_analysis.py`)

A Python-based suite that consolidates data across cohorts to produce 16 distinct figures for your TFM.

* **Key Visuals**: PCA of technical features, Log-scaled Heatmaps (Global and Cohort-specific), Panel Landscape Audits, and Systematic Amplicon Failure reports.
* **Statistical Rigour**: Implements StandardScaler and PCA to map the technical variance between Breast and Endometrial samples.

### 4. Variant Calling & Filtering (`04_variant_calling.sh`)

Utilises the Ion Torrent Variant Caller (TVC) or standard callers to identify mutations within the target regions.

* **Hard Filtering**: Applies stringency filters for Depth () and Quality ().
* **Cohort Merging**: Merges individual VCFs into a multi-sample cohort file for manual review.

---

## 🚀 Getting Started

### Prerequisites

* **Bioinformatics Tools**: `samtools`, `bcftools`, `mosdepth`.
* **Python Environment**: `pandas`, `seaborn`, `matplotlib`, `scikit-learn`.

### Execution Sequence

1. **Prepare BAMs**:
```bash
./01_check_and_index.sh -i /path/to/bams

```


2. **Run Harmonised QC**:
```bash
./02_dna_qc.sh -i ./raw_data -o ./qc_results -r genome.fa -b GSDMB_targets.bed

```


3. **Generate Visualisations**:
```bash
python3 03_qc_analysis.py

```


4. **Call Variants**:
```bash
./04_variant_calling.sh -m pass_manifest.txt -r genome.fa -b targets.bed -o ./variants

```



---

## 📊 Validated QC Thresholds

The pipeline employs the following strict criteria for diagnostic-grade analysis:

| Metric | Threshold |
| --- | --- |
| **Total Reads** |  120k |
| **Mapped Percentage** | 90% |
| **On-Target Percentage** | 75% |
| **Mean Coverage** | 200x |
| **Uniformity ()** | 95% |
| **Amplicon Floor** | 50x |
