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

* **Purpose:** Acts as the primary technical gatekeeper for the pipeline.
* **Coordinate Harmonisation:** Automatically detects chromosome naming conventions and normalises BED coordinates (e.g., ensuring `17` vs `chr17` consistency) to ensure compatibility across genomic resources.
* **Tools:** Utilises `mosdepth` for rapid depth calculation and `samtools` for alignment statistics.
* **Automated Sorting:** Segregates samples into `pass_bams/` and `fail_bams/` directories based on customisable technical thresholds such as mean coverage, mapping rate, and uniformity.
* **Detailed Reporting:** Generates a comprehensive `qc_summary.tsv` which identifies the "worst-performing amplicon" for every sample to assist with troubleshooting.

---

### 3. Post-QC Analysis & Zero-Coverage Profiling (`02.5_more-qc.sh`)

* **Purpose:** Provides a bridge between raw QC and downstream variant calling by streamlining sample management and identifying "blind spots" in the assay.
* **Manifest Generation:** Automatically parses the QC summary to create a validated list of paths for only those samples that passed technical checks, preventing "failed" data from entering the analysis stream.
* **Zero-Coverage Mapping:** Performs a high-resolution, per-base analysis across every target region to count exactly how many bases within an amplicon have zero coverage.
* **Cohort-Level Comparison:** Compiles a master matrix (the "Zero-Cov" report) that aligns samples side-by-side, making it easy to spot systemic amplicon failures versus sample-specific issues.
* **Optimisation:** Uses a "compute-once" strategy with temporary `mosdepth` files to ensure the per-base analysis is efficient even for large cohorts.

### 3. Cohort Analysis & Visualisation (`03_qc_visualisation.py`)

* **Purpose:** Consolidates all samples into a single research dashboard.
* **Key Figures:** Generates **13 figures**, including PCA of technical variance, log10 coverage heatmaps (red-highlighted failing samples), and systemic failure reports using genomic coordinates.
* **Standards:** Uses British spellings (*normalise*, *tumour*, *colour*) for all clinical labeling.

### 4. Unified Amplicon Technical Audit (`04_technical_audit.py`)

* **Purpose:** Acts as a dual-mode diagnostic tool that combines deep sequence analysis with cohort-level performance reporting.
* **Operational Modes:**
    * **Comprehensive Mode:** Performs a full technical audit by correlating sequence complexity (GC content, homopolymers) with actual coverage data.
    * **Quick Mode:** Generates rapid Excel-based reports highlighting the "Top 10 Worst Amplicons" per cohort for immediate review.
* **Sequence Complexity Profiling:** Calculates critical metrics for Ion Torrent chemistry, including GC skew, sliding-window GC extremes, and the identification of self-complementary sequences (hairpin risk).
* **Technical "Red Flag" Detection:** Automatically scans for long homopolymer runs and Simple Sequence Repeats (SSRs) that often lead to signal droop or alignment errors.
* **Automated Normalisation:** Internally handles chromosome naming inconsistencies (e.g., converting '17' to 'chr17') to ensure seamless data integration across different genomic builds.

### 5. Variant Calling (`05_variant_calling.sh`)

* **Purpose:** Clinical-grade mutation identification.
* **Engine:** Ion Torrent Variant Caller (TVC).
* **Post-Processing:** Uses `bcftools` for normalisation (splitting multi-allelic sites) and hard-filtering variants based on validated thresholds.

### 6. Variant Annotation (`06_annotation.sh`)

* **Purpose:** Adds biological and functional context to identified variants.
* **Engine:** **Ensembl Variant Effect Predictor (VEP)** using the GRCh38 assembly.
* **Key Features:** * Configured for **offline high-speed processing** using local cache directories.
    * Includes environment-specific "Secret Sauce" fixes for Perl library paths to ensure stability in WSL environments.
* **Data Extraction:** Automatically parses complex VEP headers to extract Chromosome, Position, Gene Symbol, Consequence, and Protein Change (HGVSp).

### 7. Clinical Report Integration (`07_merge-annotations.py`)

* **Purpose:** Consolidates technical QC data and biological variants into a single, researcher-friendly Excel report.
* **Core Logic:**
    * **Impact Prioritisation:** Automatically ranks variants by biological severity (High > Moderate > Low > Modifier).
    * **Robust Parsing:** Uses case-insensitive column matching and dynamic VCF header extraction to handle diverse annotation outputs.
* **Final Deliverable:** Appends a `Biological_Annotations` sheet to the master project report (`GSDMB_Annotated_Report_Fixed.xlsx`).

### 8. Multi-Gene Mutational Landscape (`08_mapping.py`)

* **Purpose:** Generates a "clean" global map of all identified variants across the targeted genomic regions.
* **Logic:** Visualises variant frequency across genomic positions (Chr17 Mb) for all four cohorts (Breast/Endometrium, Tumour/Normal) simultaneously.
* **Key Features:** * **Impact Encoding:** Uses distinct markers (Circle for HIGH, Square for MODIFIER) and high-contrast colours to differentiate between gene symbols.
    * **Cohort Faceting:** Automatically splits views to allow direct technical and biological comparison between tissue types.

### 9. GSDMB-Exclusive Analysis (`09_gsdmb-only.py`)

* **Purpose:** Provides a high-resolution, focused view of mutations specifically within the *GSDMB* gene boundaries.
* **Logic:** Filters the master dataset for *GSDMB* symbols and switches the X-axis to base-pair (bp) precision for exact exon-level localisation.
* **Key Features:** Utilises a "flare" colour palette to highlight high-impact areas, providing a cleaner visual for identifying specific *GSDMB* haplotypes.

### 10. Variant Statistics & Distribution (`10_variant-stats.py`)

* **Purpose:** Quantifies the mutational burden and classifies variants by type and impact for statistical and clinical reporting.
* **Logic:** Generates both log-scaled and linear bar plots to compare "Total Records" (every occurrence) versus "Unique Variants" across tissue groups.
* **Key Features:** * **Automated Categorisation:** Groups findings into four essential clinical plots (Total/Unique x Log/Linear).
    * **Tabular Exports:** Produces raw CSV tables detailing precise counts for every **Impact** level (HIGH/MODERATE/LOW) and **Consequence** type (Missense/Synonymous).

### 11. Global SNP Discovery & gnomAD Comparison (`11_SNPs.py`)

* **Purpose:** Benchmarks study variants against the **gnomAD Non-Finnish European (NFE)** population to distinguish common polymorphisms from rare or novel mutations.
* **Logic:** Extracts rsIDs and correlates cohort Allele Frequencies (AF) with ancestry-matched population data, ensuring novel variants (missing from gnomAD) are preserved in the analysis.
* **Key Features:**
    * **Inclusive 1% Filter:** Automatically identifies common SNPs (>1% NFE frequency) while explicitly retaining novel/rare variants for further investigation.
    * **Identity Correlation:** Generates scatter plots with an "Identity Line" to visualise deviations between study frequencies and global European averages.
    * **GSDMB Targeting:** Specifically highlights *GSDMB* variants with black target rings to distinguish them from other panel targets.
    * **Automated Reporting:** Produces a comparative Excel catalogue partitioned by cohort, tissue type, and gene-specific subsets.
 
### 12. Statistical Enrichment Analysis (`12_stats-enrichment.py`)

* **Purpose:** Performs comparative statistical testing between tissue types to identify variants significantly enriched in tumour versus normal samples.
* **Logic:** Utilises the filtered SNP list from Script 11 to construct 2x2 contingency tables and applies **Fisher’s Exact Test** to calculate the significance (p-value) and Odds Ratio for each variant.
* **Key Features:**
    * **Fisher's Exact Testing:** Robustly compares variant prevalence across cohort groups, even with small sample sizes.
    * **Enrichment Visualisation:** Generates **Volcano Plots** mapping the Odds Ratio against $-log10(p-value)$ to highlight statistically significant outcomes.
    * **Standardised Metrics:** Calculates Tumour vs. Normal frequencies (%) and log-scaled Odds Ratios to identify somatic-like enrichment patterns in germline-relevant SNPs.
    * **Targeted Results:** Exports a detailed statistical report (Excel) sorted by significance, allowing for rapid identification of top-tier candidates for further functional studies.
 

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

# 3. VEP
conda activate vep_env

Rscript 19b_haplo_stats.R

```

### Full Execution Sequence

```bash
# Step 1: Check Index
./01_check_and_index.sh -i ./raw_data/bams

# Step 2: Run QC and sort samples
./02_dna_qc.sh -i ./raw_data/bams -o ./results/dna_qc -r hg38.fa -b GSDMB_targets.bed -t 8

# Step 2.5: Generate manifest of PASS samples and map zero-coverage bases
./02.5_more-qc.sh \
    -q ./results/dna_qc \
    -m ./results/manifests/pass_samples.txt \
    -z ./results/zero_coverage_analysis \
    -t 8

# Step 3: Analysis & Plots
python3 03_qc_visualisation.py

# Step 4: Final Technical Audit
python3 04_technical_audit.py both

# Step 5: Call Variants
./05_variant_calling.sh -m ./results/qc_pass_manifest.txt -r hg38.fa -b GSDMB_targets.bed -o ./variants

# Step 6: Annotate Variants
./06_annotation.sh

# Step 7: Final Report Consolidation
python3 07_merge-annotations.py

# Step 8: Generate global mutational landscape
python3 08_mapping.py

# Step 9: Generate GSDMB-specific mutation map
python3 09_gsdmb-only.py

# Step 10: Calculate variant statistics and generate clinical bar plots
python3 10_variant-stats.py

# Step 11: Compare study SNPs with global population data (gnomAD)
python3 11_SNPs.py

# Step 12: Statistical Enrichment and Tumour vs. Normal Comparison
python3 12_stats-enrichment.py


```



