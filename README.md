# GSDMB Targeted Sequencing Analysis Pipeline

A reproducible analysis pipeline for targeted DNA sequencing of the **GSDMB locus and surrounding Chr17 panel**, including QC, variant calling, functional annotation, descriptive visualisation, SNP enrichment testing, clinical association analysis, and haplotype-based modelling.

Developed as part of a Master's thesis (TFM), this repository implements an end-to-end workflow spanning sequencing quality control, variant annotation, haplotype analysis, clinical harmonisation, and SNP-phenotype association modelling.

## What This Repository Does

The pipeline supports the following stages:

- technical validation of BAM files and indexing
- sequencing quality control and zero-coverage assessment
- targeted variant calling and filtering
- functional annotation with Ensembl VEP and pathogenicity predictors
- variant statistics and enrichment testing
- haplotype phasing and haplotype-based association testing
- clinical data integration and harmonisation
- SNP-clinical association analysis
- haplotype-clinical association analysis
- survival analysis and interactive visualisation outputs

## Repository Structure

```text
01_check_and_index.sh            BAM integrity checks and indexing
02_dna_qc.sh                     DNA sequencing QC pipeline
02b_more-qc.sh                   PASS manifest creation and zero-coverage summaries
02c_zerocovinfo.py               Zero-coverage Excel reporting
03_qc_visualisation.py           QC figures and panel-level visual summaries
04_technical_audit.py            Technical audit of amplicon performance
05_variant_calling.sh            Targeted variant calling workflow
06_annotation.sh                 Functional annotation with VEP
07_merge_annotations.py          Merge annotated VCFs into a canonical table
07b_variant_qc.py                Variant-level QC checks
08_mapping.py                    Global variant landscape plot
09_gsdmb_only.py                 GSDMB-only landscape plot
09b_landscape_comparison.py       Normal-versus-tumour follow-up for scripts 8 and 9
10_variant_stats.py              Descriptive consequence/impact summaries
11_SNPs.py                       Common SNP identification and benchmarking
12_stats_enrichment.py           Tumour-versus-control SNP enrichment testing
13_permutation_testing.py        Permutation-based SNP comparison
14_interactive_dashboard.py      Interactive Plotly dashboard
15_haplotypes.sh                 Haplotype phasing workflow
15_haplotype_stats.R             Haplotype statistical analysis
16_excel_harmonisation.py       Clinical-data harmonisation
17_snp_association.py            SNP-clinical association analysis
18_haplotype_association.py      Haplotype-clinical association analysis
19_1000g_haplotype_comparison.py Haplotype comparison against phased 1000 Genomes reference populations
pipeline_config.toml             Shared paths, thresholds, and grouping choices
pipeline_utils.py                Shared helper functions
pipeline_validation.py           Shared validation helpers
association_runtime.py           Shared runtime defaults for scripts 17, 18, and 19
environment.yml                  Shared conda/micromamba environment specification
requirements.txt                 Shared pip specification for the legacy tfm_env workflow
docs/notes/PIPELINE.md           Short pipeline notes and doc pointers
```

## Quick Start

If you want the shortest practical version of the workflow on this machine, use:

```bash
cd /home/gadeaalonsoj/tfm
./run_full_pipeline.sh
```

The launcher now handles the expected environment split automatically:

- shell-based BAM/QC/calling/phasing steps use `bam-steps`
- Python analysis steps prefer the active virtual environment, then `/home/gadeaalonsoj/tfm_env/bin/python`, then a repo-local `tfm_env` if present
- VEP annotation uses the repo-local `miniconda3/envs/vep_env/bin` when available, otherwise falls back to another `vep_env` on `conda`/`micromamba`, and finally to the local `ensembl-vep` checkout
- the R haplotype statistics step uses `Rscript` and the script now bootstraps missing R packages such as `haplo.stats` when needed

If you prefer to run by hand instead of using the launcher, the expected environments are documented below.

## Recommended Run Order

1. `01_check_and_index.sh`
2. `02_dna_qc.sh`
3. `02b_more-qc.sh`
4. `02c_zerocovinfo.py`
5. `03_qc_visualisation.py`
6. `04_technical_audit.py`
7. `05_variant_calling.sh`
8. `06_annotation.sh`
9. `07_merge_annotations.py`
10. `07b_variant_qc.py`
11. `08_mapping.py`
12. `09_gsdmb_only.py`
13. `09b_landscape_comparison.py`
14. `10_variant_stats.py`
15. `11_SNPs.py`
16. `12_stats_enrichment.py`
17. `13_permutation_testing.py`
18. `14_interactive_dashboard.py`
19. `15_haplotypes.sh`
20. `15_haplotype_stats.R`
21. `16_excel_harmonisation.py`
22. `17_snp_association.py`
23. `18_haplotype_association.py`
24. `19_1000g_haplotype_comparison.py` when a phased 1000 Genomes VCF and panel file are available

## Computational Environments

The original workflow used **three separate environments** because bioinformatics processing, annotation, and statistical analysis rely on different toolchains. On this analysis machine, those environments are not all located in the same place, so the actual paths are documented here.

### 1. Bioinformatics Environment: `bam-steps` (micromamba)

Used for shell-based genomic processing:

- `01_check_and_index.sh`
- `02_dna_qc.sh`
- `02b_more-qc.sh`
- `05_variant_calling.sh`
- `15_haplotypes.sh`

Create and activate it:

```bash
micromamba create -n bam-steps samtools bcftools mosdepth -c bioconda -c conda-forge
micromamba activate bam-steps
```

Deactivate when finished:

```bash
micromamba deactivate
```

Core tools in this environment:

- `samtools` for BAM integrity checking and indexing
- `bcftools` for VCF normalisation, filtering, and merging
- `mosdepth` for per-base and per-region coverage metrics

### 2. Python Analysis Environment: `tfm_env` (virtual environment)

Used for the Python-based analytical part of the pipeline:

- `02c_zerocovinfo.py`
- `03_qc_visualisation.py`
- `04_technical_audit.py`
- `07_merge_annotations.py`
- `07b_variant_qc.py`
- `08_mapping.py`
- `09_gsdmb_only.py`
- `09b_landscape_comparison.py`
- `10_variant_stats.py`
- `11_SNPs.py`
- `12_stats_enrichment.py`
- `13_permutation_testing.py`
- `14_interactive_dashboard.py`
- `16_excel_harmonisation.py`
- `17_snp_association.py`
- `18_haplotype_association.py`
- `19_1000g_haplotype_comparison.py`

Recommended creation on this machine:

```bash
cd /home/gadeaalonsoj
python3 -m venv tfm_env
source tfm_env/bin/activate
pip install -r /home/gadeaalonsoj/tfm/requirements.txt
```

Activate it later with:

```bash
source /home/gadeaalonsoj/tfm_env/bin/activate
```

Deactivate it with:

```bash
deactivate
```

Note:

- the launcher prefers an already-active virtual environment first
- if no virtual environment is active, it then prefers `/home/gadeaalonsoj/tfm_env/bin/python`
- a repo-local `tfm_env` is only used as a fallback if present and runnable

### 3. VEP Annotation Runtime: repo-local `vep_env` and `ensembl-vep`

Used specifically for:

- `06_annotation.sh`

On this machine, the annotation step is set up as a **repo-local runtime** rather than a global `conda` installation. The launcher currently prefers:

1. `/home/gadeaalonsoj/tfm/miniconda3/envs/vep_env/bin`
2. another `vep_env` found via `conda` or `micromamba`
3. the local `/home/gadeaalonsoj/tfm/ensembl-vep` checkout

The following repo-local components were found in this workspace:

- `/home/gadeaalonsoj/tfm/miniconda3/envs/vep_env`
- `/home/gadeaalonsoj/tfm/ensembl-vep`

For manual execution without the launcher, the safest lightweight setup is:

```bash
export PATH="/home/gadeaalonsoj/tfm/miniconda3/envs/vep_env/bin:$PATH"
cd /home/gadeaalonsoj/tfm
bash 06_annotation.sh
```

Important note about creation:

- this README documents how `bam-steps` and `tfm_env` are created from scratch
- the exact original creation command for the repo-local `vep_env` is **not currently preserved** in this repository
- in practice, this workspace uses the already-provisioned repo-local `vep_env` together with the local VEP resources under `.vep/` and `ensembl-vep/`

This annotation runtime is expected to provide:

- Ensembl VEP with the GRCh38 cache
- GRCh38 reference FASTA
- CADD resources
- dbNSFP resources
- any local annotation assets required by the script

Note: `06_annotation.sh` unsets `PERL5LIB` and `PERL_LOCAL_LIB_ROOT` automatically to avoid WSL Perl conflicts.

## Practical Execution by Stage

### Stage 1. Sequencing QC

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

Typical outputs:

- coverage metrics
- zero-coverage region reports
- PCA and QC plots
- technical audit summaries

### Stage 2. Variant Calling and Annotation

```bash
micromamba activate bam-steps
./05_variant_calling.sh

conda activate vep_env
./06_annotation.sh
```

Annotation resources used by the workflow include:

- gnomAD allele frequencies
- CADD
- SIFT
- PolyPhen-2
- REVEL
- DANN
- COSMIC context where available in the annotation resources

### Stage 3. Variant Processing and Enrichment

```bash
source tfm_env/bin/activate
python 07_merge_annotations.py
python 07b_variant_qc.py
python 08_mapping.py
python 09_gsdmb_only.py
python 10_variant_stats.py
python 11_SNPs.py
python 12_stats_enrichment.py
python 13_permutation_testing.py
python 14_interactive_dashboard.py
```

This stage produces:

- consolidated annotated variant tables
- variant QC summaries
- population-frequency filtered SNP summaries
- enrichment and permutation results
- interactive dashboards and descriptive plots

### Stage 4. Haplotype Analysis

```bash
micromamba activate bam-steps
./15_haplotypes.sh

Rscript 15_haplotype_stats.R
```

This stage covers:

- BEAGLE phasing
- haplotype frequency estimation
- haplotype association modelling

### Stage 5. Clinical Harmonisation and Association Modelling

```bash
source tfm_env/bin/activate
python 16_excel_harmonisation.py
python 17_snp_association.py
python 18_haplotype_association.py
python 19_1000g_haplotype_comparison.py --vcf /path/to/1000G_phased.vcf.gz --panel /path/to/1000G.panel
```

## Configuration

Shared paths and thresholds are defined in [pipeline_config.toml](pipeline_config.toml).

Important points:

- the new `19_1000g_haplotype_comparison.py` step writes to the configured `thousand_genomes_results` directory by default
- the external 1000 Genomes phased VCF and population panel remain command-line inputs because their local paths depend on where you store the IGSR data

- the current config uses absolute paths from the original analysis machine
- update these paths before reusing the pipeline elsewhere
- pooled-control behaviour is also defined there

Current grouping logic includes:

- pooled controls = `Breast normal + Endometrium normal`
- cohort-specific comparisons = tumour cohort vs pooled controls
- global comparisons = all tumours vs all controls
- cancer-risk comparisons in scripts `17` and `18` also use the same pooled healthy control arm

## Outputs

Major outputs generated by the pipeline include:

- QC tables and QC summary figures
- annotated Excel reports
- variant landscape plots
- SNP benchmarking and enrichment results
- permutation-test summaries
- interactive HTML dashboards
- harmonised clinical master sheets
- SNP-clinical association workbooks and figures
- phased haplotype tables and haplotype association results
- threshold-specific 1000 Genomes LD-block overlap tables and haplotype concordance summaries

Most analysis outputs are written to the configured results directory, currently:

```text
/home/gadeaalonsoj/tfm/gsdmb_final_results
```

## Reproducibility Notes

The repository now includes a shared reproducibility layer:

- [pipeline_config.toml](pipeline_config.toml) centralises paths and thresholds
- [pipeline_utils.py](pipeline_utils.py) centralises shared transformations
- [pipeline_validation.py](pipeline_validation.py) centralises common integrity checks
- [association_runtime.py](association_runtime.py) keeps scripts 17 and 18 aligned
- [environment.yml](environment.yml) provides a shared conda/micromamba environment definition
- [requirements.txt](requirements.txt) supports the legacy `tfm_env` virtual-environment workflow

Additional details are described in [REPRODUCIBILITY.md](docs/notes/REPRODUCIBILITY.md).

## Data Availability

This repository contains code and workflow logic. Raw sequencing data, patient-level clinical data, and any restricted metadata should only be shared if ethics approval, institutional policy, and privacy requirements allow it.

If the repository is made public, it is usually best to:

- exclude patient-identifiable data
- exclude large derived result folders unless they are intentionally released
- provide only de-identified example inputs where appropriate

## Use and Reuse

This repository is intended to document the analytical workflow clearly enough for academic review, reproducibility, and supervised reuse. When reused on another machine or dataset, update local paths in `pipeline_config.toml` and verify the required external bioinformatics tools are available.
