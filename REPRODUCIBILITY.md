# Reproducibility Guide

This document explains how to rerun the reorganised workflow reproducibly.

## 1. Launcher profiles

The main launcher is now profile-aware.

### Core rerun

```bash
cd /path/to/tfm
./run_full_pipeline.sh
```

Equivalent explicit form:

```bash
./run_full_pipeline.sh --profile core
```

### Exploratory / validation rerun

```bash
./run_full_pipeline.sh --profile exploratory
```

### Utility / reporting rerun

```bash
./run_full_pipeline.sh --profile utility
```

### Full rerun

```bash
./run_full_pipeline.sh --profile full
```

### Bounded reruns

```bash
./run_full_pipeline.sh --profile core --from 11 --to 18
./run_full_pipeline.sh --profile exploratory --from 24 --to 26
```

## 2. Environment strategy

The original three-runtime split still applies.

### `bam-steps`

Used for:

- `01_check_and_index.sh`
- `02_dna_qc.sh`
- `02b_more-qc.sh`
- `05_variant_calling.sh`
- `15_haplotypes.sh`

### `tfm_env` or another explicit Python interpreter

Used for the Python-based core, exploratory, and utility scripts.

### `vep_env` or the repo-local VEP runtime

Used for:

- `06_annotation.sh`

## 3. Core reproducibility boundary

The canonical reproducible thesis workflow is:

1. DNA BAM validation and QC
2. Variant calling and annotation
3. Canonical SNP backbone generation and primary SNP enrichment
4. Haplotype phasing and principal haplotype statistics
5. Clinical master harmonisation
6. SNP-clinical and haplotype-clinical modelling
7. RNA QC gating
8. Workbook-first isoform analysis
9. BAM-derived panel-gene quantification
10. Final RNA/SNP/haplotype integration

Scripts outside that chain are available but should be described as optional in methods and supervisor-facing workflow summaries.

## 4. Important fallback in stage 23

Stage `23` no longer requires stage `22` as a hard prerequisite.

- if `analysis_results/22_haplotype_first_interpretation/22_Haplotype_First_Interpretation.xlsx` exists, stage `23` uses it
- otherwise stage `23` derives significant haplotypes directly from the stage-15 association workbook

This keeps the core rerun self-contained.

## 5. External-data caveats

These steps are optional and depend on local external resources:

- stage `19` requires the phased 1000 Genomes VCF and panel
- stage `24` requires the phased 1000 Genomes VCF and panel and applies female-only filtering when panel sex metadata are available
- stage `24` cannot apply age matching because the local 1000 Genomes panel does not contain age metadata

## 6. Minimal verification

A lightweight syntax check after editing is:

```bash
bash -n run_full_pipeline.sh .run_full_pipeline_clean.sh 01_check_and_index.sh 02_dna_qc.sh 02b_more-qc.sh 05_variant_calling.sh 06_annotation.sh 15_haplotypes.sh
python -m py_compile 07_merge_annotations.py 11_SNPs.py 12_stats_enrichment.py 16_excel_harmonisation.py 17_snp_association.py 18_haplotype_association.py 19b_rna_qc.py 20_objective2_isoform_expression_association.py 20b_panel_gene_expression.py 23_rna_integration.py 24_external_cross_validation.py 25_rs11078928_rs869402_haplotype_focus.py 26_snp_functional_interpretation.py association_runtime.py objective2_gene_report_exports.py objective2_rna_qc_utils.py pipeline_utils.py pipeline_validation.py figure_style.py
Rscript -e "parse(file='15_haplotype_stats.R')"
```

That does not replace full data reruns, but it catches syntax drift quickly.
