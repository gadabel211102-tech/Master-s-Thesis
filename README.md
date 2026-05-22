# Reproducible Cancer Genomics Workflow for `GSDMB` and `17q12-q21`

This repository contains the analysis code for my Master's thesis on the
clinical significance of polymorphic variation around `GSDMB` and the
`17q12-q21` region in breast and endometrial cancer.

## What this repository demonstrates

- End-to-end biomedical data analysis across Bash, Python, and R
- Workflow orchestration for a multi-stage genomics pipeline
- Variant interpretation around a focused disease-relevant locus
- Integration of sequencing, haplotype, clinical, and RNA-derived results
- Careful documentation of analytical boundaries and non-causal interpretation

## Privacy and data availability

This public repository does not include raw or participant-level study data.

- No raw sequencing files such as BAMs, VCFs, or full clinical source workbooks
- No patient-identifiable outputs or participant-level result exports
- No generated `analysis_results/`, local runtime environments, or machine-local
  reference bundles
- Public-facing configuration should be based on
  `pipeline_config.example.toml`, not local machine-specific overrides


## Study design at a glance

The active project is built around an explicitly unpaired study design:

- healthy and tumour cohorts are analysed as separate groups rather than
  matched tumour-normal pairs
- tumour-versus-healthy comparisons are used for risk or enrichment questions
- tumour-only analyses are used for clinicopathological associations
- RNA integration is a downstream objective built on the DNA and haplotype
  backbone



## Canonical computational core

The default launcher now runs the smallest coherent thesis core:

1. `01_check_and_index.sh`
2. `02_dna_qc.sh`
3. `02b_more-qc.sh`
4. `05_variant_calling.sh`
5. `05b_force_genotype_union_sites.py`
6. `06_annotation.sh`
7. `07_merge_annotations.py`
8. `11_SNPs.py`
9. `12_stats_enrichment.py`
10. `15_haplotypes.sh`
11. `15_haplotype_stats.R`
12. `16_excel_harmonisation.py`
13. `17_snp_association.py`
14. `18_haplotype_association.py`
15. `19b_rna_qc.py`
16. `20_objective2_isoform_expression_association.py`
17. `20b_panel_gene_expression.py`
18. `23_rna_integration.py`



## Recommended thesis story

For the thesis text, the work is best described in this order:

1. Technical QC and sequencing validity
2. Variant calling and annotation
3. Variant landscape and regional mapping
4. `GSDMB`-focused descriptive analyses
5. Common SNP summary
6. SNP tumour-versus-healthy enrichment
7. Haplotype construction and haplotype results
8. Tumour clinical association analyses
9. RNA or expression integration
10. Supplementary sensitivity, interpretation, and dashboard layers


## Repository grouping

The repository is now best understood in five groups.

### Core preprocessing and discovery

- `01_check_and_index.sh`
- `02_dna_qc.sh`
- `02b_more-qc.sh`
- `05_variant_calling.sh`
- `05b_force_genotype_union_sites.py`
- `06_annotation.sh`
- `07_merge_annotations.py`

### Core descriptive and results-support layers

- `08_mapping.py`
- `09_gsdmb_only.py`
- `10_variant_stats.py`
- `11_SNPs.py`
- `12_stats_enrichment.py`
- `15_haplotypes.sh`
- `15_haplotype_stats.R`

### Core association layers

- `16_excel_harmonisation.py`
- `17_snp_association.py`
- `18_haplotype_association.py`
- `19b_rna_qc.py`
- `20_objective2_isoform_expression_association.py`
- `20b_panel_gene_expression.py`
- `23_rna_integration.py`

### Supplementary or exploratory analyses

- `02c_zerocovinfo.py`
- `03_qc_visualisation.py`
- `04_technical_audit.py`
- `04b_amplicon_failure_visualisation.py`
- `04b_amplicon_failure_panel.py`
- `07b_variant_qc.py`
- `09b_landscape_comparison.py`
- `13_permutation_testing.py`
- `14_interactive_dashboard.py`
- `17b_rare_variant_association.py`
- `19_1000g_haplotype_comparison.py`
- `21_core_haplotype_comparison.py`
- `22_haplotype_first_interpretation.py`
- `24_external_cross_validation.py`
- `25_rs11078928_rs869402_haplotype_focus.py`
- `26_snp_functional_interpretation.py`

### Utility, audit, and assembly helpers

These scripts are useful, but they are not part of the canonical rerun launched by `run_full_pipeline.sh`.

- `08b_bed_vs_detected_excel.py`
- `27_sample_similarity_audit.py`
- `28_collect_significant_survival_curves.py`

## How to run

### Canonical thesis rerun

```bash
cd /path/to/tfm
./run_full_pipeline.sh
```

Equivalent explicit form:

```bash
./run_full_pipeline.sh --profile core
```

### Supplementary reruns

```bash
./run_full_pipeline.sh --profile exploratory
./run_full_pipeline.sh --profile utility
./run_full_pipeline.sh --profile full
```

### Bounded reruns

```bash
./run_full_pipeline.sh --profile core --from 11 --to 18
./run_full_pipeline.sh --profile exploratory --from 13 --to 26
```


## Key documentation

- [docs/notes/PIPELINE.md](docs/notes/PIPELINE.md): execution order, story order, and repository grouping
- [docs/notes/SCRIPT_MAP.md](docs/notes/SCRIPT_MAP.md): one-line purpose and status for each numbered stage
- [docs/notes/RESULTS_GUIDE.md](docs/notes/RESULTS_GUIDE.md): which outputs are main-text, supplementary, or secondary
- [docs/notes/METHODS_ALIGNMENT.md](docs/notes/METHODS_ALIGNMENT.md): thesis-facing wording guardrails
- [docs/notes/REPRODUCIBILITY.md](docs/notes/REPRODUCIBILITY.md): rerun strategy and environment notes
- [docs/notes/GITHUB_REPOSITORY_GUIDE.md](docs/notes/GITHUB_REPOSITORY_GUIDE.md): what belongs in a public-facing repository
- [script_classification_table.tsv](script_classification_table.tsv): concise machine-readable classification of scripts

