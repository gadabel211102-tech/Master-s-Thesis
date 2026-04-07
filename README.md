# GSDMB Thesis Workflow

This repository contains the reproducible workflow for the GSDMB-focused Master's thesis analyses. The project now distinguishes explicitly between:

1. `CORE PIPELINE`
2. `EXPLORATORY / VALIDATION ANALYSES`
3. `UTILITY / REPORTING / DIAGNOSTIC SCRIPTS`

The default launcher now runs the smallest coherent reproducible core needed for the two main thesis objectives only:

1. Genetic, epidemiological, and SNP-haplotype comparison in breast and endometrial cohorts
2. RNA, isoform, and clinicopathological integration

## What Changed

The repository was audited against the current file state, not an older script list. The default launcher no longer treats every available downstream script as part of the canonical workflow.

- `./run_full_pipeline.sh` now runs the `core` profile by default.
- Exploratory and validation branches remain available through explicit profiles.
- Utility, reporting, and diagnostic scripts remain available but sit outside the default scientific rerun.
- Stage `23` no longer requires stage `22` to exist; if the stage-22 interpretation workbook is absent, stage `23` falls back to the significant haplotypes already present in the stage-15 association workbook.
- The former breast paired-validation branch was removed from the runnable pipeline after confirming that the breast cohort was ultimately processed as unpaired tumour and normal BAM sets rather than tumour-matched pairs.
- Stage `24` now performs external comparison against female-only 1000 Genomes only; the earlier TCGA branch was removed from the active workflow because the required TCGA germline genotype files were controlled-access.

## Workflow Split

### Core pipeline

These steps are the default reproducible thesis workflow:

1. `01_check_and_index.sh`
2. `02_dna_qc.sh`
3. `02b_more-qc.sh`
4. `05_variant_calling.sh`
5. `06_annotation.sh`
6. `07_merge_annotations.py`
7. `11_SNPs.py`
8. `12_stats_enrichment.py`
9. `15_haplotypes.sh`
10. `15_haplotype_stats.R`
11. `16_excel_harmonisation.py`
12. `17_snp_association.py`
13. `18_haplotype_association.py`
14. `19b_rna_qc.py`
15. `20_objective2_isoform_expression_association.py`
16. `20b_panel_gene_expression.py`
17. `23_rna_integration.py`

This sequence produces the core analytical objects used in the thesis:

- DNA pass manifests and filtered cohort VCFs
- Canonical annotated variant workbook
- Common-SNP backbone and primary enrichment results
- Phased haplotype matrix and haplotype statistics workbook
- Harmonised clinical master workbook
- SNP-clinical and haplotype-clinical association workbooks
- RNA QC manifest
- Isoform-expression and BAM-derived gene-expression workbooks
- Final RNA/SNP/haplotype integration workbook

### Exploratory and validation analyses

These scripts remain available but are not part of the default scientific rerun:

- `09b_landscape_comparison.py`
- `13_permutation_testing.py`
- `19_1000g_haplotype_comparison.py`
- `22_haplotype_first_interpretation.py`
- `24_external_cross_validation.py`
- `25_rs11078928_rs869402_haplotype_focus.py`
- `26_snp_functional_interpretation.py`

These are best described as sensitivity checks, external comparison, biological interpretation, or cohort-specific follow-up layers.

### Utility, reporting, and diagnostics

These scripts help with QC review, descriptive visualisation, collaborator outputs, dashboards, or preparation work, but are not required to reproduce the principal thesis results:

- `02c_zerocovinfo.py`
- `03_qc_visualisation.py`
- `04_technical_audit.py`
- `04b_amplicon_failure_visualisation.py`
- `04b_amplicon_failure_panel.py`
- `07b_variant_qc.py`
- `08_mapping.py`
- `09_gsdmb_only.py`
- `10_variant_stats.py`
- `14_interactive_dashboard.py`
- `21_core_haplotype_comparison.py`
- shared helper modules such as `pipeline_utils.py`, `pipeline_validation.py`, `association_runtime.py`, and the RNA utility modules

## How To Run

### Run the core pipeline only

```bash
cd /home/gadeaalonsoj/tfm
./run_full_pipeline.sh
```

That is equivalent to:

```bash
./run_full_pipeline.sh --profile core
```

### Run only exploratory / validation analyses

```bash
./run_full_pipeline.sh --profile exploratory
```

### Run only utility / reporting scripts

```bash
./run_full_pipeline.sh --profile utility
```

### Run everything in step order

```bash
./run_full_pipeline.sh --profile full
```

### Run a subset of a profile

```bash
./run_full_pipeline.sh --profile utility --from 04 --to 14
./run_full_pipeline.sh --profile exploratory --from 24 --to 26
```

## Script Classification Table

The full current audit is recorded in [script_classification_table.tsv](script_classification_table.tsv).

Columns:

- `script_name`
- `current_path`
- `category`
- `role_summary`
- `required_by_core`
- `recommended_default_run`

## Repository Files That Anchor The Workflow

- [`run_full_pipeline.sh`](run_full_pipeline.sh): profile-aware launcher
- [`pipeline_config.toml`](pipeline_config.toml): shared paths and thresholds
- [`pipeline_config.example.toml`](pipeline_config.example.toml): public-safe template for GitHub or new machines
- [`script_classification_table.tsv`](script_classification_table.tsv): current script audit and classification
- [`docs/notes/PIPELINE.md`](docs/notes/PIPELINE.md): workflow overview and execution logic
- [`docs/notes/REPRODUCIBILITY.md`](docs/notes/REPRODUCIBILITY.md): environment and rerun guidance
- [`docs/notes/METHODS_ALIGNMENT.md`](docs/notes/METHODS_ALIGNMENT.md): thesis-facing methods guidance
- [`docs/notes/GITHUB_REPOSITORY_GUIDE.md`](docs/notes/GITHUB_REPOSITORY_GUIDE.md): what to upload, what to keep local, and a suggested public repo layout

## Outputs

Most generated outputs are written under:

```text
/home/gadeaalonsoj/tfm/analysis_results
```

Representative core output locations:

- `analysis_results/07_annotated_variants`
- `analysis_results/11_common_snps`
- `analysis_results/12_snp_enrichment`
- `analysis_results/15_haplotype_phasing`
- `analysis_results/15_haplotype_statistics`
- `analysis_results/17_snp_clinical_associations`
- `analysis_results/18_haplotype_clinical_associations`
- `analysis_results/19b_objective2_rna_qc`
- `analysis_results/20_isoform_expression_associations`
- `analysis_results/20b_panel_gene_expression`
- `analysis_results/23_rna_integration`

## Environment Strategy

The workflow still uses the original three-runtime split:

- `bam-steps` for shell-based BAM/QC/calling/phasing stages
- `tfm_env` or another explicit Python interpreter for the Python analyses
- `vep_env` or the repo-local VEP runtime for annotation

See [docs/notes/REPRODUCIBILITY.md](docs/notes/REPRODUCIBILITY.md) for the exact rerun guidance.

## GitHub Preparation

For a public or supervisor-facing GitHub repository, do not upload raw cohort BAMs, generated `analysis_results`, private clinical workbooks, or local environments. The recommended upload plan is documented in [docs/notes/GITHUB_REPOSITORY_GUIDE.md](docs/notes/GITHUB_REPOSITORY_GUIDE.md).

The tracked `pipeline_config.toml` in this working directory still reflects the local machine used during analysis. For GitHub, upload [`pipeline_config.example.toml`](pipeline_config.example.toml) as the reusable template and keep any machine-specific override in an untracked local file such as `pipeline_config.local.toml` or via the `PIPELINE_CONFIG` environment variable.

## Thesis Methods Mapping

For the thesis Methods section, describe the core workflow as the canonical analysis path. Describe external comparison, interpretation, dashboards, and technical audit layers as secondary or supplementary analyses rather than part of the mandatory default pipeline.

The detailed wording guidance lives in [docs/notes/METHODS_ALIGNMENT.md](docs/notes/METHODS_ALIGNMENT.md).
