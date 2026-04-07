# Workflow Overview

## Core principle

The repository is organised around a small canonical thesis workflow plus optional branches.

- `CORE` means required to reproduce the principal thesis-ready analytical objects.
- `EXPLORATORY_VALIDATION` means optional sensitivity, external comparison, interpretation, or follow-up analysis.
- `UTILITY_REPORTING` means diagnostics, descriptive plotting, dashboards, helper scripts, or collaborator exports.

## Default launcher behaviour

[`run_full_pipeline.sh`](../../run_full_pipeline.sh) now defaults to `--profile core`.

Supported profiles:

- `core`
- `exploratory`
- `utility`
- `full`

The launcher also supports `--from` and `--to` so optional branches can be rerun in a bounded range.

## Core execution order

| Step | Script | Why it stays in core | Main output(s) |
| --- | --- | --- | --- |
| 01 | `01_check_and_index.sh` | Validates DNA BAM inputs before analysis starts | Indexed BAMs |
| 02 | `02_dna_qc.sh` | Defines pass/fail QC state for DNA samples | QC summaries |
| 02b | `02b_more-qc.sh` | Builds the PASS manifests used by downstream calling | PASS manifests |
| 05 | `05_variant_calling.sh` | Produces the filtered cohort VCFs | Filtered cohort VCFs |
| 06 | `06_annotation.sh` | Adds the annotation layer used by the SNP branch | Annotated VCFs |
| 07 | `07_merge_annotations.py` | Builds the canonical variant workbook | `GSDMB_Annotated_Variants.xlsx` |
| 11 | `11_SNPs.py` | Defines the common-SNP backbone | Common-SNP workbook |
| 12 | `12_stats_enrichment.py` | Runs the principal tumour-versus-control SNP testing | SNP enrichment workbook |
| 15 | `15_haplotypes.sh` | Produces the phased genotype matrix | `phased_genotypes.tsv` |
| 15R | `15_haplotype_stats.R` | Produces the principal haplotype results workbook | `GSDMB_Haplotype_Results.xlsx` |
| 16 | `16_excel_harmonisation.py` | Builds the harmonised clinical master | `MASTER_SNP_plus_clinical_HARMONISED.xlsx` |
| 17 | `17_snp_association.py` | Principal SNP-clinical modelling | SNP association workbook |
| 18 | `18_haplotype_association.py` | Principal haplotype-clinical modelling | Haplotype association workbook |
| 19b | `19b_rna_qc.py` | Defines which RNA rows can enter objective 2 | `RNA_QC_Manifest.tsv` |
| 20 | `20_objective2_isoform_expression_association.py` | Workbook-first isoform and RNA association layer | Stage-20 workbook |
| 20b | `20b_panel_gene_expression.py` | BAM-derived panel-gene quantification for objective 2 | Stage-20b workbook |
| 23 | `23_rna_integration.py` | Final integration of RNA, SNP, and haplotype signals | Stage-23 integration workbook |

## Important dependency rule

Stage `23` is considered part of the core workflow. It now behaves as follows:

- if stage `22` exists, stage `23` uses the curated stage-22 significant haplotypes
- if stage `22` is absent, stage `23` falls back to the significant haplotypes already available in the stage-15 association workbook

That keeps the interpretation layer optional while preserving a clean core RNA integration rerun.

## Exploratory / validation branch

These steps remain available but do not run by default:

| Step | Script | Rationale |
| --- | --- | --- |
| 09b | `09b_landscape_comparison.py` | Follow-up review of descriptive landscape differences |
| 13 | `13_permutation_testing.py` | Sensitivity analysis for the SNP enrichment branch |
| 19 | `19_1000g_haplotype_comparison.py` | External haplotype context against 1000 Genomes |
| 22 | `22_haplotype_first_interpretation.py` | Interpretation layer for significant haplotypes |
| 24 | `24_external_cross_validation.py` | External validation against female-only 1000 Genomes |
| 25 | `25_rs11078928_rs869402_haplotype_focus.py` | Targeted follow-up on a specific haplotype question |
| 26 | `26_snp_functional_interpretation.py` | SNP biological-context synthesis |

The previously drafted stage `27` paired breast validation branch is not part of the active workflow because the breast BAM material was not generated as a tumour-matched paired set.

## Utility / reporting branch

| Step | Script | Rationale |
| --- | --- | --- |
| 02c | `02c_zerocovinfo.py` | Formats zero-coverage summaries into a report |
| 03 | `03_qc_visualisation.py` | QC plots and sample tables |
| 04 | `04_technical_audit.py` | Amplicon technical diagnostics |
| 04bV / 04bP | `04b_*` scripts | Detailed failure visualisation and overview panel |
| 07b | `07b_variant_qc.py` | Variant-level QC assurance |
| 08 / 09 / 10 | descriptive mapping/statistics scripts | Presentation-facing landscape and composition outputs |
| 14 | `14_interactive_dashboard.py` | Interactive dashboard |
| 21 | `21_core_haplotype_comparison.py` | Collaborator-facing summary export |

## Recommended manual run logic

### Canonical thesis rerun

```bash
./run_full_pipeline.sh
```

### Optional exploratory rerun

```bash
./run_full_pipeline.sh --profile exploratory
```

### Optional utility rerun

```bash
./run_full_pipeline.sh --profile utility
```

### Everything

```bash
./run_full_pipeline.sh --profile full
```

## Output map

Core scientific outputs live mainly in:

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

Optional outputs remain in their existing folders so older paths are preserved.
