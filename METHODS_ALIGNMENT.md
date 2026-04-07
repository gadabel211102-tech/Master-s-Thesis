# Methods Alignment Note

## What belongs in the main Methods workflow

The main thesis Methods section should describe the canonical workflow as:

1. DNA BAM validation and QC (`01`, `02`, `02b`)
2. Variant calling and annotation (`05`, `06`, `07`)
3. Common-SNP definition and principal tumour-versus-control SNP analysis (`11`, `12`)
4. Haplotype phasing and principal haplotype statistics (`15`, `15R`)
5. Clinical harmonisation (`16`)
6. SNP-clinical and haplotype-clinical association analyses (`17`, `18`)
7. RNA QC gating (`19b`)
8. Workbook-first isoform analysis plus BAM-derived panel-gene quantification (`20`, `20b`)
9. Final RNA/SNP/haplotype integration (`23`)

That is the cleanest defensible description of the main reproducible thesis workflow.

## What should be described as secondary or supplementary analyses

These scripts should be described as follow-up, sensitivity, validation, or interpretation layers rather than as part of the canonical default pipeline:

- `09b_landscape_comparison.py`
- `13_permutation_testing.py`
- `19_1000g_haplotype_comparison.py`
- `22_haplotype_first_interpretation.py`
- `24_external_cross_validation.py`
- `25_rs11078928_rs869402_haplotype_focus.py`
- `26_snp_functional_interpretation.py`

Suggested wording pattern:

- "Secondary validation analyses included external haplotype comparison and targeted haplotype follow-up."
- "Biological interpretation layers were generated after the main association analyses and were not part of the default discovery workflow."

Do not describe a breast paired tumour-normal validation branch in the active Methods workflow, because the breast BAM material was ultimately processed as unpaired tumour and normal cohorts rather than confirmed matched pairs.

Describe stage `24` as external comparison against female-only 1000 Genomes rather than as a TCGA-backed branch, because the required TCGA germline genotype files were controlled-access and not available for this project.

## What should stay out of the main Methods workflow narrative

These scripts are useful for reproducibility, review, and presentation, but they should not inflate the scientific Methods description unless a specific figure or diagnostic section requires them:

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

Suggested wording pattern:

- "Additional QC, descriptive visualisation, dashboard, and export utilities were used for review and reporting but were not required for the principal analytical rerun."

## Why this split is preferable

This separation improves:

- reproducibility, because the default launcher matches the true canonical analysis
- supervisor review, because primary and secondary layers are no longer mixed together
- thesis writing, because the Methods section can stay focused on the main scientific workflow
- future reruns, because the expensive optional branches are explicit rather than implicit
- interpretation, because the reader can see what is discovery, what is validation, and what is reporting support
