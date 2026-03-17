# Reproducibility Guide

This document explains how to run the pipeline reproducibly and what should be adjusted before sharing or reusing the repository.

## 1. Environment

The Python environment is defined in [environment.yml](environment.yml).

Create it with:

```bash
conda env create -f environment.yml
conda activate tfm-gsdmb
```

The pipeline also depends on external bioinformatics tools that must be installed separately, including `samtools`, `bcftools`, `mosdepth`, `vep`, Ion Torrent calling tools, `beagle`, and `Rscript`.

## 2. Configuration

Shared runtime settings are defined in [pipeline_config.toml](pipeline_config.toml).

The most important fields are:

- `paths`: where inputs, manifests, results, and harmonised files live
- `thresholds`: cutoffs used across SNP and haplotype analyses
- `grouping`: tumour/control definitions and pooled-normal behaviour

Before running on another machine, update all absolute paths.

## 3. Shared Logic

To keep scripts consistent, common operations were centralised in:

- [pipeline_utils.py](pipeline_utils.py)
- [pipeline_validation.py](pipeline_validation.py)
- [association_runtime.py](association_runtime.py)

This reduces the risk that one script applies different labelling, frequency, or control-group logic than another.

## 4. Control Definition

The current comparison design uses pooled controls for stronger statistical power.

That means:

- `Breast tumour` is compared against `Breast normal + Endometrium normal`
- `Endometrium tumour` is compared against `Breast normal + Endometrium normal`
- `Global` compares all tumours against all controls

User-facing outputs now label this arm as `Control`.

## 5. Sample-Size Bias Control

Plots used for cohort comparison were updated so they do not overstate differences caused purely by unequal sample counts.

Examples include:

- carrier percentages instead of raw carrier counts in mapping plots
- density or percentage-based summaries in comparison plots
- cohort-normalised frequencies where direct comparison is intended

Descriptive count plots were left as counts only when the goal is to describe composition rather than compare prevalence.

## 6. Validation

Several scripts now include lightweight validation steps such as:

- input file existence checks
- required-column checks
- non-empty filtered dataset checks
- expected tissue-label checks
- percentage bounds checks

These checks are intentionally simple but useful for catching broken inputs early.

## 7. Recommended Practice Before Publishing

Before pushing this repository to GitHub, review the following:

- remove or ignore patient-level data files
- review whether result directories should be versioned
- replace local absolute paths if you want others to run the code
- add a project-specific license
- add maintainer contact and formal citation details to the README

## 8. Minimal Verification

A practical lightweight verification routine is:

```bash
python -m py_compile 07_merge_annotations.py 08_mapping.py 09_gsdmb_only.py 10_variant_stats.py 11_SNPs.py 12_stats_enrichment.py 13_permutation_testing.py 14_interactive_dashboard.py association_runtime.py pipeline_utils.py pipeline_validation.py
bash -n 01_check_and_index.sh 05_variant_calling.sh 06_annotation.sh
Rscript -e "parse(file='15_haplotype_stats.R')"
```

That does not replace full reruns on data, but it helps catch syntax issues before sharing.
