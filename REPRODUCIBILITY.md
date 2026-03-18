# Reproducibility Guide

This document explains how to run the pipeline reproducibly and how the repository maps onto the three environments used in the original analysis workflow.

## 1. Environment Strategy

The original project used three distinct environments:

- `bam-steps` via **micromamba** for shell-based bioinformatics processing
- `tfm_env` via **Python virtualenv** for the Python analysis scripts
- `vep_env` via **conda** for Ensembl VEP annotation

For GitHub sharing, this repository now includes two reusable dependency files:

- [environment.yml](environment.yml) for a shared conda/micromamba environment
- [requirements.txt](requirements.txt) for the historical `tfm_env` virtual-environment workflow

## 2. Exact Environment Activation

### `bam-steps`

Create it once:

```bash
micromamba create -n bam-steps samtools bcftools mosdepth -c bioconda -c conda-forge
```

Activate it:

```bash
micromamba activate bam-steps
```

Deactivate it:

```bash
micromamba deactivate
```

Typical scripts run here:

- `01_check_and_index.sh`
- `02_dna_qc.sh`
- `02b_more-qc.sh`
- `05_variant_calling.sh`
- `15_haplotypes.sh`

### `tfm_env`

Create it once in the original style:

```bash
python3 -m venv tfm_env
source tfm_env/bin/activate
pip install -r requirements.txt
```

Activate it later:

```bash
source tfm_env/bin/activate
```

Deactivate it:

```bash
deactivate
```

Typical scripts run here:

- `02c_zerocovinfo.py`
- `03_qc_visualisation.py`
- `04_technical_audit.py`
- `07_merge_annotations.py`
- `07b_variant_qc.py`
- `08_mapping.py`
- `09_gsdmb_only.py`
- `10_variant_stats.py`
- `11_SNPs.py`
- `12_stats_enrichment.py`
- `13_permutation_testing.py`
- `14_interactive_dashboard.py`
- `16_excel_harmonisation.py`
- `17_snp_association.py`
- `18_haplotype_association.py`
- `19_1000g_haplotype_comparison.py`

### `vep_env`

Activate it with:

```bash
conda activate vep_env
```

Deactivate it with:

```bash
conda deactivate
```

Typical script run here:

- `06_annotation.sh`

This environment must provide Ensembl VEP and its local annotation resources.

## 3. Shared Dependency Files

### `requirements.txt`

Use this if you want to recreate the Python-analysis environment in the same style as the original `tfm_env` workflow.

```bash
python3 -m venv tfm_env
source tfm_env/bin/activate
pip install -r requirements.txt
```


The `19_1000g_haplotype_comparison.py` step additionally requires:

- an indexed phased 1000 Genomes VCF or BCF
- a matching 1000 Genomes sample panel containing sample, population, and super-population labels
- the `pysam` dependency, now included in both `requirements.txt` and `environment.yml`

### `environment.yml`

Use this if you prefer a single reproducible environment through conda or micromamba.

With conda:

```bash
conda env create -f environment.yml
conda activate tfm-gsdmb
```

With micromamba:

```bash
micromamba create -f environment.yml
micromamba activate tfm-gsdmb
```

This shared environment is convenient for documentation and portability, but it does not replace the fact that the original workflow used the three-environment split described above.

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

## 6. Shared Logic

To keep scripts consistent, common operations were centralised in:

- [pipeline_utils.py](pipeline_utils.py)
- [pipeline_validation.py](pipeline_validation.py)
- [association_runtime.py](association_runtime.py)

This reduces the risk that one script applies different labelling, frequency, or control-group logic than another.

## 7. Validation

Several scripts now include lightweight validation steps such as:

- input file existence checks
- required-column checks
- non-empty filtered dataset checks
- expected tissue-label checks
- percentage bounds checks

These checks are intentionally simple but useful for catching broken inputs early.

## 8. Minimal Verification

A practical lightweight verification routine is:

```bash
python -m py_compile 07_merge_annotations.py 08_mapping.py 09_gsdmb_only.py 10_variant_stats.py 11_SNPs.py 12_stats_enrichment.py 13_permutation_testing.py 14_interactive_dashboard.py association_runtime.py pipeline_utils.py pipeline_validation.py 17_snp_association.py 18_haplotype_association.py 19_1000g_haplotype_comparison.py
bash -n 01_check_and_index.sh 05_variant_calling.sh 06_annotation.sh
Rscript -e "parse(file='15_haplotype_stats.R')"
```

That does not replace full reruns on data, but it helps catch syntax issues before sharing.
