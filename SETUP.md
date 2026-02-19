# Environment Setup Guide

Step-by-step instructions for reproducing the three environments used in this pipeline.

---

## Prerequisites

- Linux or WSL2 (Ubuntu 22.04+ recommended)
- [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) installed
- [Conda](https://docs.conda.io/en/latest/miniconda.html) installed (for VEP environment)
- Python 3.11+
- R 4.3+ (for haplotype analysis)

---

## 1. Bioinformatics Environment (Micromamba)

Used for: BAM QC, variant calling, haplotype phasing.

```bash
# Create environment from file
micromamba env create -f environment.yml

# Or create manually
micromamba create -n bam-steps samtools bcftools mosdepth -c bioconda -c conda-forge

# Activate
micromamba activate bam-steps
```

**Verify installation:**
```bash
samtools --version
bcftools --version
mosdepth --version
```

---

## 2. Visualisation Environment (Python venv)

Used for: all Python analysis scripts (02c, 03–16).

```bash
# Create virtual environment
python3 -m venv tfm_env

# Activate
source tfm_env/bin/activate   # Linux/macOS
# tfm_env\Scripts\activate    # Windows

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt
```

**Verify installation:**
```bash
python -c "import pandas, seaborn, matplotlib, sklearn, plotly, scipy; print('All OK')"
```

---

## 3. VEP Annotation Environment (Conda)

Used for: `06_annotation.sh` only.

This environment requires a working Ensembl VEP installation. Refer to the [official VEP documentation](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html) for full installation instructions.

```bash
# Activate the existing VEP environment
conda activate vep_env

# Verify
vep --help
```

### Required VEP Resources

Download and place in your VEP cache directory (default: `~/.vep` or `$ROOT_DIR/.vep`):

| File | Source |
|---|---|
| VEP cache (`homo_sapiens_vep_112_GRCh38`) | [Ensembl FTP](https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/) |
| Reference FASTA (`Homo_sapiens.GRCh38.dna.toplevel.fa.gz`) | [Ensembl FTP](https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna_index/) |
| CADD scores (`CADD_GRCh38_whole_genome_SNVs.tsv.gz`) | [CADD website](https://cadd.gs.washington.edu/download) |
| dbNSFP (`dbNSFP5.3.1a_grch38.gz`) | [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP) |

---

## 4. R Environment

Used for: `17_haplo_stats.r`.

```r
# Install CRAN packages
install.packages(c("ggplot2", "dplyr", "tidyr", "openxlsx"))

# Install haplo.stats
install.packages("haplo.stats")
```

**Verify:**
```r
library(haplo.stats)
library(ggplot2)
cat("R environment OK\n")
```

---

## WSL2-Specific Notes

- The VEP annotation script (`06_annotation.sh`) automatically unsets conflicting Perl environment variables (`PERL5LIB`, `PERL_LOCAL_LIB_ROOT`) that commonly cause issues in WSL2. No manual action is required.
- If Micromamba is not found after installation, ensure your shell config (`.bashrc` / `.zshrc`) has been sourced: `source ~/.bashrc`.
