#!/usr/bin/env python3
"""
Script 20b: BAM-derived Panel Gene Expression Quantification and Association Analysis
======================================================================================

PURPOSE
-------
Quantifies expression of all 63 genes in the Ion RNA AmpliSeq panel directly
from the RNA BAM files using samtools bedcov (read-depth per amplicon), then:

  1. Builds a normalised expression matrix (CPM + GAPDH/RPS18 ratio normalisation)
  2. Associates SNP genotype / haplotype with total GSDMB BAM expression and
     all other panel genes
  3. Correlates qPCR isoform fractions (G1-G4 from script 20 Excel) with
     BAM-derived expression of immune/pyroptosis panel genes â€” the key
     biological question: does G4 (pyroptotic isoform) fraction correlate
     with immune cell infiltration or downstream effectors?
  4. Compares BAM-derived total GSDMB against the GSDMB column in the
     qPCR Excel as a cross-validation step
  5. Produces publication-quality figures for all analyses

DESIGN RATIONALE
----------------
GSDMB isoforms (G1/G2/G3/G3b/G4) CANNOT be quantified from these BAMs because
the panel has only one GSDMB amplicon (AMPL6874839) that spans shared exonic
sequence across all isoforms. Isoform quantification requires the dedicated
qPCR assay (script 20 Excel). This script therefore focuses exclusively on:
  - Total GSDMB expression (single amplicon, whole-gene signal)
  - All 62 other panel genes

PANEL GENE GROUPS (from IAD258499_4_DataSheet.csv)
---------------------------------------------------
  Pyroptosis:     GSDMB, GSDMA, GSDMC, GSDMD, DFNA5(GSDME), CASP1, CASP3,
                  CASP7, IL1B, IL18
  Chr17q12 locus: ORMDL3, ZPBP2, IKZF3, ERBB2, GRB7, PGAP3, MIEN1, STARD3,
                  PNMT, PPP1R1B, TCAP, LRRC3C, PSMD3, GSDMA
  Immune markers: CD3E, CD4, CD8A, CD14, CD33, CD68, CD163, CD274(PD-L1),
                  CD279(PD-1), CD80, CD86, FOXP3, CTLA4, MS4A1, PTPRC,
                  NCAM1, ITGA4, ITGAM, CEACAM8, FUT4, MRC1
  Cytokines:      IFNG, TNF, IL6, IL13, IL15, IL25, IL2RA, TGFB1, TNFSF14,
                  TSLP, ELANE, MPO, GZMA, GZMB
  Housekeeping:   GAPDH, RPS18, PCNA, MKI67

INPUTS
------
  --rna-bed          IAD258499_4_DataSheet.csv (panel design, from Ion AmpliSeq)
  --bam-dirs         Comma-separated root dirs to search for RNA BAMs
                     e.g. breast/tumour/rna,breast/normal/rna,...
                     Defaults to the standard project layout under --project-root
  --project-root     Base project directory (default: current repository root)
  --rna-qc-manifest  RNA_QC_Manifest.tsv from script 19b (filters to QC-passed
                     BAMs; if absent, all BAMs found are processed)
  --expr-xlsx        qPCR isoform Excel from script 20 (Sheet3 + eva sheets)
                     Used for cross-validation and isoform-vs-panelgene correlations
  --master           Harmonised clinical master from script 16b
  --variant-workbook GSDMB_Annotated_Report_Fixed.xlsx (for SNP genotypes)
  --haplotype-input  phased_genotypes.tsv from script 15
  --out-dir          Output directory

OUTPUTS
-------
  20b_Panel_Gene_Expression.xlsx  â€” main results workbook with sheets:
    raw_counts          â€” raw bedcov read counts per gene per sample
    normalised_cpm      â€” counts per million (library-size normalised)
    normalised_ratio    â€” GAPDH+RPS18 ratio-normalised expression
    sample_manifest     â€” per-sample QC and metadata
    gsdmb_crossval      â€” BAM total GSDMB vs qPCR GSDMB column comparison
    snp_gene_assoc      â€” SNP dosage vs each panel gene (linear regression)
    haplotype_gene_assocâ€” haplotype dosage vs each panel gene
    isoform_immune_corr â€” qPCR isoform fractions vs immune/pyroptosis genes
    clinical_gene_assoc â€” panel gene expression vs clinical variables
    gene_group_summary  â€” median expression per gene group per cohort

  20b_GSDMB_CrossValidation.png   â€” scatter: BAM GSDMB vs qPCR GSDMB
  20b_Expression_Heatmap.png      â€” normalised expression heatmap (all genes)
  20b_Gene_Group_Distributions.pngâ€” boxplots per gene group per cohort
  20b_Isoform_Immune_Heatmap.png  â€” isoform fraction vs immune gene correlations
  20b_SNP_PanelGene_Heatmap.png   â€” SNP vs panel gene association heatmap
  20b_Haplotype_PanelGene_Heatmap.png
  20b_Chr17_Locus_Expression.png  â€” chr17q12 locus gene expression per cohort

USAGE
-----
  python3 20b_panel_gene_expression.py

  # With explicit paths:
  python3 20b_panel_gene_expression.py \\
      --project-root /path/to/tfm \\
      --rna-bed /path/to/IAD258499_4_DataSheet.csv \\
      --expr-xlsx /path/to/expression_workbook.xlsx \\
      --master /path/to/MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx \\
      --variant-workbook /path/to/GSDMB_Annotated_Report_Fixed.xlsx \\
      --haplotype-input /path/to/phased_genotypes.tsv \\
      --out-dir /path/to/output/

DEPENDENCIES
------------
  samtools (must be on PATH)
  python packages: pandas, numpy, scipy, matplotlib, seaborn (optional)
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from association_runtime import script20_defaults
from pipeline_utils import build_genomic_variant_id_series, combine_gnomad_nfe, common_nfe_variant_mask, get_paths
from sample_identity_utils import attach_analysis_sample_ids, build_analysis_sample_map

try:
    import seaborn as sns
    _HAS_SNS = True
except ImportError:
    _HAS_SNS = False

warnings.filterwarnings("ignore")

# â”€â”€ DEFAULT PATHS â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

PATHS = get_paths()
SCRIPT20_DEFAULTS = script20_defaults()
ROOT = PATHS["project_root"]
DEFAULTS = {
    "project_root": ROOT,
    "rna_bed": SCRIPT20_DEFAULTS["rna_bed"],
    "expr_xlsx": SCRIPT20_DEFAULTS["expr_xlsx"],
    "master": SCRIPT20_DEFAULTS["master"],
    "variant_workbook": SCRIPT20_DEFAULTS["variant_workbook"],
    "haplotype_input": SCRIPT20_DEFAULTS["haplotype_input"],
    "rna_qc_manifest": SCRIPT20_DEFAULTS["rna_qc_manifest"],
    "out_dir": PATHS["results_dir"] / "20b_panel_gene_expression",
}

# RNA BAM search directories relative to project root
RNA_BAM_SUBDIRS = [
    ("Breast_Tumour",       "breast/tumour/rna"),
    ("Breast_Normal",       "breast/normal/rna"),
    ("Endometrial_Tumour",  "endometrium/tumour/rna"),
    ("Endometrial_Normal",  "endometrium/normal/rna"),
]

# Tissue label â†’ (cohort, tissue, analysis_group, role)
TMAP = {
    "breast tumor":         ("Breast",        "Tumour",  "Breast_Tumour",        "primary_tumour"),
    "breast tumour":        ("Breast",        "Tumour",  "Breast_Tumour",        "primary_tumour"),
    "breast normal":        ("Breast",        "Healthy", "Breast_Normal",        "context_control"),
    "endometrial cancer":   ("Endometrial",   "Tumour",  "Endometrial_Tumour",   "primary_tumour"),
    "endometrial normal":   ("Endometrial",   "Healthy", "Endometrial_Normal",   "context_control"),
    "endometrial_normal":   ("Endometrial",   "Healthy", "Endometrial_Normal",   "context_control"),
}

# Gene groups for biological interpretation
GENE_GROUPS = {
    "Pyroptosis": [
        "GSDMB", "GSDMA", "GSDMC", "GSDMD", "DFNA5",
        "CASP1", "CASP3", "CASP7", "IL1B", "IL18",
    ],
    "Chr17q12_Locus": [
        "ORMDL3", "ZPBP2", "IKZF3", "ERBB2", "GRB7",
        "PGAP3", "MIEN1", "STARD3", "PNMT", "PPP1R1B",
        "TCAP", "LRRC3C", "PSMD3",
    ],
    "Immune_Markers": [
        "CD3E", "CD4", "CD8A", "CD14", "CD33", "CD68",
        "CD163", "CD274", "PDCD1", "CD80", "CD86", "FOXP3",
        "CTLA4", "MS4A1", "PTPRC", "NCAM1", "ITGA4",
        "ITGAM", "CEACAM8", "FUT4", "MRC1",
    ],
    "Cytokines_Inflammation": [
        "IFNG", "TNF", "IL6", "IL13", "IL15", "IL25",
        "IL2RA", "TGFB1", "TNFSF14", "TSLP", "ELANE",
        "MPO", "GZMA", "GZMB",
    ],
    "Housekeeping": ["GAPDH", "RPS18", "PCNA", "MKI67"],
}

# Genes used for normalisation
HOUSEKEEPING_GENES = ["GAPDH", "RPS18"]

# Minimum read count to consider a sample usable for a given gene
MIN_READS_PER_GENE = 10

# Minimum samples with expression data to run an association test
MIN_N_FOR_ASSOC = 10

# FDR threshold
FDR_THRESHOLD = 0.10

RNA_GATE_COLS = [
    "RNA_QC_Analysis_Ready",
    "RNA_QC_Exploratory_Ready",
    "RNA_QC_Final_Status",
    "RNA_QC_Exclusion_Reason",
    "RNA_QC_Eligibility_Note",
]
META_COLS = [
    "snp_code",
    "bam_path",
    "bam_group",
    "cohort",
    "tissue",
    "analysis_role",
    "RNA_QC_Analysis_Ready",
    "RNA_QC_Exploratory_Ready",
    "RNA_QC_Final_Status",
    "RNA_QC_Exclusion_Reason",
    "RNA_QC_Eligibility_Note",
    "RNA_BAM_Usage_Class",
    "RNA_QC_Gate_Used",
]

# qPCR isoform columns from script 20 Excel
ISOFORM_COLS = ["G1", "G2", "G3", "G3b", "G4",
                "pyroptotic_isoform_fraction", "non_pyroptotic_isoform_fraction",
                "G4_vs_total", "G2_vs_total",
                "exon6_containing_fraction", "exon6_lacking_fraction",
                "exon6_balance_log2", "exon7_containing_fraction",
                "exon7_lacking_fraction", "exon7_balance_log2"]

# Clinical variables (reuse from script 17/20 definitions)
CLINICAL_VARS_BREAST = {
    "BREAST_GRADE_NUMERIC":     ("continuous", "Tumour grade"),
    "BREAST_KI67_NUMERIC":      ("continuous", "KI67"),
    "BREAST_ER_BIN":            ("binary",     "ER positive"),
    "BREAST_PR_BIN":            ("binary",     "PR positive"),
    "BREAST_RECURRENCE_DERIVED":("binary",     "Recurrence/progression"),
    "BREAST_EXITUS_DERIVED":    ("binary",     "Death"),
    "canon__age":               ("continuous", "Age"),
}
CLINICAL_VARS_ENDO = {
    "ENDO_FIGO_NUMERIC":        ("continuous", "FIGO stage"),
    "ENDO_GRADE_NUMERIC":       ("continuous", "Tumour grade"),
    "ENDO_RISK_ORDINAL":        ("continuous", "Risk of recurrence"),
    "ENDO_LVSI_BIN":            ("binary",     "LVSI"),
    "ENDO_MYOINV_BIN":          ("binary",     "Myometrial invasion â‰¥50%"),
    "ENDO_MSI_BIN":             ("binary",     "MSI-H"),
    "ENDO_PD_BIN":              ("binary",     "Disease progression"),
    "ENDO_EXITUS_BIN":          ("binary",     "Death"),
    "ENDO_ER_BIN":              ("binary",     "ER positive"),
    "ENDO_PR_BIN":              ("binary",     "PR positive"),
    "ENDO_NEEC_BIN":            ("binary",     "Non-endometrioid histology"),
    "canon__age":               ("continuous", "Age"),
}


# â”€â”€ UTILITIES â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def _s(v) -> str:
    return "" if pd.isna(v) else str(v).strip()


def _nk(v) -> str:
    return re.sub(r"[^A-Z0-9]+", "", _s(v).upper())


def _sx(v) -> Optional[str]:
    """Extract canonical SNP code from a BAM filename stem."""
    t = _s(v)
    pattens = [
        (r"^(SNP_(?:AT|EN|MN|MT-T)_\d+)",         lambda m: m.group(1).upper()),
        (r"^SNP_DNA_(EN|MN|AT)_(\d+)",             lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^SNP_DNA_MT[-_]T_(\d+)",                lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^DNA_SNP_(EN|MN|AT)_(\d+)",             lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^DNA_(AT|EN|MN)_(\d+)",                 lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^DNA_(?:SNP_)?MT[-_]T_(\d+)",           lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^MAMAH2_MT[-_]T_(\d+)",                lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^MAMAH2_MT[-_]N_(\d+)",                lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^RNA_SNP_(AT|EN|MN|MT-T|MT-N)[-_]?(\d+)", lambda m: f"SNP_{'MT-T' if m.group(1).upper() == 'MT-N' else m.group(1).upper()}_{m.group(2)}"),
        (r"^SNP_RNA_(AT|EN|MN|MT-T|MT-N)_(\d+)", lambda m: f"SNP_{'MT-T' if m.group(1).upper() == 'MT-N' else m.group(1).upper()}_{m.group(2)}"),
        (r"^SNP_(AT|EN|MN|MT-T|MT-N)_RNA_(\d+)", lambda m: f"SNP_{'MT-T' if m.group(1).upper() == 'MT-N' else m.group(1).upper()}_{m.group(2)}"),
        (r"^RNA_(AT|EN|MN)_(\d+)",                 lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^RNA_MT[-_]T_(\d+)",                    lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^RNA_MT[-_]N_(\d+)",                    lambda m: f"SNP_MT-T_{m.group(1)}"),
    ]
    for pat, fmt in pattens:
        m = re.match(pat, t, re.I)
        if m:
            return fmt(m)
    return None


def _bh_fdr(pvals: pd.Series) -> pd.Series:
    """Benjamini-Hochberg FDR correction."""
    a = np.asarray(pvals, dtype=float)
    out = np.full(a.shape, np.nan)
    m = np.isfinite(a)
    if not m.any():
        return pd.Series(out, index=pvals.index)
    p = a[m]
    o = np.argsort(p)
    r = p[o]
    n = len(r)
    adj = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, r[i] * n / (i + 1))
        adj[i] = min(prev, 1.0)
    b = np.empty(n)
    b[o] = adj
    out[m] = b
    return pd.Series(out, index=pvals.index)


def _add_fdr(df: pd.DataFrame, p_col: str,
             group_cols: List[str]) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    df["FDR"] = np.nan
    for _, ix in df.groupby(group_cols, dropna=False).groups.items():
        df.loc[ix, "FDR"] = _bh_fdr(
            pd.to_numeric(df.loc[ix, p_col], errors="coerce")
        ).values
    df["Nominal_Sig"] = pd.to_numeric(df[p_col], errors="coerce") < 0.05
    df["FDR_Sig"] = pd.to_numeric(df["FDR"], errors="coerce") < FDR_THRESHOLD
    return df


def _deduplicate_qc_manifest(qc_manifest: pd.DataFrame) -> pd.DataFrame:
    if qc_manifest is None or qc_manifest.empty or "snp_code" not in qc_manifest.columns:
        return pd.DataFrame()
    keep_cols = [
        c for c in [
            "snp_code",
            "analysis_group",
            "analysis_role",
            "cohort",
            "tumour_normal",
            "RNA_BAM_Group",
            "RNA_BAM_Name",
            "RNA_BAM_Path",
        ] + RNA_GATE_COLS
        if c in qc_manifest.columns
    ]
    dedup = qc_manifest[keep_cols].copy()
    dedup["snp_code"] = dedup["snp_code"].astype(str).str.strip()
    dedup = dedup[
        dedup["snp_code"].notna()
        & dedup["snp_code"].ne("")
        & dedup["snp_code"].ne("nan")
    ].copy()
    sort_cols = [c for c in ["RNA_QC_Analysis_Ready", "RNA_QC_Exploratory_Ready"] if c in dedup.columns]
    if sort_cols:
        dedup = dedup.sort_values(sort_cols, ascending=False)
    dedup_keys = [c for c in ["snp_code", "RNA_BAM_Path", "RNA_BAM_Name", "analysis_group", "analysis_role"] if c in dedup.columns]
    if not dedup_keys:
        dedup_keys = ["snp_code"]
    return dedup.drop_duplicates(dedup_keys).reset_index(drop=True)


# â”€â”€ PANEL DESIGN â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def load_panel_design(panel_path: Path) -> pd.DataFrame:
    """
    Parse either the AmpliSeq CSV datasheet (genomic coordinates) or the
    transcript-coordinate Designed BED into a clean amplicon table.

    Retuns columns: gene, amplicon_id, chrom, start, end, length
    """
    rows = []
    if panel_path.suffix.lower() == ".bed":
        with open(panel_path) as fh:
            for line in fh:
                line = line.strip()
                if (not line
                        or line.startswith("#")
                        or line.startswith("track")):
                    continue
                parts = line.split("	")
                if len(parts) < 6:
                    continue
                chrom, start, end, amplicon, _, gene = parts[:6]
                if (not amplicon or amplicon == "."
                        or not chrom or chrom == "."):
                    continue
                try:
                    s, e = int(start), int(end)
                except ValueError:
                    continue
                rows.append({
                    "gene": gene.strip(),
                    "amplicon_id": amplicon.strip(),
                    "chrom": chrom.strip(),
                    "start": s,
                    "end": e,
                    "length": e - s,
                    "coordinate_space": "transcript",
                })
    else:
        with open(panel_path) as fh:
            header_seen = False
            for line in fh:
                line = line.strip()
                if line.startswith("#"):
                    continue
                if not header_seen:
                    header_seen = True
                    continue  # skip column header line
                parts = line.split(",")
                if len(parts) < 9:
                    continue
                gene = parts[1].strip()
                amplicon = parts[4].strip()
                chrom = parts[5].strip()
                start = parts[7].strip()
                end = parts[8].strip()
                if (not amplicon or amplicon == "."
                        or not chrom or chrom == "."):
                    continue
                try:
                    s, e = int(start), int(end)
                except ValueError:
                    continue
                rows.append({
                    "gene": gene,
                    "amplicon_id": amplicon,
                    "chrom": chrom,
                    "start": s,
                    "end": e,
                    "length": e - s,
                    "coordinate_space": "genomic",
                })

    df = pd.DataFrame(rows).drop_duplicates("amplicon_id").reset_index(drop=True)
    if df.empty:
        raise ValueError(f"No usable panel intervals found in {panel_path}")
    print(f"  Panel design: {len(df)} amplicons covering "
          f"{df['gene'].nunique()} genes")
    return df


# â”€â”€ BAM DISCOVERY â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def discover_bams(project_root: Path,
                  qc_manifest: Optional[pd.DataFrame],
                  strict_only: bool = False) -> pd.DataFrame:
    """
    Find RNA BAMs under the standard project directory layout.

    When a QC manifest is available, the default behaviour is to quantify both
    strict analysis-ready and exploratory-ready BAMs. Strict-only mode remains
    available but is opt-in so the workbook can preserve the distinction.
    """
    qc_samples = _deduplicate_qc_manifest(qc_manifest)
    qc_ready = pd.DataFrame()
    passed_codes: Optional[set] = None
    gate_used = "exploratory_inclusive"
    if strict_only:
        gate_used = "strict_only"
    if not qc_samples.empty:
        ready_col = "RNA_QC_Analysis_Ready" if strict_only else (
            "RNA_QC_Exploratory_Ready" if "RNA_QC_Exploratory_Ready" in qc_samples.columns else "RNA_QC_Analysis_Ready"
        )
        if ready_col in qc_samples.columns:
            qc_ready = qc_samples.loc[qc_samples[ready_col].fillna(False).astype(bool)].copy()
            passed_codes = set(qc_ready["snp_code"].dropna().astype(str))
            gate_label = "strict analysis-ready" if strict_only else "strict + exploratory"
            print(f"  QC manifest: restricting to {len(passed_codes)} {gate_label} RNA samples")
        else:
            print("  QC manifest loaded but no ready-flag column found; using all BAMs")
    rows = []
    for group_label, subdir in RNA_BAM_SUBDIRS:
        bam_dir = project_root / subdir
        if not bam_dir.exists():
            print(f"  [WARN] BAM directory not found: {bam_dir}")
            continue
        parts = group_label.split("_")
        cohort = parts[0]
        tissue = parts[1] if len(parts) > 1 else "Unknown"
        role = "primary_tumour" if tissue == "Tumour" else "context_control"
        for bam in sorted(bam_dir.glob("*.bam")):
            code = _sx(bam.stem)
            if passed_codes is not None and code not in passed_codes:
                continue
            if not qc_ready.empty:
                matched = qc_ready[qc_ready["snp_code"].astype(str).eq(str(code))].copy()
                if matched.empty:
                    continue
                if "RNA_BAM_Path" in matched.columns and matched["RNA_BAM_Path"].notna().any():
                    path_match = matched["RNA_BAM_Path"].astype(str).eq(str(bam))
                    if path_match.any():
                        matched = matched.loc[path_match]
                    else:
                        matched = pd.DataFrame()
                if matched.empty and "RNA_BAM_Name" in qc_ready.columns:
                    matched = qc_ready[
                        qc_ready["snp_code"].astype(str).eq(str(code))
                        & qc_ready["RNA_BAM_Name"].astype(str).eq(bam.name)
                    ].copy()
                if matched.empty and "analysis_group" in qc_ready.columns:
                    matched = qc_ready[
                        qc_ready["snp_code"].astype(str).eq(str(code))
                        & qc_ready["analysis_group"].astype(str).eq(group_label)
                    ].copy()
                if matched.empty:
                    continue
            rows.append({
                "snp_code": code,
                "bam_path": str(bam),
                "bam_name": bam.name,
                "bam_group": group_label,
                "cohort": cohort,
                "tissue": tissue,
                "analysis_role": role,
                "RNA_QC_Gate_Used": gate_used,
            })
    df = pd.DataFrame(rows)
    if df.empty:
        print("  [WARN] No RNA BAMs found under project root.")
        return df
    if not qc_samples.empty:
        qc_merge = qc_samples.copy()
        merge_keys = ["snp_code"]
        if "RNA_BAM_Name" in qc_merge.columns:
            qc_merge = qc_merge.rename(columns={"RNA_BAM_Name": "bam_name"})
            merge_keys = ["snp_code", "bam_name"]
        merge_cols = merge_keys + [c for c in ["analysis_group", "analysis_role", "cohort", "tumour_normal"] + RNA_GATE_COLS if c in qc_merge.columns and c not in merge_keys]
        df = df.merge(qc_merge[merge_cols].drop_duplicates(merge_keys), on=merge_keys, how="left", suffixes=("", "_manifest"))
        for col in ["cohort", "analysis_role"]:
            manifest_col = f"{col}_manifest"
            if manifest_col in df.columns:
                df[col] = df[col].fillna(df[manifest_col])
                df = df.drop(columns=[manifest_col])
        df["RNA_QC_Analysis_Ready"] = df.get("RNA_QC_Analysis_Ready", False).fillna(False).astype(bool)
        df["RNA_QC_Exploratory_Ready"] = df.get("RNA_QC_Exploratory_Ready", False).fillna(df["RNA_QC_Analysis_Ready"]).astype(bool)
        df["RNA_QC_Final_Status"] = df.get("RNA_QC_Final_Status", pd.Series(index=df.index, dtype=object)).fillna("Included_without_manifest")
        df["RNA_QC_Exclusion_Reason"] = df.get("RNA_QC_Exclusion_Reason", pd.Series(index=df.index, dtype=object)).fillna("")
        df["RNA_QC_Eligibility_Note"] = df.get("RNA_QC_Eligibility_Note", pd.Series(index=df.index, dtype=object)).fillna("")
        df["RNA_BAM_Usage_Class"] = np.where(
            df["RNA_QC_Analysis_Ready"],
            "Strict",
            np.where(df["RNA_QC_Exploratory_Ready"], "Exploratory", "Excluded")
        )
    else:
        for col in RNA_GATE_COLS:
            if col not in df.columns:
                df[col] = np.nan
        df["RNA_BAM_Usage_Class"] = "Unfiltered"
    print(f"  Found {len(df)} RNA BAMs across {df['bam_group'].nunique()} groups")
    if "RNA_BAM_Usage_Class" in df.columns:
        counts = df["RNA_BAM_Usage_Class"].value_counts(dropna=False).to_dict()
        print(f"  BAM usage classes: {counts}")
    return df


# â”€â”€ QUANTIFICATION â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def _make_bed_file(panel: pd.DataFrame, bed_out: Path) -> None:
    """Write a temporary BED file for samtools bedcov."""
    with open(bed_out, "w") as fh:
        for _, row in panel.iterrows():
            fh.write(
                f"{row['chrom']}\t{row['start']}\t{row['end']}\t"
                f"{row['gene']}\n"
            )


def _run_bedcov(bam_path: str, bed_path: Path) -> Dict[str, float]:
    """
    Run samtools bedcov to get total read depth per region.
    Retuns dict mapping gene -> mean depth (total_bases / amplicon_length).
    samtools bedcov outputs: chr start end name total_bases
    """
    result = subprocess.run(
        ["samtools", "bedcov", str(bed_path), bam_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools bedcov failed for {bam_path}:\n{result.stderr[:300]}"
        )
    counts: Dict[str, float] = {}
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        gene       = parts[3]
        total_bases = float(parts[4])
        length     = int(parts[2]) - int(parts[1])
        # Mean depth = total bases covered / amplicon length
        counts[gene] = total_bases / length if length > 0 else 0.0
    return counts


def quantify_all_bams(
    bam_df: pd.DataFrame,
    panel: pd.DataFrame,
    out_dir: Path,
) -> pd.DataFrame:
    """
    Run samtools bedcov for every BAM and return a raw counts matrix.
    Columns: snp_code, bam_path, bam_group, cohort, tissue, analysis_role,
             + one column per gene (mean depth)
    Uses a cached TSV if present and --recompute not requested.
    """
    cache_path = out_dir / "20b_raw_counts_cache.tsv"
    if cache_path.exists():
        cached = pd.read_csv(cache_path, sep="\t")
        cached_paths = set(cached.get("bam_path", pd.Series(dtype=object)).dropna().astype(str))
        requested_paths = set(bam_df.get("bam_path", pd.Series(dtype=object)).dropna().astype(str))
        if cached_paths == requested_paths and len(cached_paths) == len(bam_df):
            print(f"  Loading cached raw counts from {cache_path}")
            return cached
        print("  Raw-count cache sample set is stale; recomputing to include the current BAM selection")

    if not shutil.which("samtools"):
        raise EnvironmentError(
            "samtools not found on PATH. Please install samtools or ensure "
            "it is available before running this script."
        )

    # Write BED once
    bed_tmp = out_dir / "_tmp_panel.bed"
    _make_bed_file(panel, bed_tmp)

    genes = panel["gene"].tolist()
    rows = []
    n = len(bam_df)
    for i, (_, bam_row) in enumerate(bam_df.iterrows(), 1):
        bam_path = bam_row["bam_path"]
        print(f"  [{i:3d}/{n}] quantifying {Path(bam_path).name} ...",
              end="\r", flush=True)
        try:
            counts = _run_bedcov(bam_path, bed_tmp)
        except Exception as e:
            print(f"\n  [WARN] {Path(bam_path).name}: {e}")
            counts = {}
        row = {col: bam_row[col] for col in META_COLS if col in bam_row.index}
        row.update({
            "snp_code": bam_row["snp_code"],
            "bam_path": bam_path,
            "bam_group": bam_row["bam_group"],
            "cohort": bam_row["cohort"],
            "tissue": bam_row["tissue"],
            "analysis_role": bam_row["analysis_role"],
        })
        for gene in genes:
            row[gene] = counts.get(gene, np.nan)
        rows.append(row)

    print()  # newline after progress
    bed_tmp.unlink(missing_ok=True)

    df = pd.DataFrame(rows)
    df.to_csv(cache_path, sep="\t", index=False)
    print(f"  Raw counts saved to {cache_path}")
    return df


# â”€â”€ NORMALISATION â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def normalise(raw: pd.DataFrame, panel: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Two normalisation strategies:

    1. CPM (counts per million) â€” library-size normalised mean depth.
       Library size proxy = sum of all gene mean depths per sample.

    2. Housekeeping-ratio normalisation â€” each gene divided by the geometric
       mean of GAPDH and RPS18 depths. This controls for RNA integrity and
       input variation better than CPM for targeted panels.

    Retuns (cpm_df, ratio_df) â€” same shape as raw, only gene columns updated.
    """
    meta_cols = [c for c in META_COLS if c in raw.columns]
    genes = [c for c in raw.columns if c not in meta_cols]

    # â”€â”€ CPM â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    gene_mat = raw[genes].apply(pd.to_numeric, errors="coerce")
    lib_size = gene_mat.sum(axis=1).replace(0, np.nan)
    cpm_mat = gene_mat.div(lib_size, axis=0) * 1e6
    cpm = raw[meta_cols].copy()
    for g in genes:
        cpm[g] = cpm_mat[g]

    # â”€â”€ Housekeeping ratio â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    hk_available = [g for g in HOUSEKEEPING_GENES if g in gene_mat.columns]
    if hk_available:
        hk_vals = gene_mat[hk_available].apply(pd.to_numeric, errors="coerce")
        # Geometric mean across housekeeping genes per sample
        log_hk = np.log1p(hk_vals)
        geo_mean_hk = np.expm1(log_hk.mean(axis=1)).replace(0, np.nan)
    else:
        print("  [WARN] No housekeeping genes found; ratio normalisation "
              "falls back to CPM.")
        geo_mean_hk = lib_size / 1e6

    ratio_mat = gene_mat.div(geo_mean_hk, axis=0)
    ratio = raw[meta_cols].copy()
    for g in genes:
        ratio[g] = ratio_mat[g]

    print(f"  Normalised {len(raw)} samples Ã— {len(genes)} genes")
    return cpm, ratio


def build_bam_usage_summary(raw: pd.DataFrame) -> pd.DataFrame:
    if raw.empty:
        return pd.DataFrame({"Note": ["No RNA BAM rows were quantified."]})
    summary = raw.copy()
    if "RNA_BAM_Usage_Class" not in summary.columns:
        summary["RNA_BAM_Usage_Class"] = "Unfiltered"
    if "RNA_QC_Final_Status" not in summary.columns:
        summary["RNA_QC_Final_Status"] = "Unavailable"
    out = (
        summary.groupby(["bam_group", "analysis_role", "RNA_BAM_Usage_Class", "RNA_QC_Final_Status"], dropna=False)
        .agg(
            quantified_bams=("snp_code", "count"),
            unique_samples=("snp_code", pd.Series.nunique),
        )
        .reset_index()
        .sort_values(["analysis_role", "bam_group", "RNA_BAM_Usage_Class"])
    )
    return out


def build_normalisation_audit(raw: pd.DataFrame, strict_only: bool) -> pd.DataFrame:
    gate_label = "strict analysis-ready only" if strict_only else "strict plus exploratory RNA BAMs"
    rows = [
        {
            "Layer": "BAM selection",
            "Current_Method": gate_label,
            "Notes": "Gate provenance is carried per BAM row via RNA_QC_Analysis_Ready, RNA_QC_Exploratory_Ready, RNA_QC_Final_Status, and RNA_BAM_Usage_Class.",
        },
        {
            "Layer": "raw_counts sheet",
            "Current_Method": "samtools bedcov mean amplicon depth per gene",
            "Notes": "The sheet label is retained for compatibility, but values are depth-derived panel quantification rather than literal molecule counts.",
        },
        {
            "Layer": "normalised_cpm sheet",
            "Current_Method": "CPM-like scaling across panel mean-depth totals",
            "Notes": "This is an intenal targeted-panel normalisation, not TPM/FPKM.",
        },
        {
            "Layer": "normalised_ratio sheet",
            "Current_Method": "Housekeeping ratio normalisation using GAPDH and RPS18",
            "Notes": "Stage 23 uses this layer for BAM-derived RNA integration.",
        },
        {
            "Layer": "Excel RNA units",
            "Current_Method": "Source workbook values treated as already-normalised assay output and preserved as provided",
            "Notes": "No extra repo-level rescaling is applied to the Excel gene columns. BAM-vs-Excel checks therefore remain cross-platform concordance checks unless the compared workbook columns are known to be in directly compatible units.",
        },
    ]
    if not raw.empty and "RNA_BAM_Usage_Class" in raw.columns:
        counts = raw["RNA_BAM_Usage_Class"].value_counts(dropna=False).to_dict()
        rows.append({
            "Layer": "BAM usage counts",
            "Current_Method": str(counts),
            "Notes": "These counts reflect the BAM rows actually quantified into this workbook.",
        })
    return pd.DataFrame(rows)


# â”€â”€ LOAD ISOFORM DATA (from script 20 Excel) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def load_isoform_data(expr_path: Path) -> pd.DataFrame:
    """
    Load isoform expression from the script 20 qPCR Excel workbook.
    Extracts G1-G4, GSDMB (total qPCR), and derived isoform fractions.
    Retuns one row per sample with snp_code + isoform columns.
    """
    try:
        xls = pd.ExcelFile(expr_path)
    except Exception as e:
        print(f"  [WARN] Could not open expression workbook: {e}")
        return pd.DataFrame()

    # Try the cleaned expression sheet first (written by script 20 if run)
    if "expression_cleaned" in xls.sheet_names:
        df = pd.read_excel(expr_path, sheet_name="expression_cleaned")
        iso_cols = [c for c in ISOFORM_COLS if c in df.columns]
        keep = ["snp_code"] + iso_cols
        keep = [c for c in keep if c in df.columns]
        if "snp_code" in df.columns and iso_cols:
            print(f"  Loaded isoform data from 'expression_cleaned' sheet: "
                  f"{len(df)} rows, {len(iso_cols)} isoform columns")
            return df[keep].dropna(subset=["snp_code"]).copy()

    # Fall back to Sheet3 (raw workbook)
    if "Sheet3" not in xls.sheet_names:
        print("  [WARN] Expression workbook has neither 'expression_cleaned' "
              "nor 'Sheet3'; skipping isoform data.")
        return pd.DataFrame()

    df = pd.read_excel(expr_path, sheet_name="Sheet3")
    for c in ["G1", "G2", "G3", "G3b", "G4", "GSDMB"]:
        df[c] = pd.to_numeric(df.get(c), errors="coerce")

    # Compute derived columns. The legacy pyroptotic/non-pyroptotic columns are
    # retained for backward compatibility; exon-defined variables are the
    # biologically preferred interpretation layer.
    df["pyroptotic_isoform_fraction"] = df[["G3", "G3b", "G4"]].sum(
        axis=1, min_count=1)
    df["non_pyroptotic_isoform_fraction"] = df[["G1", "G2"]].sum(
        axis=1, min_count=1)
    df["isoform_total"] = df[["G1", "G2", "G3", "G3b", "G4"]].sum(
        axis=1, min_count=1)
    df["G4_vs_total"] = np.where(
        df["isoform_total"] > 0, df["G4"] / df["isoform_total"], np.nan)
    df["G2_vs_total"] = np.where(
        df["isoform_total"] > 0, df["G2"] / df["isoform_total"], np.nan)
    canonical_total = df[["G1", "G2", "G3", "G3b", "G4"]].sum(axis=1, min_count=1)
    exon6_containing = df[["G3", "G3b", "G4"]].sum(axis=1, min_count=1)
    exon6_lacking = df[["G1", "G2"]].sum(axis=1, min_count=1)
    exon7_containing = df[["G1", "G3", "G3b"]].sum(axis=1, min_count=1)
    exon7_lacking = df[["G2", "G4"]].sum(axis=1, min_count=1)
    df["exon6_containing_fraction"] = np.where(
        canonical_total > 0, exon6_containing / canonical_total, np.nan)
    df["exon6_lacking_fraction"] = np.where(
        canonical_total > 0, exon6_lacking / canonical_total, np.nan)
    df["exon6_balance_log2"] = np.log2(
        (pd.to_numeric(exon6_containing, errors="coerce") + 1e-6) /
        (pd.to_numeric(exon6_lacking, errors="coerce") + 1e-6))
    df["exon7_containing_fraction"] = np.where(
        canonical_total > 0, exon7_containing / canonical_total, np.nan)
    df["exon7_lacking_fraction"] = np.where(
        canonical_total > 0, exon7_lacking / canonical_total, np.nan)
    df["exon7_balance_log2"] = np.log2(
        (pd.to_numeric(exon7_containing, errors="coerce") + 1e-6) /
        (pd.to_numeric(exon7_lacking, errors="coerce") + 1e-6))

    # Derive snp_code from sample name columns
    for raw_col in ["CODIGO JC", "NOMBRE DE LA MUESTRA"]:
        if raw_col in df.columns:
            df["snp_code"] = df[raw_col].map(_sx)
            if df["snp_code"].notna().sum() > 0:
                break

    iso_cols = [c for c in ISOFORM_COLS if c in df.columns]
    keep = ["snp_code"] + iso_cols
    keep = [c for c in keep if c in df.columns]
    result = df[keep].dropna(subset=["snp_code"]).copy()
    print(f"  Loaded isoform data from 'Sheet3': {len(result)} rows, "
          f"{len(iso_cols)} isoform columns")
    return result


# â”€â”€ LOAD GENETIC DATA (SNP genotypes + haplotypes) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def load_genetic_data(
    variant_wb_path: Path,
    haplotype_path: Path,
    min_nfe: float = 0.01,
    min_hap_freq: float = 0.02,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load SNP genotype long table and haplotype dosage table.
    Mirrors the logic from script 20 load_backbone / load_phased.

    Retuns:
        snp_long   â€” snp_code Ã— Variant_ID with Dosage
        hap_long   â€” snp_code Ã— Haplotype_ID with Dosage
        hap_freqs  â€” haplotype frequency table
    """
    # â”€â”€ SNP backbone â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    try:
        vdf = pd.read_excel(variant_wb_path, sheet_name="Biological_Annotations")
    except Exception as e:
        print(f"  [WARN] Could not load variant workbook: {e}")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    # Identify gnomAD NFE AF columns
    if "gnomADe_NFE_AF" in vdf.columns or "gnomADg_NFE_AF" in vdf.columns:
        vdf = combine_gnomad_nfe(vdf)
    nfe_col = "gnomAD_NFE_AF_combined" if "gnomAD_NFE_AF_combined" in vdf.columns else next(
        (c for c in vdf.columns if re.search(r"gnomad.*nfe.*af|nfe.*af", str(c), re.I)),
        None,
    )
    if nfe_col is None:
        print("  [WARN] No gnomAD NFE AF column found; using all variants")
        vdf["_nfe"] = 1.0
        nfe_col = "_nfe"
        common_mask = pd.Series(True, index=vdf.index)
    else:
        vdf[nfe_col] = pd.to_numeric(vdf[nfe_col], errors="coerce")
        common_mask = common_nfe_variant_mask(vdf, min_nfe) if nfe_col == "gnomAD_NFE_AF_combined" else vdf[nfe_col].gt(min_nfe).fillna(False)
    vdf = vdf[common_mask].copy()

    # Derive formal snp_code plus analysis_sample_id. Raw runs are collapsed
    # only when their phased genotype columns are exactly identical.
    vdf = attach_analysis_sample_ids(
        vdf,
        raw_col="Sample",
        phased_path=haplotype_path,
        extract_snp_code=_sx,
    )
    rsid_col = next(
        (c for c in vdf.columns if "existing_variation" in c.lower()), None)
    if rsid_col:
        vdf["rsID"] = vdf[rsid_col].astype(str).str.extract(r"(rs\d+)")
    else:
        vdf["rsID"] = np.nan
    vdf["POS"] = pd.to_numeric(vdf.get("POS", vdf.get("Pos")), errors="coerce")
    vdf["REF"] = vdf.get("REF", pd.Series("N", index=vdf.index)).astype(str).str.upper()
    vdf["ALT"] = vdf.get("ALT", pd.Series("N", index=vdf.index)).astype(str).str.upper()
    vdf["SNP_Label"] = (vdf["POS"].astype("Int64").astype(str)
                        + "_" + vdf["REF"] + ">" + vdf["ALT"])
    vdf["Variant_ID"] = build_genomic_variant_id_series(
        vdf["Variant_Key"] if "Variant_Key" in vdf.columns else None,
        vdf["CHROM"] if "CHROM" in vdf.columns else None,
        vdf["POS"],
        vdf["REF"],
        vdf["ALT"],
    )

    # Normalise GT
    def _dose(gt):
        t = _s(gt).replace("|", "/")
        if t in {"", ".", "./.", ".|."}:
            return np.nan
        try:
            return float(sum(int(x) for x in t.split("/")))
        except Exception:
            return np.nan

    vdf["Dosage"] = vdf.get("GT", pd.Series(np.nan, index=vdf.index)).map(_dose)
    snp_long = (vdf[["analysis_sample_id", "snp_code", "Variant_ID", "SNP_Label",
                      "Dosage", nfe_col]]
                .dropna(subset=["analysis_sample_id", "snp_code", "Variant_ID"])
                .drop_duplicates(["analysis_sample_id", "Variant_ID"])
                .rename(columns={nfe_col: "gnomAD_NFE_AF"})
                .copy())
    print(f"  SNP backbone: {snp_long['Variant_ID'].nunique()} variants, "
          f"{snp_long['analysis_sample_id'].nunique()} analysis samples")

    # â”€â”€ Haplotype data â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    try:
        phased = pd.read_csv(haplotype_path, sep="\t", dtype=str)
    except Exception as e:
        print(f"  [WARN] Could not load haplotype file: {e}")
        return snp_long, pd.DataFrame(), pd.DataFrame()

    phased["POS"] = pd.to_numeric(phased["POS"], errors="coerce")
    snp_labels_in_backbone = set(snp_long["SNP_Label"].dropna())
    phased["REF"] = phased["REF"].astype(str).str.upper()
    phased["ALT"] = phased["ALT"].astype(str).str.upper()
    phased["SNP_Label"] = (phased["POS"].astype("Int64").astype(str)
                           + "_" + phased["REF"] + ">" + phased["ALT"])
    phased = phased[phased["SNP_Label"].isin(snp_labels_in_backbone)].copy()

    sample_cols = [c for c in phased.columns
                   if c not in {"CHROM", "POS", "ID", "REF", "ALT", "SNP_Label"}]
    identity_map = build_analysis_sample_map(sample_cols, haplotype_path, _sx)
    representative_cols = (
        identity_map.sort_values(["analysis_sample_id", "raw_sample_name"])
        .drop_duplicates("analysis_sample_id")[["raw_sample_name", "analysis_sample_id", "snp_code"]]
    )

    # Build haplotype strings from one representative raw column per retained
    # analysis sample ID.
    hap_rows = []
    for _, id_row in representative_cols.iterrows():
        col = id_row["raw_sample_name"]
        code = id_row["snp_code"]
        analysis_sample_id = id_row["analysis_sample_id"]
        h1_parts, h2_parts = [], []
        for _, row in phased.iterrows():
            gt = _s(row[col])
            if gt in {"", ".", ".|.", "./.", "nan"}:
                h1_parts.append("N")
                h2_parts.append("N")
                continue
            sep = "|" if "|" in gt else "/" if "/" in gt else None
            if sep:
                a, b = gt.split(sep, 1)
                h1_parts.append(a)
                h2_parts.append(b)
            else:
                h1_parts.append("N")
                h2_parts.append("N")
        hap_rows.append({
            "analysis_sample_id": analysis_sample_id,
            "snp_code": code,
            "Sample_phased": col,
            "hap1": "".join(h1_parts),
            "hap2": "".join(h2_parts),
        })

    hap_df = pd.DataFrame(hap_rows)

    # Count haplotype frequencies
    all_haps: Dict[str, int] = {}
    for col in ["hap1", "hap2"]:
        for hp in hap_df[col].dropna():
            if hp and "N" not in hp:
                all_haps[hp] = all_haps.get(hp, 0) + 1
    total_chrom = sum(all_haps.values())
    hap_freqs = pd.DataFrame([
        {"Haplotype": k, "Count": v,
         "Global_Freq": v / total_chrom if total_chrom else np.nan}
        for k, v in sorted(all_haps.items(), key=lambda kv: (-kv[1], kv[0]))
    ])
    hap_freqs = hap_freqs[
        hap_freqs["Global_Freq"] >= min_hap_freq
    ].reset_index(drop=True)
    hap_freqs["Haplotype_ID"] = [f"H{i+1}" for i in range(len(hap_freqs))]

    # Build haplotype dosage long table
    hap_long_rows = []
    for _, hf_row in hap_freqs.iterrows():
        for _, s_row in hap_df.iterrows():
            ok = int("N" not in str(s_row["hap1"])
                     and "N" not in str(s_row["hap2"]))
            dose = (np.nan if not ok else
                    int(str(s_row["hap1"]) == hf_row["Haplotype"])
                    + int(str(s_row["hap2"]) == hf_row["Haplotype"]))
            hap_long_rows.append({
                "analysis_sample_id": s_row["analysis_sample_id"],
                "snp_code":     s_row["snp_code"],
                "Haplotype_ID": hf_row["Haplotype_ID"],
                "Haplotype":    hf_row["Haplotype"],
                "Global_Freq":  hf_row["Global_Freq"],
                "Dosage":       dose,
                "Carrier":      np.nan if not ok else int(dose >= 1),
            })
    hap_long = pd.DataFrame(hap_long_rows)
    print(f"  Haplotypes: {len(hap_freqs)} retained "
          f"(freq â‰¥ {min_hap_freq:.0%}), "
          f"{hap_long['analysis_sample_id'].nunique()} analysis samples")
    return snp_long, hap_long, hap_freqs


# â”€â”€ CROSS-VALIDATION: BAM GSDMB vs qPCR GSDMB â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def cross_validate_gsdmb(
    ratio: pd.DataFrame,
    isoform: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merge BAM-derived total GSDMB (ratio-normalised mean depth) with the
    GSDMB total expression column from the qPCR Excel, then compute
    Spearman correlation as a cross-validation metric.
    """
    if "GSDMB" not in ratio.columns or isoform.empty:
        return pd.DataFrame({"Note": ["Cross-validation skipped: missing data"]})
    if "GSDMB" not in isoform.columns:
        return pd.DataFrame({"Note": ["No GSDMB column in isoform workbook"]})

    bam_gsdmb = ratio[["snp_code", "cohort", "tissue", "GSDMB"]].copy()
    bam_gsdmb = bam_gsdmb.rename(columns={"GSDMB": "GSDMB_BAM"})
    qpcr_gsdmb = isoform[["snp_code", "GSDMB"]].rename(
        columns={"GSDMB": "GSDMB_qPCR"})

    merged = bam_gsdmb.merge(qpcr_gsdmb, on="snp_code", how="inner")
    merged = merged.dropna(subset=["GSDMB_BAM", "GSDMB_qPCR"])

    rows = []
    for label, sub in [("All", merged)] + [
        (grp, merged[merged["cohort"] == grp])
        for grp in merged["cohort"].unique()
    ]:
        if len(sub) < 5:
            continue
        rho, p = stats.spearmanr(sub["GSDMB_BAM"], sub["GSDMB_qPCR"],
                                 nan_policy="omit")
        rows.append({
            "Group":    label,
            "N":        len(sub),
            "Spearman_rho": round(rho, 4),
            "P_Value":  round(p, 4),
            "Note":     "Cross-platform concordance: BAM-derived ratio-normalised depth vs Excel total GSDMB source value",
        })
    result = pd.DataFrame(rows)
    print(f"  Cross-validation: {len(result)} group comparisons")
    return merged, result


# â”€â”€ ASSOCIATION: SNP / HAPLOTYPE vs PANEL GENES â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def assoc_genetic_vs_genes(
    expr_df: pd.DataFrame,
    genetic_long: pd.DataFrame,
    id_col: str,
    genes: List[str],
    label: str,
) -> pd.DataFrame:
    """
    Linear regression: genetic dosage (SNP or haplotype) â†’ gene expression.
    Runs per cohort and globally across all primary tumour samples.
    """
    rows = []
    cohorts = {
        "Breast":       expr_df[expr_df["cohort"] == "Breast"].copy(),
        "Endometrial":  expr_df[expr_df["cohort"] == "Endometrial"].copy(),
        "All_Primary":  expr_df.copy(),
    }
    for cname, cdf in cohorts.items():
        if cdf.empty:
            continue
        codes = set(cdf["snp_code"].dropna())
        sub = genetic_long[genetic_long["snp_code"].isin(codes)].copy()
        for key, gdf in sub.groupby(id_col, dropna=False):
            merge_cols = [c for c in ["analysis_sample_id", "snp_code", "Dosage"] if c in gdf.columns]
            m = cdf.merge(gdf[merge_cols].drop_duplicates([c for c in ["analysis_sample_id"] if c in merge_cols]),
                          on="snp_code", how="inner")
            if m.empty:
                continue
            for gene in genes:
                if gene not in m.columns:
                    continue
                t = m[["Dosage", gene]].apply(
                    pd.to_numeric, errors="coerce").dropna()
                if (len(t) < MIN_N_FOR_ASSOC
                        or t["Dosage"].nunique() < 2
                        or int((t["Dosage"] > 0).sum()) < 3):
                    continue
                slope, inter, r, p, se = stats.linregress(
                    t["Dosage"], t[gene])
                gene_group = next(
                    (grp for grp, genes_in_grp in GENE_GROUPS.items()
                     if gene in genes_in_grp),
                    "Other"
                )
                rows.append({
                    "Analysis_Type":  label,
                    "Cohort":         cname,
                    id_col:           key,
                    "Gene":           gene,
                    "Gene_Group":     gene_group,
                    "N":              len(t),
                    "Carrier_Count":  int((t["Dosage"] > 0).sum()),
                    "Slope":          round(slope, 6),
                    "R":              round(r, 4),
                    "P_Value":        p,
                })

    if not rows:
        return pd.DataFrame()
    result = pd.DataFrame(rows)
    return _add_fdr(result, "P_Value", ["Analysis_Type", "Cohort"])


# â”€â”€ ASSOCIATION: ISOFORM FRACTION vs IMMUNE/PYROPTOSIS GENES â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def assoc_isoform_vs_genes(
    ratio: pd.DataFrame,
    isoform: pd.DataFrame,
    genes: List[str],
) -> pd.DataFrame:
    """
    Spearman correlation between each qPCR isoform endpoint and each
    BAM-derived panel gene. This is the key biological analysis â€” does
    exon-defined isoform balance and individual isoform fractions correlate
    with downstream immune gene expression?

    Runs per cohort and globally.
    """
    if isoform.empty:
        return pd.DataFrame(
            {"Note": ["No isoform data available for correlation analysis"]})

    iso_endpoints = [c for c in ISOFORM_COLS if c in isoform.columns]
    if not iso_endpoints:
        return pd.DataFrame(
            {"Note": ["No isoform columns found in isoform data"]})

    merged = ratio.merge(
        isoform[["snp_code"] + iso_endpoints],
        on="snp_code", how="inner"
    )

    rows = []
    cohorts = {
        "Breast":       merged[merged["cohort"] == "Breast"].copy(),
        "Endometrial":  merged[merged["cohort"] == "Endometrial"].copy(),
        "All_Primary":  merged.copy(),
    }
    for cname, cdf in cohorts.items():
        if cdf.empty:
            continue
        for iso in iso_endpoints:
            for gene in genes:
                if gene not in cdf.columns:
                    continue
                t = cdf[[iso, gene]].apply(
                    pd.to_numeric, errors="coerce").dropna()
                if len(t) < MIN_N_FOR_ASSOC:
                    continue
                rho, p = stats.spearmanr(t[iso], t[gene],
                                         nan_policy="omit")
                gene_group = next(
                    (grp for grp, gg in GENE_GROUPS.items()
                     if gene in gg),
                    "Other"
                )
                rows.append({
                    "Analysis_Type": "Isoform_vs_PanelGene",
                    "Cohort":        cname,
                    "Isoform":       iso,
                    "Gene":          gene,
                    "Gene_Group":    gene_group,
                    "N":             len(t),
                    "Spearman_rho":  round(rho, 4),
                    "P_Value":       p,
                })

    if not rows:
        return pd.DataFrame(
            {"Note": ["No isoform vs gene tests met the minimum thresholds"]})
    result = pd.DataFrame(rows)
    return _add_fdr(result, "P_Value", ["Analysis_Type", "Cohort"])


# â”€â”€ ASSOCIATION: PANEL GENES vs CLINICAL VARIABLES â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def assoc_genes_vs_clinical(
    ratio: pd.DataFrame,
    master: pd.DataFrame,
    genes: List[str],
) -> pd.DataFrame:
    """
    Test whether BAM-derived panel gene expression differs by clinical
    variable. Uses Mann-Whitney U (binary) or Spearman (continuous).
    """
    # Merge expression with clinical master
    master_na = (master[master["nucleic_acid"] == "RNA"]
                  .drop_duplicates("snp_code")
                  .copy())
    merged = ratio.merge(master_na, on="snp_code", how="left")

    rows = []
    for cname, clin_vars in [("Breast", CLINICAL_VARS_BREAST),
                               ("Endometrial", CLINICAL_VARS_ENDO)]:
        cdf = merged[merged["cohort"] == cname].copy()
        if cdf.empty:
            continue
        for gene in genes:
            if gene not in cdf.columns:
                continue
            for cv, (typ, label) in clin_vars.items():
                if cv not in cdf.columns:
                    continue
                t = cdf[[gene, cv]].apply(
                    pd.to_numeric, errors="coerce").dropna()
                if len(t) < MIN_N_FOR_ASSOC:
                    continue
                gene_group = next(
                    (grp for grp, gg in GENE_GROUPS.items()
                     if gene in gg),
                    "Other"
                )
                row = {
                    "Analysis_Type": "PanelGene_vs_Clinical",
                    "Cohort":        cname,
                    "Gene":          gene,
                    "Gene_Group":    gene_group,
                    "Clin_Var":      cv,
                    "Clin_Label":    label,
                    "Clin_Type":     typ,
                    "N":             len(t),
                    "P_Value":       np.nan,
                    "Effect":        np.nan,
                }
                if typ == "continuous":
                    rho, p = stats.spearmanr(t[gene], t[cv],
                                             nan_policy="omit")
                    row.update({"Test": "Spearman", "P_Value": p,
                                "Effect": round(rho, 4)})
                elif typ == "binary":
                    t[cv] = pd.to_numeric(t[cv], errors="coerce")
                    t = t[t[cv].isin([0, 1])]
                    if t[cv].nunique() < 2:
                        continue
                    a = t.loc[t[cv] == 1, gene].dropna()
                    b = t.loc[t[cv] == 0, gene].dropna()
                    if len(a) < 3 or len(b) < 3:
                        continue
                    _, p = stats.mannwhitneyu(a, b, altenative="two-sided")
                    row.update({"Test": "Mann-Whitney U", "P_Value": p,
                                "Effect": round(a.median() - b.median(), 4),
                                "N": len(t)})
                else:
                    continue
                rows.append(row)

    if not rows:
        return pd.DataFrame()
    result = pd.DataFrame(rows)
    return _add_fdr(result, "P_Value", ["Analysis_Type", "Cohort"])


# â”€â”€ GENE GROUP SUMMARY â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def gene_group_summary(ratio: pd.DataFrame, genes: List[str]) -> pd.DataFrame:
    """Median expression per gene per cohort/tissue group."""
    rows = []
    for group_label, gdf in ratio.groupby(
            ["cohort", "tissue", "bam_group"], dropna=False):
        cohort, tissue, group = group_label
        for gene in genes:
            if gene not in gdf.columns:
                continue
            vals = pd.to_numeric(gdf[gene], errors="coerce").dropna()
            gene_group = next(
                (grp for grp, gg in GENE_GROUPS.items() if gene in gg),
                "Other"
            )
            rows.append({
                "cohort":     cohort,
                "tissue":     tissue,
                "bam_group":  group,
                "gene":       gene,
                "gene_group": gene_group,
                "n_samples":  len(vals),
                "median":     round(vals.median(), 4) if len(vals) else np.nan,
                "mean":       round(vals.mean(), 4) if len(vals) else np.nan,
                "std":        round(vals.std(), 4) if len(vals) else np.nan,
                "pct_zero":   round((vals == 0).mean() * 100, 1) if len(vals) else np.nan,
            })
    return pd.DataFrame(rows).sort_values(
        ["gene_group", "gene", "cohort", "tissue"])


# â”€â”€ FIGURES â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def _style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.grid(True, linestyle=":", color="#cccccc", alpha=0.5)
    ax.set_axisbelow(True)


def plot_crossval(merged_cv: pd.DataFrame, out_dir: Path) -> None:
    """Scatter plot: BAM GSDMB vs qPCR GSDMB with Spearman rho annotation."""
    if merged_cv.empty or "GSDMB_BAM" not in merged_cv.columns:
        return

    cohorts = merged_cv["cohort"].dropna().unique()
    colors = {"Breast": "#AD1457", "Endometrial": "#00695C"}
    fig, ax = plt.subplots(figsize=(7, 6))

    for cohort in cohorts:
        sub = merged_cv[merged_cv["cohort"] == cohort].dropna(
            subset=["GSDMB_BAM", "GSDMB_qPCR"])
        if sub.empty:
            continue
        ax.scatter(sub["GSDMB_BAM"], sub["GSDMB_qPCR"],
                   label=cohort, color=colors.get(cohort, "#888888"),
                   alpha=0.75, s=55, edgecolors="white", linewidths=0.5)

    # Overall correlation line
    t = merged_cv[["GSDMB_BAM", "GSDMB_qPCR"]].apply(
        pd.to_numeric, errors="coerce").dropna()
    if len(t) >= 5:
        rho, p = stats.spearmanr(t["GSDMB_BAM"], t["GSDMB_qPCR"])
        m, b = np.polyfit(t["GSDMB_BAM"], t["GSDMB_qPCR"], 1)
        xr = np.linspace(t["GSDMB_BAM"].min(), t["GSDMB_BAM"].max(), 100)
        ax.plot(xr, m * xr + b, color="#333333", linewidth=1.5,
                linestyle="--", alpha=0.7)
        ax.text(0.05, 0.95,
                f"Spearman Ï = {rho:.3f}\np = {p:.3e} (n={len(t)})",
                transform=ax.transAxes, va="top", ha="left",
                fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="#cccccc", alpha=0.9))

    ax.set_xlabel("BAM-derived GSDMB (ratio-normalised depth)", fontsize=10)
    ax.set_ylabel("qPCR GSDMB total expression", fontsize=10)
    ax.set_title("Cross-validation: BAM vs qPCR GSDMB expression",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    _style_ax(ax)
    plt.tight_layout()
    fig.savefig(out_dir / "20b_GSDMB_CrossValidation.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: 20b_GSDMB_CrossValidation.png")


def plot_heatmap(df: pd.DataFrame, row_col: str, col_col: str,
                 val_col: str, title: str, out_path: Path,
                 max_rows: int = 30, max_cols: int = 20) -> None:
    """Generic signed -log10(p) heatmap."""
    if df.empty or val_col not in df.columns:
        return
    d = df.copy()
    d["_lp"] = -np.log10(pd.to_numeric(d[val_col],
                                        errors="coerce").clip(lower=1e-300))
    effect_col = next((c for c in ["Spearman_rho", "Slope", "R", "Effect"]
                       if c in d.columns), None)
    sign = (pd.to_numeric(d[effect_col], errors="coerce").fillna(0)
            if effect_col else pd.Series(1, index=d.index))
    d["_signed"] = d["_lp"] * np.sign(sign)

    top_rows = (d.groupby(row_col)["_lp"].max()
                .sort_values(ascending=False)
                .head(max_rows).index)
    top_cols = (d.groupby(col_col)["_lp"].max()
                .sort_values(ascending=False)
                .head(max_cols).index)
    d = d[d[row_col].isin(top_rows) & d[col_col].isin(top_cols)]
    pivot = d.pivot_table(index=row_col, columns=col_col,
                          values="_signed", aggfunc="max")
    if pivot.empty:
        return

    fig_w = max(10, 0.7 * len(pivot.columns) + 4)
    fig_h = max(4,  0.4 * len(pivot.index)   + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    if _HAS_SNS:
        sns.heatmap(pivot, cmap="coolwarm", center=0,
                    linewidths=0.4, linecolor="white", ax=ax,
                    cbar_kws={"label": "signed -log10(P value)"})
    else:
        im = ax.imshow(pivot.values, aspect="auto", cmap="coolwarm")
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_yticks(range(len(pivot.index)))
        ax.set_xticklabels(pivot.columns, rotation=40, ha="right", fontsize=8)
        ax.set_yticklabels(pivot.index, fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02,
                     label="signed -log10(P value)")

    ax.set_title(title, fontsize=12, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


def plot_gene_group_distributions(
    ratio: pd.DataFrame,
    genes: List[str],
    out_dir: Path,
) -> None:
    """
    One subplot per gene group showing per-cohort expression distributions
    as box + strip plots.
    """
    group_order = ["Breast_Tumour", "Endometrial_Tumour",
                   "Breast_Normal", "Endometrial_Normal"]
    colors = {
        "Breast_Tumour":       "#AD1457",
        "Endometrial_Tumour":  "#00695C",
        "Breast_Normal":       "#F48FB1",
        "Endometrial_Normal":  "#80CBC4",
    }
    group_names = [g for g in GENE_GROUPS.keys()
                   if g != "Housekeeping"]

    fig, axes = plt.subplots(
        2, 2, figsize=(18, 12), constrained_layout=True)
    axes_flat = axes.flatten()

    for ax, group_name in zip(axes_flat, group_names):
        group_genes = [g for g in GENE_GROUPS[group_name] if g in genes]
        if not group_genes:
            ax.set_visible(False)
            continue

        # Melt to long format for this gene group
        plot_cols = ["bam_group"] + group_genes
        plot_cols = [c for c in plot_cols if c in ratio.columns]
        melt = ratio[plot_cols].melt(
            id_vars="bam_group", var_name="gene", value_name="expression")
        melt["expression"] = pd.to_numeric(
            melt["expression"], errors="coerce")
        melt = melt.dropna()
        melt["bam_group"] = pd.Categorical(
            melt["bam_group"], categories=group_order, ordered=True)

        if _HAS_SNS:
            sns.boxplot(data=melt, x="gene", y="expression",
                        hue="bam_group", hue_order=group_order,
                        palette=colors, ax=ax, fliersize=0,
                        width=0.7, linewidth=0.8)
            ax.legend(title="Cohort", fontsize=7,
                      title_fontsize=7, loc="upper right")
        else:
            # Fallback: one box per gene, all groups combined
            gene_data = [
                melt.loc[melt["gene"] == g, "expression"].values
                for g in group_genes
            ]
            ax.boxplot(gene_data, labels=group_genes,
                       patch_artist=True,
                       boxprops=dict(facecolor="#d9e8f5"))

        ax.set_title(f"{group_name.replace('_', ' ')} genes",
                     fontweight="bold", fontsize=11)
        ax.set_xlabel("")
        ax.set_xticklabels(group_genes, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Normalised expression", fontsize=9)
        ax.set_yscale("symlog", linthresh=0.1)
        _style_ax(ax)

    fig.suptitle("Panel gene expression by biological group and cohort",
                 fontsize=14, fontweight="bold")
    fig.savefig(out_dir / "20b_Gene_Group_Distributions.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: 20b_Gene_Group_Distributions.png")


def plot_chr17_locus(ratio: pd.DataFrame, out_dir: Path) -> None:
    """
    Bar chart of median expression for chr17q12 locus genes per cohort.
    Highlights GSDMB.
    """
    locus_genes = [g for g in GENE_GROUPS["Chr17q12_Locus"]
                   if g in ratio.columns]
    if not locus_genes:
        return

    groups = ratio["bam_group"].dropna().unique()
    colors = {
        "Breast_Tumour":      "#AD1457",
        "Endometrial_Tumour": "#00695C",
        "Breast_Normal":      "#F48FB1",
        "Endometrial_Normal": "#80CBC4",
    }
    x = np.arange(len(locus_genes))
    width = 0.8 / max(len(groups), 1)

    fig, ax = plt.subplots(figsize=(14, 5))
    for i, group in enumerate(sorted(groups)):
        sub = ratio[ratio["bam_group"] == group]
        medians = [
            pd.to_numeric(sub[g], errors="coerce").median()
            for g in locus_genes
        ]
        offset = (i - len(groups) / 2) * width + width / 2
        bars = ax.bar(x + offset, medians, width * 0.9,
                      label=group,
                      color=colors.get(group, "#888888"),
                      alpha=0.85, edgecolor="white", linewidth=0.5)

    # Highlight GSDMB
    if "GSDMB" in locus_genes:
        gsdmb_idx = locus_genes.index("GSDMB")
        ax.axvspan(gsdmb_idx - 0.5, gsdmb_idx + 0.5,
                   alpha=0.08, color="#c62828", zorder=0)
        ax.text(gsdmb_idx, ax.get_ylim()[1] * 0.98, "GSDMB",
                ha="center", va="top", fontsize=8,
                color="#c62828", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(locus_genes, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Median normalised expression", fontsize=10)
    ax.set_title("Chr17q12 locus gene expression by cohort",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right", title="Cohort")
    _style_ax(ax)
    plt.tight_layout()
    fig.savefig(out_dir / "20b_Chr17_Locus_Expression.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: 20b_Chr17_Locus_Expression.png")


def plot_expression_heatmap(ratio: pd.DataFrame,
                            genes: List[str],
                            out_dir: Path) -> None:
    """
    Sample Ã— gene expression heatmap, log1p scaled, samples grouped by cohort.
    """
    if ratio.empty:
        return
    gene_cols = [g for g in genes if g in ratio.columns]
    if not gene_cols:
        return

    mat = ratio[gene_cols].apply(pd.to_numeric, errors="coerce")
    mat = np.log1p(mat.fillna(0))

    fig_w = max(14, 0.25 * len(gene_cols) + 4)
    fig_h = max(6,  0.12 * len(ratio)      + 2)
    fig_h = min(fig_h, 30)  # cap height

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    if _HAS_SNS:
        # Row colors by group
        row_colors = ratio["bam_group"].map({
            "Breast_Tumour":      "#AD1457",
            "Endometrial_Tumour": "#00695C",
            "Breast_Normal":      "#F48FB1",
            "Endometrial_Normal": "#80CBC4",
        }).fillna("#888888")
        sns.heatmap(mat, cmap="viridis", ax=ax,
                    xticklabels=gene_cols,
                    yticklabels=False,
                    cbar_kws={"label": "log1p(normalised expression)"})
        ax.set_xticklabels(gene_cols, rotation=45, ha="right", fontsize=7)
    else:
        im = ax.imshow(mat.values, aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(gene_cols)))
        ax.set_xticklabels(gene_cols, rotation=45, ha="right", fontsize=7)
        ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02,
                     label="log1p(normalised expression)")

    ax.set_title("Panel gene expression heatmap (all samples)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out_dir / "20b_Expression_Heatmap.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: 20b_Expression_Heatmap.png")


# â”€â”€ MAIN â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="20b: BAM-derived panel gene quantification and association",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--project-root",   default=str(DEFAULTS["project_root"]))
    ap.add_argument("--rna-bed",        default=str(DEFAULTS["rna_bed"]))
    ap.add_argument("--expr-xlsx",      default=str(DEFAULTS["expr_xlsx"]))
    ap.add_argument("--master",         default=str(DEFAULTS["master"]))
    ap.add_argument("--variant-workbook", default=str(DEFAULTS["variant_workbook"]))
    ap.add_argument("--haplotype-input",  default=str(DEFAULTS["haplotype_input"]))
    ap.add_argument("--rna-qc-manifest",  default=str(DEFAULTS["rna_qc_manifest"]))
    ap.add_argument("--out-dir",        default=str(DEFAULTS["out_dir"]))
    ap.add_argument("--recompute",      action="store_true",
                    help="Ignore cached raw counts and rerun samtools bedcov")
    ap.add_argument("--min-nfe-af",     type=float, default=0.01)
    ap.add_argument("--min-hap-freq",   type=float, default=0.02)
    ap.add_argument("--strict-rna-qc", action="store_true",
                    help="Restrict BAM-derived quantification to strict analysis-ready RNA samples only.")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root)
    rna_bed     = Path(args.rna_bed)
    expr_path    = Path(args.expr_xlsx)
    master_path  = Path(args.master)
    var_path     = Path(args.variant_workbook)
    hap_path     = Path(args.haplotype_input)
    manifest_path = Path(args.rna_qc_manifest)
    out_dir      = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Delete cache if recompute requested
    if args.recompute:
        cache = out_dir / "20b_raw_counts_cache.tsv"
        if cache.exists():
            cache.unlink()
            print("  Cache deleted; will rerun samtools bedcov")

    print("\n" + "=" * 70)
    print("SCRIPT 20b: BAM-DERIVED PANEL GENE EXPRESSION ANALYSIS")
    print("=" * 70)

    # â”€â”€ 1. Panel design â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    print("\n[1/9] Loading panel design...")
    if not rna_bed.exists():
        raise FileNotFoundError(f"RNA BED not found: {rna_bed}")
    panel = load_panel_design(rna_bed)
    genes = panel["gene"].tolist()

    # â”€â”€ 2. QC manifest â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    print("\n[2/9] Loading RNA QC manifest...")
    qc_manifest = None
    if manifest_path.exists():
        qc_manifest = pd.read_csv(manifest_path, sep="\t")
        print(f"  Loaded manifest: {len(qc_manifest)} rows")
    else:
        print("  No QC manifest found; all BAMs will be processed")

    # â”€â”€ 3. Discover BAMs â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    print("\n[3/9] Discovering RNA BAMs...")
    bam_df = discover_bams(project_root, qc_manifest, strict_only=args.strict_rna_qc)
    if bam_df.empty:
        print("  No BAMs found. Exiting.")
        return

    # â”€â”€ 4. Quantify â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    print("\n[4/9] Quantifying gene expression (samtools bedcov)...")
    raw = quantify_all_bams(bam_df, panel, out_dir)

    # â”€â”€ 5. Normalise â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    print("\n[5/9] Normalising expression...")
    cpm, ratio = normalise(raw, panel)

    bam_usage_summary = build_bam_usage_summary(raw)
    normalisation_audit = build_normalisation_audit(raw, args.strict_rna_qc)

    # Restrict to primary tumour samples for all association analyses
    primary_ratio = ratio[ratio["analysis_role"] == "primary_tumour"].copy()
    print(f"  Primary tumour samples for association analyses: "
          f"{len(primary_ratio)}")

    # â”€â”€ 6. Load supporting data â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    print("\n[6/9] Loading isoform, genetic, and clinical data...")
    isoform = load_isoform_data(expr_path) if expr_path.exists() else pd.DataFrame()
    snp_long, hap_long, hap_freqs = load_genetic_data(
        var_path, hap_path,
        min_nfe=args.min_nfe_af,
        min_hap_freq=args.min_hap_freq,
    ) if var_path.exists() and hap_path.exists() else (
        pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    )
    master = pd.DataFrame()
    if master_path.exists():
        try:
            master = pd.read_excel(master_path,
                                   sheet_name="harmonised_plus_canon")
            master["snp_code"] = master["snp_code"].astype(str).str.strip()
            master["nucleic_acid"] = master["nucleic_acid"].astype(str).str.upper()
            print(f"  Master loaded: {len(master)} rows")
        except Exception as e:
            print(f"  [WARN] Could not load master: {e}")

    # â”€â”€ 7. Analyses â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    print("\n[7/9] Running association analyses...")

    # Cross-validation
    cv_merged = pd.DataFrame()
    cv_stats  = pd.DataFrame({"Note": ["No cross-validation data available"]})
    if not isoform.empty and "GSDMB" in isoform.columns:
        cv_result = cross_validate_gsdmb(primary_ratio, isoform)
        if isinstance(cv_result, tuple):
            cv_merged, cv_stats = cv_result
        else:
            cv_stats = cv_result
        print(f"  Cross-validation: {len(cv_stats)} group results")

    # SNP vs panel genes
    snp_gene = pd.DataFrame()
    if not snp_long.empty:
        snp_gene = assoc_genetic_vs_genes(
            primary_ratio, snp_long, "Variant_ID", genes,
            "SNP_vs_PanelGene")
        print(f"  SNP vs panel gene: {len(snp_gene)} tests, "
              f"{int(snp_gene['Nominal_Sig'].sum()) if not snp_gene.empty else 0} nominal")

    # Haplotype vs panel genes
    hap_gene = pd.DataFrame()
    if not hap_long.empty:
        hap_gene = assoc_genetic_vs_genes(
            primary_ratio, hap_long, "Haplotype_ID", genes,
            "Haplotype_vs_PanelGene")
        print(f"  Haplotype vs panel gene: {len(hap_gene)} tests, "
              f"{int(hap_gene['Nominal_Sig'].sum()) if not hap_gene.empty else 0} nominal")

    # Isoform fractions vs panel genes (key biological analysis)
    isoform_gene = pd.DataFrame()
    if not isoform.empty:
        isoform_gene = assoc_isoform_vs_genes(
            primary_ratio, isoform, genes)
        n_tests = len(isoform_gene) if not isoform_gene.empty and "P_Value" in isoform_gene.columns else 0
        n_nom   = int(isoform_gene["Nominal_Sig"].sum()) if n_tests else 0
        print(f"  Isoform vs panel gene: {n_tests} tests, {n_nom} nominal")

    # Panel genes vs clinical
    gene_clin = pd.DataFrame()
    if not master.empty:
        gene_clin = assoc_genes_vs_clinical(
            primary_ratio, master, genes)
        n_tests = len(gene_clin) if not gene_clin.empty and "P_Value" in gene_clin.columns else 0
        n_nom   = int(gene_clin["Nominal_Sig"].sum()) if n_tests else 0
        print(f"  Panel gene vs clinical: {n_tests} tests, {n_nom} nominal")

    # Gene group summary
    summary = gene_group_summary(ratio, genes)

    # â”€â”€ 8. Save outputs â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    print("\n[8/9] Saving results workbook...")
    xlsx_path = out_dir / "20b_Panel_Gene_Expression.xlsx"

    def _safe(df, note="No data available"):
        return df if (not df.empty and "Note" not in df.columns) else pd.DataFrame({"Note": [note]})

    sample_manifest = raw[[c for c in META_COLS if c in raw.columns]].copy()

    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        raw.to_excel(writer, sheet_name="raw_counts", index=False)
        cpm.to_excel(writer, sheet_name="normalised_cpm", index=False)
        ratio.to_excel(writer, sheet_name="normalised_ratio", index=False)
        sample_manifest.to_excel(writer, sheet_name="sample_manifest", index=False)
        bam_usage_summary.to_excel(writer, sheet_name="bam_usage_summary", index=False)
        normalisation_audit.to_excel(writer, sheet_name="normalisation_audit", index=False)
        summary.to_excel(writer, sheet_name="gene_group_summary", index=False)
        panel.to_excel(writer, sheet_name="panel_design", index=False)
        cv_stats.to_excel(writer, sheet_name="gsdmb_crossval", index=False)
        (snp_gene if not snp_gene.empty else
         pd.DataFrame({"Note": ["No SNP vs gene tests met thresholds"]}
                      )).to_excel(writer, sheet_name="snp_gene_assoc", index=False)
        (hap_gene if not hap_gene.empty else
         pd.DataFrame({"Note": ["No haplotype vs gene tests met thresholds"]}
                      )).to_excel(writer, sheet_name="haplotype_gene_assoc", index=False)
        (isoform_gene if (not isoform_gene.empty
                          and "P_Value" in isoform_gene.columns) else
         pd.DataFrame({"Note": ["No isoform vs gene tests met thresholds"]}
                      )).to_excel(writer, sheet_name="isoform_immune_corr", index=False)
        (gene_clin if (not gene_clin.empty
                       and "P_Value" in gene_clin.columns) else
         pd.DataFrame({"Note": ["No gene vs clinical tests met thresholds"]}
                      )).to_excel(writer, sheet_name="clinical_gene_assoc", index=False)
        if not hap_freqs.empty:
            hap_freqs.to_excel(writer, sheet_name="haplotype_freqs", index=False)

    print(f"  Saved: {xlsx_path}")

    # â”€â”€ 9. Figures â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    print("\n[9/9] Generating figures...")

    if not cv_merged.empty:
        plot_crossval(cv_merged, out_dir)

    plot_expression_heatmap(ratio, genes, out_dir)
    plot_gene_group_distributions(ratio, genes, out_dir)
    plot_chr17_locus(ratio, out_dir)

    if not snp_gene.empty and "P_Value" in snp_gene.columns:
        plot_heatmap(
            snp_gene, "Variant_ID", "Gene", "P_Value",
            "SNP vs panel gene associations",
            out_dir / "20b_SNP_PanelGene_Heatmap.png",
        )
    if not hap_gene.empty and "P_Value" in hap_gene.columns:
        plot_heatmap(
            hap_gene, "Haplotype_ID", "Gene", "P_Value",
            "Haplotype vs panel gene associations",
            out_dir / "20b_Haplotype_PanelGene_Heatmap.png",
        )
    if (not isoform_gene.empty
            and "P_Value" in isoform_gene.columns):
        plot_heatmap(
            isoform_gene, "Isoform", "Gene", "P_Value",
            "qPCR isoform fraction vs BAM panel gene expression\n"
            "(signed Spearman rho, -log10(P value) scaled)",
            out_dir / "20b_Isoform_Immune_Heatmap.png",
        )

    # â”€â”€ Summary â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    print("\n" + "=" * 70)
    print("SCRIPT 20b COMPLETE")
    print("=" * 70)
    print(f"  RNA BAMs processed        : {len(raw)}")
    print(f"  Genes quantified          : {len(genes)}")
    print(f"  Primary tumour samples    : {len(primary_ratio)}")
    if "RNA_BAM_Usage_Class" in raw.columns:
        usage_counts = raw["RNA_BAM_Usage_Class"].value_counts(dropna=False)
        print(f"  Strict BAM rows           : {int(usage_counts.get('Strict', 0))}")
        print(f"  Exploratory BAM rows      : {int(usage_counts.get('Exploratory', 0))}")
    print(f"    Breast                  : {int((primary_ratio['cohort'] == 'Breast').sum())}")
    print(f"    Endometrial             : {int((primary_ratio['cohort'] == 'Endometrial').sum())}")
    if not snp_gene.empty and "Nominal_Sig" in snp_gene.columns:
        print(f"  SNP vs gene (nominal)     : {int(snp_gene['Nominal_Sig'].sum())}")
        print(f"  SNP vs gene (FDR<{FDR_THRESHOLD:.0%})    : {int(snp_gene['FDR_Sig'].sum())}")
    if not hap_gene.empty and "Nominal_Sig" in hap_gene.columns:
        print(f"  Haplotype vs gene (nom.)  : {int(hap_gene['Nominal_Sig'].sum())}")
        print(f"  Haplotype vs gene (FDR)   : {int(hap_gene['FDR_Sig'].sum())}")
    if (not isoform_gene.empty
            and "Nominal_Sig" in isoform_gene.columns):
        print(f"  Isoform vs gene (nom.)    : {int(isoform_gene['Nominal_Sig'].sum())}")
        print(f"  Isoform vs gene (FDR)     : {int(isoform_gene['FDR_Sig'].sum())}")
    if not gene_clin.empty and "Nominal_Sig" in gene_clin.columns:
        print(f"  Gene vs clinical (nom.)   : {int(gene_clin['Nominal_Sig'].sum())}")
        print(f"  Gene vs clinical (FDR)    : {int(gene_clin['FDR_Sig'].sum())}")
    print(f"  Results workbook          : {xlsx_path}")
    print(f"  Output directory          : {out_dir}")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()







