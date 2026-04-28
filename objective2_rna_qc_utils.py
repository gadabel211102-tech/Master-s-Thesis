from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline_utils import get_thresholds

try:
    import pysam
except ImportError:
    pysam = None

try:
    import seaborn as sns
except ImportError:
    sns = None

RNA_THRESHOLDS = get_thresholds()
RNA_PANEL_LOW_DEPTH = float(RNA_THRESHOLDS.get("rna_panel_depth_target", 20.0))
RNA_PANEL_EXPLORATORY_TARGET_PCT = float(RNA_THRESHOLDS.get("rna_panel_exploratory_target_pct", 0.50))
RNA_PANEL_STRICT_TARGET_PCT = float(RNA_THRESHOLDS.get("rna_panel_strict_target_pct", 0.70))
RNA_GENE_EXPLORATORY_MIN_READS = int(RNA_THRESHOLDS.get("rna_gene_exploratory_min_reads", 75))
RNA_GENE_STRICT_MIN_READS = int(RNA_THRESHOLDS.get("rna_gene_strict_min_reads", 100))
RNA_GENE_EXPLORATORY_TARGET_PCT = float(RNA_THRESHOLDS.get("rna_gene_exploratory_target_pct", 0.60))
RNA_GENE_STRICT_TARGET_PCT = float(RNA_THRESHOLDS.get("rna_gene_strict_target_pct", 0.70))
RNA_GENE_THRESHOLD_SCOPE = str(RNA_THRESHOLDS.get("rna_gene_threshold_scope", "all_panel"))
RNA_MACHINE_HARD_FAIL_READS = int(RNA_THRESHOLDS.get("rna_machine_hard_fail_reads", 10000))
RNA_MACHINE_WARN_READS = int(RNA_THRESHOLDS.get("rna_machine_warn_reads", 50000))
RNA_MACHINE_WARN_Q20_PCT = float(RNA_THRESHOLDS.get("rna_machine_warn_q20_pct", 92.0))
RNA_MACHINE_HARD_FAIL_Q20_PCT = float(RNA_THRESHOLDS.get("rna_machine_hard_fail_q20_pct", 85.0))
DNA_BASELINE_COLUMNS = [
    "DNA_BAM_Group",
    "DNA_BAM_Path",
    "DNA_BAM_Name",
    "DNA_BAM_Size_MB",
    "DNA_Baseline_Present",
]
RNA_QC_COLUMNS = [
    "RNA_BAM_Group",
    "RNA_BAM_Path",
    "RNA_BAM_Name",
    "RNA_BAM_Size_MB",
    "RNA_BAI_Present",
    "RNA_BAI_Created",
    "RNA_Coverage_QC_Status",
    "RNA_Panel_Target_Count",
    "RNA_Panel_Targets_Covered",
    "RNA_Panel_Target_Coverage_Pct",
    "RNA_Panel_Targets_GE20",
    "RNA_Panel_Targets_GE20_Pct",
    "RNA_Panel_Low_Depth_Targets",
    "RNA_Panel_Mean_Depth",
    "RNA_Panel_Median_Depth",
    "RNA_Panel_Targets_With_Read_Counts",
    "RNA_Panel_Targets_GE75_Reads",
    "RNA_Panel_Targets_GE75_Reads_Pct",
    "RNA_Panel_Targets_GE100_Reads",
    "RNA_Panel_Targets_GE100_Reads_Pct",
    "RNA_Panel_Min_Target_Reads",
    "RNA_GSDMB_Target_Count",
    "RNA_GSDMB_Targets_Covered",
    "RNA_GSDMB_Mean_Depth",
    "RNA_GSDMB_Min_Depth",
    "RNA_Machine_Reads",
    "RNA_Machine_Q20_Pct",
    "RNA_Machine_QC_Flag",
    "RNA_Coverage_Risk_Flag",
    "RNA_Coverage_Risk_Level",
    "RNA_Coverage_Risk_Note",
]


def load_rna_panel_design(path: str | Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("track") or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            transcript, start, end, amplicon_id, score, gene_symbol = parts[:6]
            rows.append(
                {
                    "transcript": transcript,
                    "start": int(start),
                    "end": int(end),
                    "amplicon_id": amplicon_id,
                    "score": score,
                    "gene_symbol": gene_symbol,
                }
            )
    design = pd.DataFrame(rows)
    if design.empty:
        return design
    design["target_length"] = design["end"] - design["start"]
    design["region_label"] = design["transcript"] + ":" + design["start"].astype(str) + "-" + design["end"].astype(str)
    return design


def ensure_bam_index(bam_path: str | Path) -> tuple[bool, bool, str]:
    bam = Path(bam_path)
    index_path = bam.with_suffix(".bam.bai")
    if index_path.exists():
        return True, False, "Index already present"
    if shutil.which("samtools") is None:
        return False, False, "samtools not available for BAM indexing"
    proc = subprocess.run(["samtools", "index", "-@", "2", str(bam)], capture_output=True, text=True)
    if proc.returncode != 0:
        message = proc.stderr.strip() or proc.stdout.strip() or "samtools index failed"
        return False, False, message
    return index_path.exists(), True, "Index created"


def run_bedcov(bed_path: str | Path, bam_path: str | Path, panel_design: pd.DataFrame) -> pd.DataFrame:
    proc = subprocess.run(["samtools", "bedcov", str(bed_path), str(bam_path)], capture_output=True, text=True)
    if proc.returncode != 0:
        message = proc.stderr.strip() or proc.stdout.strip() or "samtools bedcov failed"
        raise RuntimeError(message)
    coverage_values: list[float] = []
    for line in proc.stdout.splitlines():
        if not line.strip() or line.startswith("track"):
            continue
        parts = line.split("	")
        coverage_values.append(float(parts[-1]))
    if len(coverage_values) != len(panel_design):
        raise RuntimeError(f"bedcov returned {len(coverage_values)} rows for {len(panel_design)} design targets")
    out = panel_design.copy()
    out["covered_base_sum"] = coverage_values
    out["mean_depth"] = np.where(out["target_length"] > 0, out["covered_base_sum"] / out["target_length"], np.nan)
    return out


def run_target_read_counts(bam_path: str | Path, panel_design: pd.DataFrame) -> list[int]:
    counts: list[int] = []
    bam = str(bam_path)
    if pysam is not None:
        try:
            with pysam.AlignmentFile(bam, "rb") as handle:
                for row in panel_design.itertuples(index=False):
                    counts.append(int(handle.count(contig=str(row.transcript), start=int(row.start), end=int(row.end), read_callback="all")))
            return counts
        except Exception as exc:
            raise RuntimeError(f"pysam read counting failed: {exc}") from exc
    for row in panel_design.itertuples(index=False):
        region = f"{row.transcript}:{int(row.start) + 1}-{int(row.end)}"
        proc = subprocess.run(["samtools", "view", "-c", bam, region], capture_output=True, text=True)
        if proc.returncode != 0:
            message = proc.stderr.strip() or proc.stdout.strip() or f"samtools view -c failed for {region}"
            raise RuntimeError(message)
        try:
            counts.append(int((proc.stdout or "0").strip() or "0"))
        except ValueError as exc:
            raise RuntimeError(f"Unexpected read-count output for {region}: {proc.stdout!r}") from exc
    return counts

def classify_rna_coverage(
    panel_ge20_pct: float,
    panel_target_count: int,
    targets_with_read_counts: int,
    panel_targets_ge75_pct: float,
    panel_targets_ge100_pct: float,
) -> tuple[bool, str, str]:
    if panel_target_count <= 0:
        return True, "Unavailable", "The RNA panel design did not yield any measurable panel-gene targets to validate."
    if targets_with_read_counts <= 0:
        return True, "High", "No measurable per-gene read support was detected across the RNA panel targets."
    if targets_with_read_counts < panel_target_count:
        return True, "High", "Per-gene read counts were missing for part of the RNA panel target set; interpret these samples very cautiously."
    if pd.isna(panel_ge20_pct) or panel_ge20_pct < RNA_PANEL_EXPLORATORY_TARGET_PCT:
        return True, "High", f"Too few panel targets reached {RNA_PANEL_LOW_DEPTH:.0f}x depth ({panel_ge20_pct:.1%} < {RNA_PANEL_EXPLORATORY_TARGET_PCT:.0%})."
    if pd.isna(panel_targets_ge75_pct) or panel_targets_ge75_pct < RNA_GENE_EXPLORATORY_TARGET_PCT:
        return True, "High", (
            f"Too few panel targets reached {RNA_GENE_EXPLORATORY_MIN_READS} reads "
            f"({panel_targets_ge75_pct:.1%} < {RNA_GENE_EXPLORATORY_TARGET_PCT:.0%})."
        )
    if pd.isna(panel_targets_ge100_pct) or panel_targets_ge100_pct < RNA_GENE_STRICT_TARGET_PCT or panel_ge20_pct < RNA_PANEL_STRICT_TARGET_PCT:
        return True, "Moderate", (
            f"RNA coverage was usable but below the strict gate: {panel_ge20_pct:.1%} of panel targets reached "
            f"{RNA_PANEL_LOW_DEPTH:.0f}x and {panel_targets_ge100_pct:.1%} reached {RNA_GENE_STRICT_MIN_READS} reads."
        )
    return False, "Low", "RNA panel target coverage and panel-gene read uniformity passed the strict QC gate."


def summarise_bam_coverage(target_cov: pd.DataFrame) -> dict[str, object]:
    panel_target_count = len(target_cov)
    panel_depths = pd.to_numeric(target_cov["mean_depth"], errors="coerce")
    panel_targets_covered = int((panel_depths > 0).sum())
    panel_pct = panel_targets_covered / panel_target_count if panel_target_count else np.nan
    panel_targets_ge20 = int((panel_depths >= RNA_PANEL_LOW_DEPTH).sum())
    panel_ge20_pct = panel_targets_ge20 / panel_target_count if panel_target_count else np.nan
    panel_mean_depth = panel_depths.mean()
    panel_median_depth = panel_depths.median()
    panel_low_depth_targets = int((panel_depths < RNA_PANEL_LOW_DEPTH).sum())

    read_counts = pd.to_numeric(target_cov.get("target_read_count"), errors="coerce")
    panel_targets_with_read_counts = int(read_counts.notna().sum())
    panel_targets_ge75_reads = int((read_counts >= RNA_GENE_EXPLORATORY_MIN_READS).sum())
    panel_targets_ge75_pct = panel_targets_ge75_reads / panel_target_count if panel_target_count else np.nan
    panel_targets_ge100_reads = int((read_counts >= RNA_GENE_STRICT_MIN_READS).sum())
    panel_targets_ge100_pct = panel_targets_ge100_reads / panel_target_count if panel_target_count else np.nan
    panel_min_target_reads = read_counts.min() if panel_targets_with_read_counts else np.nan

    gsdmb_cov = target_cov[target_cov["gene_symbol"].astype(str).str.upper() == "GSDMB"].copy()
    gsdmb_target_count = len(gsdmb_cov)
    gsdmb_targets_covered = int((pd.to_numeric(gsdmb_cov.get("mean_depth"), errors="coerce") > 0).sum()) if gsdmb_target_count else 0
    gsdmb_mean_depth = pd.to_numeric(gsdmb_cov.get("mean_depth"), errors="coerce").mean() if gsdmb_target_count else np.nan
    gsdmb_min_depth = pd.to_numeric(gsdmb_cov.get("mean_depth"), errors="coerce").min() if gsdmb_target_count else np.nan

    risk_flag, risk_level, risk_note = classify_rna_coverage(
        panel_ge20_pct,
        panel_target_count,
        panel_targets_with_read_counts,
        panel_targets_ge75_pct,
        panel_targets_ge100_pct,
    )
    return {
        "RNA_Coverage_QC_Status": "Complete",
        "RNA_Panel_Target_Count": panel_target_count,
        "RNA_Panel_Targets_Covered": panel_targets_covered,
        "RNA_Panel_Target_Coverage_Pct": panel_pct,
        "RNA_Panel_Targets_GE20": panel_targets_ge20,
        "RNA_Panel_Targets_GE20_Pct": panel_ge20_pct,
        "RNA_Panel_Low_Depth_Targets": panel_low_depth_targets,
        "RNA_Panel_Mean_Depth": panel_mean_depth,
        "RNA_Panel_Median_Depth": panel_median_depth,
        "RNA_Panel_Targets_With_Read_Counts": panel_targets_with_read_counts,
        "RNA_Panel_Targets_GE75_Reads": panel_targets_ge75_reads,
        "RNA_Panel_Targets_GE75_Reads_Pct": panel_targets_ge75_pct,
        "RNA_Panel_Targets_GE100_Reads": panel_targets_ge100_reads,
        "RNA_Panel_Targets_GE100_Reads_Pct": panel_targets_ge100_pct,
        "RNA_Panel_Min_Target_Reads": panel_min_target_reads,
        "RNA_GSDMB_Target_Count": gsdmb_target_count,
        "RNA_GSDMB_Targets_Covered": gsdmb_targets_covered,
        "RNA_GSDMB_Mean_Depth": gsdmb_mean_depth,
        "RNA_GSDMB_Min_Depth": gsdmb_min_depth,
        "RNA_Coverage_Risk_Flag": risk_flag,
        "RNA_Coverage_Risk_Level": risk_level,
        "RNA_Coverage_Risk_Note": risk_note,
    }


def bam_validation(root: Path, expr: pd.DataFrame, rna_bed: str | Path | None, sx_func) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    roots = [
        ("Breast_Tumour", root / "breast" / "tumour" / "rna"),
        ("Breast_Normal", root / "breast" / "normal" / "rna"),
        ("Endometrial_Tumour", root / "endometrium" / "tumour" / "rna"),
        ("Endometrial_Normal", root / "endometrium" / "normal" / "rna"),
    ]
    expr_codes = set(expr["snp_code"].dropna())
    primary_codes = set(expr.loc[expr["analysis_include_primary"], "snp_code"].dropna())
    rows: list[dict[str, object]] = []
    for label, folder in roots:
        if not folder.exists():
            continue
        for bam in sorted(folder.glob("*.bam")):
            code = sx_func(bam.stem)
            rows.append(
                {
                    "RNA_BAM_Group": label,
                    "RNA_BAM_Path": str(bam),
                    "RNA_BAM_Name": bam.name,
                    "RNA_BAM_Size_MB": round(bam.stat().st_size / (1024 * 1024), 2),
                    "RNA_BAI_Present": bam.with_suffix(".bam.bai").exists(),
                    "RNA_BAI_Created": False,
                    "snp_code": code,
                    "matched_expression_sample": bool(code in expr_codes) if code else False,
                    "matched_primary_sample": bool(code in primary_codes) if code else False,
                }
            )
    inventory = pd.DataFrame(rows)
    if inventory.empty:
        note_df = pd.DataFrame({"Note": ["No RNA BAM files were found under the expected project directories."]})
        return inventory, pd.DataFrame(), note_df, pd.DataFrame(), pd.DataFrame()
    if not rna_bed or not Path(rna_bed).exists():
        inventory["RNA_Coverage_QC_Status"] = "Unavailable"
        inventory["RNA_Coverage_Risk_Flag"] = True
        inventory["RNA_Coverage_Risk_Level"] = "Unavailable"
        inventory["RNA_Coverage_Risk_Note"] = "RNA BED was not available, so only BAM inventory could be performed."
        summary = inventory.groupby("RNA_BAM_Group", dropna=False).agg(n_bams=("RNA_BAM_Name", "count"), matched_expression=("matched_expression_sample", "sum"), matched_primary=("matched_primary_sample", "sum")).reset_index()
        summary["validation_note"] = inventory["RNA_Coverage_Risk_Note"].iloc[0]
        return inventory, pd.DataFrame(), summary, pd.DataFrame(), pd.DataFrame()
    panel_design = load_rna_panel_design(rna_bed)
    if panel_design.empty:
        inventory["RNA_Coverage_QC_Status"] = "Unavailable"
        inventory["RNA_Coverage_Risk_Flag"] = True
        inventory["RNA_Coverage_Risk_Level"] = "Unavailable"
        inventory["RNA_Coverage_Risk_Note"] = "RNA BED was present but no targets could be parsed from it."
        summary = pd.DataFrame({"Note": ["RNA BED was present but no targets could be parsed from it."]})
        return inventory, pd.DataFrame(), summary, pd.DataFrame(), panel_design
    samtools_available = shutil.which("samtools") is not None
    target_cov_rows: list[pd.DataFrame] = []
    enriched_rows: list[dict[str, object]] = []
    for _, row in inventory.iterrows():
        out = row.to_dict()
        out.setdefault("RNA_Coverage_QC_Status", "Unavailable")
        out.setdefault("RNA_Coverage_Risk_Flag", True)
        out.setdefault("RNA_Coverage_Risk_Level", "Unavailable")
        out.setdefault("RNA_Coverage_Risk_Note", "RNA coverage QC was not run.")
        if not samtools_available:
            out["RNA_Coverage_Risk_Note"] = "samtools is not available in this runtime, so RNA coverage QC could not be run."
            enriched_rows.append(out)
            continue
        index_ready, index_created, index_note = ensure_bam_index(out["RNA_BAM_Path"])
        out["RNA_BAI_Created"] = bool(index_created)
        out["RNA_BAI_Present"] = bool(index_ready)
        if not index_ready:
            out["RNA_Coverage_Risk_Note"] = f"RNA BAM indexing failed before coverage QC: {index_note}"
            enriched_rows.append(out)
            continue
        try:
            cov_df = run_bedcov(rna_bed, out["RNA_BAM_Path"], panel_design)
            cov_df["target_read_count"] = run_target_read_counts(out["RNA_BAM_Path"], panel_design)
        except Exception as exc:
            out["RNA_Coverage_Risk_Note"] = f"RNA coverage QC failed: {exc}"
            enriched_rows.append(out)
            continue
        cov_df["RNA_BAM_Group"] = out["RNA_BAM_Group"]
        cov_df["RNA_BAM_Name"] = out["RNA_BAM_Name"]
        cov_df["RNA_BAM_Path"] = out["RNA_BAM_Path"]
        cov_df["snp_code"] = out["snp_code"]
        cov_df["matched_expression_sample"] = out["matched_expression_sample"]
        cov_df["matched_primary_sample"] = out["matched_primary_sample"]
        target_cov_rows.append(cov_df)
        out.update(summarise_bam_coverage(cov_df))
        enriched_rows.append(out)
    inventory = pd.DataFrame(enriched_rows)
    target_cov = pd.concat(target_cov_rows, ignore_index=True) if target_cov_rows else pd.DataFrame()
    sample_qc = pd.DataFrame()
    ranked = inventory[inventory["snp_code"].notna()].copy() if not inventory.empty else pd.DataFrame()
    if not ranked.empty:
        rank_map = {"Low": 0, "Moderate": 1, "High": 2, "Unavailable": 3}
        ranked["_risk_rank"] = ranked["RNA_Coverage_Risk_Level"].map(rank_map).fillna(-1)
        ranked["_read_rank"] = pd.to_numeric(ranked.get("RNA_Panel_Min_Target_Reads"), errors="coerce").fillna(-1)
        ranked = ranked.sort_values(["_risk_rank", "_read_rank"], ascending=[True, False])
        sample_qc = ranked.drop_duplicates("snp_code").copy()
        sample_qc = sample_qc[["snp_code"] + [col for col in RNA_QC_COLUMNS if col in sample_qc.columns]]
    summary = inventory.groupby("RNA_BAM_Group", dropna=False).agg(
        n_bams=("RNA_BAM_Name", "count"),
        matched_expression=("matched_expression_sample", "sum"),
        matched_primary=("matched_primary_sample", "sum"),
        indexed_bams=("RNA_BAI_Present", "sum"),
        new_indices=("RNA_BAI_Created", "sum"),
        completed_qc=("RNA_Coverage_QC_Status", lambda series: int((series == "Complete").sum())),
        unavailable_qc=("RNA_Coverage_QC_Status", lambda series: int((series != "Complete").sum())),
        high_risk_bams=("RNA_Coverage_Risk_Level", lambda series: int((series == "High").sum())),
        moderate_risk_bams=("RNA_Coverage_Risk_Level", lambda series: int((series == "Moderate").sum())),
        unavailable_risk_bams=("RNA_Coverage_Risk_Level", lambda series: int((series == "Unavailable").sum())),
        mean_panel_depth=("RNA_Panel_Mean_Depth", "mean"),
        mean_min_target_reads=("RNA_Panel_Min_Target_Reads", "mean"),
        mean_gsdmb_depth=("RNA_GSDMB_Mean_Depth", "mean"),
    ).reset_index()
    summary["rna_bed"] = str(rna_bed)
    summary["validation_note"] = np.where(summary["completed_qc"] > 0, "BED-aware RNA coverage QC completed with transcript-target coverage metrics.", "Only BAM inventory could be completed for this group.")
    return inventory, target_cov, summary, sample_qc, panel_design


def attach_rna_qc(expr: pd.DataFrame, sample_qc: pd.DataFrame) -> pd.DataFrame:
    out = expr.copy()
    if not sample_qc.empty:
        out = out.merge(sample_qc, on="snp_code", how="left")
    for col in RNA_QC_COLUMNS:
        if col not in out.columns:
            out[col] = np.nan
    out["RNA_Coverage_Risk_Flag"] = out["RNA_Coverage_Risk_Flag"].fillna(False).astype(bool)
    mask = out["master_join_success"] & out["snp_code"].notna() & out["RNA_BAM_Path"].isna()
    out.loc[mask, "RNA_Coverage_QC_Status"] = "Missing_BAM"
    out.loc[mask, "RNA_Coverage_Risk_Flag"] = True
    out.loc[mask, "RNA_Coverage_Risk_Level"] = "Missing_BAM"
    out.loc[mask, "RNA_Coverage_Risk_Note"] = "No matched RNA BAM was found for this expression sample under the expected project directories."
    for col in ["RNA_BAM_Group", "RNA_BAM_Path", "RNA_BAM_Name", "RNA_Coverage_QC_Status", "RNA_Coverage_Risk_Level", "RNA_Coverage_Risk_Note", "RNA_Machine_QC_Flag"]:
        out[col] = out[col].fillna("")
    return out


def attach_dna_baseline(expr: pd.DataFrame, dna_baseline: pd.DataFrame) -> pd.DataFrame:
    out = expr.copy()
    if not dna_baseline.empty:
        dedup = (
            dna_baseline.sort_values(["DNA_Baseline_Present", "DNA_BAM_Size_MB"], ascending=[False, False])
            .drop_duplicates("snp_code")
            .copy()
        )
        keep_cols = ["snp_code"] + [col for col in DNA_BASELINE_COLUMNS if col in dedup.columns]
        out = out.merge(dedup[keep_cols], on="snp_code", how="left")
    for col in DNA_BASELINE_COLUMNS:
        if col not in out.columns:
            out[col] = np.nan
    out["DNA_Baseline_Present"] = out["DNA_Baseline_Present"].fillna(False).astype(bool)
    for col in ["DNA_BAM_Group", "DNA_BAM_Path", "DNA_BAM_Name"]:
        out[col] = out[col].fillna("")
    return out


RNA_MANIFEST_DECISION_COLUMNS = [
    "RNA_QC_Analysis_Ready",
    "RNA_QC_Exploratory_Ready",
    "RNA_QC_Final_Status",
    "RNA_QC_Exclusion_Reason",
    "RNA_QC_Eligibility_Note",
]


def classify_rna_analysis_row(row: pd.Series) -> tuple[bool, bool, str, str, str]:
    if not bool(row.get("master_join_success", False)) or pd.isna(row.get("snp_code")):
        return False, False, "Excluded", "Unmatched to harmonised master", "The expression row could not be reconciled to a harmonised RNA sample."

    qc_status = str(row.get("RNA_Coverage_QC_Status", "") or "").strip()
    risk_level = str(row.get("RNA_Coverage_Risk_Level", "") or "").strip()
    risk_note = str(row.get("RNA_Coverage_Risk_Note", "") or "").strip()
    bam_path = str(row.get("RNA_BAM_Path", "") or "").strip()
    bai_present = bool(row.get("RNA_BAI_Present", False))
    dna_baseline_present = bool(row.get("DNA_Baseline_Present", False))
    panel_target_count = pd.to_numeric(pd.Series([row.get("RNA_Panel_Target_Count")]), errors="coerce").iloc[0]
    targets_with_read_counts = pd.to_numeric(pd.Series([row.get("RNA_Panel_Targets_With_Read_Counts")]), errors="coerce").iloc[0]
    panel_min_target_reads = pd.to_numeric(pd.Series([row.get("RNA_Panel_Min_Target_Reads")]), errors="coerce").iloc[0]

    if not dna_baseline_present:
        return False, False, "Excluded", "No DNA baseline BAM", "This sample was not present in the current DNA pass-manifest baseline."
    if not bam_path or qc_status == "Missing_BAM":
        return False, False, "Excluded", "Missing matched RNA BAM", "A DNA baseline BAM was present, but no matched RNA BAM was available for this sample."
    if not bai_present:
        return False, False, "Excluded", "BAM/index failure", risk_note or "RNA BAM indexing was not available for downstream QC."
    if qc_status != "Complete":
        lowered = risk_note.lower()
        if "bed" in lowered or "target" in lowered:
            return False, False, "Excluded", "Missing/invalid RNA BED QC context", risk_note or "RNA target design context was not available."
        if "index" in lowered or "samtools" in lowered:
            return False, False, "Excluded", "BAM/index failure", risk_note or "RNA QC could not be completed because BAM indexing failed."
        return False, False, "Excluded", "RNA QC incomplete", risk_note or "RNA coverage QC did not complete."
    if pd.isna(panel_target_count) or panel_target_count <= 0:
        return False, False, "Excluded", "Missing/invalid RNA BED QC context", risk_note or "The RNA panel design did not yield a measurable panel-gene target set."
    if pd.isna(targets_with_read_counts) or targets_with_read_counts < panel_target_count:
        return False, False, "Excluded", "Missing/invalid panel-gene read-count context", risk_note or "Per-gene read counts were not available across the full RNA panel target set."
    if risk_level == "High":
        return False, True, "Exploratory", "", risk_note or "RNA panel-wide read uniformity was low; retain as exploratory because across-panel non-uniformity may reflect biology or panel bias as well as technical weakness."
    if risk_level == "Moderate":
        return False, True, "Exploratory", "", risk_note or "RNA panel coverage was usable but below the strict analysis-ready uniformity threshold; treat as exploratory."
    if risk_level == "Unavailable":
        return False, False, "Excluded", "Missing/invalid RNA BED QC context", risk_note or "RNA target coverage context was unavailable for this sample."
    return True, True, "Analysis-ready", "", risk_note or "RNA BAM, panel-wide depth, and per-gene read counts passed the strict QC gate."


def build_analysis_manifest(expr_with_qc: pd.DataFrame) -> pd.DataFrame:
    manifest = expr_with_qc.copy()
    decisions = manifest.apply(classify_rna_analysis_row, axis=1, result_type="expand")
    decisions.columns = RNA_MANIFEST_DECISION_COLUMNS
    manifest = pd.concat([manifest, decisions], axis=1)
    manifest["RNA_QC_Analysis_Ready"] = manifest["RNA_QC_Analysis_Ready"].fillna(False).astype(bool)
    manifest["RNA_QC_Exploratory_Ready"] = manifest["RNA_QC_Exploratory_Ready"].fillna(False).astype(bool)
    manifest["RNA_QC_Final_Status"] = manifest["RNA_QC_Final_Status"].fillna("Excluded")
    manifest["RNA_QC_Exclusion_Reason"] = manifest["RNA_QC_Exclusion_Reason"].fillna("")
    manifest["RNA_QC_Eligibility_Note"] = manifest["RNA_QC_Eligibility_Note"].fillna("")
    return manifest


def summarise_analysis_manifest(manifest: pd.DataFrame) -> pd.DataFrame:
    if manifest.empty:
        return pd.DataFrame(columns=["group", "total_rows", "matched_rows", "qc_passed_rows", "qc_exploratory_rows", "qc_exploratory_only_rows", "qc_excluded_rows", "pass_percentage", "exploratory_percentage"])
    if "DNA_Baseline_Present" in manifest.columns:
        manifest = manifest[manifest["DNA_Baseline_Present"].fillna(False)].copy()
    if manifest.empty:
        return pd.DataFrame(columns=["group", "total_rows", "matched_rows", "qc_passed_rows", "qc_exploratory_rows", "qc_exploratory_only_rows", "qc_excluded_rows", "pass_percentage", "exploratory_percentage"])
    rows = []
    groups = [g for g in ["Breast_Tumour", "Endometrial_Tumour", "Breast_Normal", "Endometrial_Normal"] if g in manifest["analysis_group"].astype(str).unique()]
    for group in groups + ["Overall"]:
        sub = manifest if group == "Overall" else manifest[manifest["analysis_group"].astype(str) == group]
        total = len(sub)
        matched = int(sub["master_join_success"].fillna(False).sum()) if total else 0
        passed = int(sub["RNA_QC_Analysis_Ready"].fillna(False).sum()) if total else 0
        exploratory = int(sub["RNA_QC_Exploratory_Ready"].fillna(False).sum()) if total else 0
        exploratory_only = max(exploratory - passed, 0)
        excluded = total - exploratory
        rows.append({
            "group": group,
            "total_rows": total,
            "matched_rows": matched,
            "qc_passed_rows": passed,
            "qc_exploratory_rows": exploratory,
            "qc_exploratory_only_rows": exploratory_only,
            "qc_excluded_rows": excluded,
            "pass_percentage": (passed / total * 100.0) if total else 0.0,
            "exploratory_percentage": (exploratory / total * 100.0) if total else 0.0,
        })
    return pd.DataFrame(rows)


def summarise_exclusion_reasons(manifest: pd.DataFrame) -> pd.DataFrame:
    if manifest.empty:
        return pd.DataFrame(columns=["group", "exclusion_reason", "sample_count", "excluded_total", "excluded_percentage"])
    if "DNA_Baseline_Present" in manifest.columns:
        manifest = manifest[manifest["DNA_Baseline_Present"].fillna(False)].copy()
    if manifest.empty:
        return pd.DataFrame(columns=["group", "exclusion_reason", "sample_count", "excluded_total", "excluded_percentage"])
    rows = []
    groups = [g for g in ["Breast_Tumour", "Endometrial_Tumour", "Breast_Normal", "Endometrial_Normal"] if g in manifest["analysis_group"].astype(str).unique()]
    excluded = manifest[manifest["RNA_QC_Final_Status"].astype(str) == "Excluded"].copy()
    for group in groups + ["Overall"]:
        sub = excluded if group == "Overall" else excluded[excluded["analysis_group"].astype(str) == group]
        total = len(sub)
        counts = sub["RNA_QC_Exclusion_Reason"].replace("", "Unspecified").value_counts() if total else pd.Series(dtype=int)
        if counts.empty:
            rows.append({"group": group, "exclusion_reason": "No excluded samples", "sample_count": 0, "excluded_total": total, "excluded_percentage": 0.0})
        else:
            for reason, count in counts.items():
                rows.append({
                    "group": group,
                    "exclusion_reason": reason,
                    "sample_count": int(count),
                    "excluded_total": total,
                    "excluded_percentage": (count / total * 100.0) if total else 0.0,
                })
    return pd.DataFrame(rows)


def summarise_machine_qc(manifest: pd.DataFrame) -> pd.DataFrame:
    if manifest.empty or "RNA_Machine_QC_Flag" not in manifest.columns:
        return pd.DataFrame(columns=["group", "machine_qc_flag", "sample_count", "median_machine_reads", "median_machine_q20_pct"])
    work = manifest.copy()
    if "DNA_Baseline_Present" in work.columns:
        work = work[work["DNA_Baseline_Present"].fillna(False)].copy()
    work = work[work["RNA_Machine_QC_Flag"].astype(str).ne("")].copy()
    if work.empty:
        return pd.DataFrame(columns=["group", "machine_qc_flag", "sample_count", "median_machine_reads", "median_machine_q20_pct"])
    rows = []
    groups = [g for g in ["Breast_Tumour", "Endometrial_Tumour", "Breast_Normal", "Endometrial_Normal"] if g in work["analysis_group"].astype(str).unique()]
    for group in groups + ["Overall"]:
        sub = work if group == "Overall" else work[work["analysis_group"].astype(str) == group]
        for flag, chunk in sub.groupby(sub["RNA_Machine_QC_Flag"].astype(str), dropna=False):
            rows.append({
                "group": group,
                "machine_qc_flag": flag,
                "sample_count": int(len(chunk)),
                "median_machine_reads": float(pd.to_numeric(chunk.get("RNA_Machine_Reads"), errors="coerce").median()) if len(chunk) else np.nan,
                "median_machine_q20_pct": float(pd.to_numeric(chunk.get("RNA_Machine_Q20_Pct"), errors="coerce").median()) if len(chunk) else np.nan,
            })
    return pd.DataFrame(rows)


def plot_rna_gate_summary(summary_df: pd.DataFrame, out_path: str | Path) -> None:
    if summary_df.empty:
        return
    plot_df = summary_df.copy()
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    xpos = np.arange(len(plot_df))
    exploratory_only = plot_df.get("qc_exploratory_only_rows", pd.Series(0, index=plot_df.index))
    ax.bar(xpos, plot_df["qc_passed_rows"], color="#2e8b57", width=0.58, label="Analysis-ready")
    ax.bar(xpos, exploratory_only, bottom=plot_df["qc_passed_rows"], color="#e9a03b", width=0.58, label="Exploratory")
    ax.bar(xpos, plot_df["qc_excluded_rows"], bottom=plot_df["qc_passed_rows"] + exploratory_only, color="#d95d39", width=0.58, label="Excluded")
    for idx, row in plot_df.iterrows():
        strict_label = f"strict {row['pass_percentage']:.1f}%"
        exploratory_pct = row.get("exploratory_percentage", row["pass_percentage"])
        exp_label = f"usable {exploratory_pct:.1f}%"
        ax.text(idx, row["total_rows"] + max(plot_df["total_rows"].max() * 0.02, 0.2), f"{strict_label}\n{exp_label}", ha="center", va="bottom", fontsize=8.0)
    ax.set_xticks(xpos)
    ax.set_xticklabels([str(x).replace("_", "\n") for x in plot_df["group"]], fontsize=8.5)
    ax.set_ylabel("Rows")
    ax.set_title("Objective 2 RNA QC tier summary", fontweight="bold")
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.04))
    ax.grid(axis="y", linestyle=":", alpha=0.3)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_rna_exclusion_summary(exclusion_df: pd.DataFrame, out_path: str | Path) -> None:
    if exclusion_df.empty:
        return
    plot_df = exclusion_df[exclusion_df["group"] == "Overall"].copy()
    fig, ax = plt.subplots(figsize=(8.5, 5.2), constrained_layout=True)
    if plot_df.empty or plot_df["sample_count"].sum() == 0:
        ax.text(0.5, 0.5, "No excluded RNA samples", ha="center", va="center", fontsize=11, fontweight="bold")
        ax.axis("off")
    else:
        wedges, _ = ax.pie(plot_df["sample_count"], startangle=90, wedgeprops={"width": 0.42, "edgecolor": "white"})
        labels = [f"{row.exclusion_reason} ({int(row.sample_count)})" for row in plot_df.itertuples()]
        ax.legend(wedges, labels, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
        ax.set_title("Objective 2 RNA exclusion reasons", fontweight="bold")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_rna_qc(sample_qc: pd.DataFrame, out_dir: str | Path) -> None:
    if sample_qc.empty:
        return
    plot_df = sample_qc[sample_qc["RNA_Coverage_QC_Status"].astype(str) == "Complete"].copy()
    if plot_df.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    order = [group for group in ["Breast_Tumour", "Endometrial_Tumour", "Breast_Normal", "Endometrial_Normal"] if group in plot_df["RNA_BAM_Group"].unique()]
    metrics = [
        ("RNA_Panel_Targets_GE20_Pct", f"Panel targets >={RNA_PANEL_LOW_DEPTH:.0f}x"),
        ("RNA_Panel_Min_Target_Reads", "Minimum panel-gene read count"),
    ]
    for ax, (metric, title) in zip(axes.flat, metrics):
        sub = plot_df[["RNA_BAM_Group", metric]].dropna().copy()
        if sub.empty:
            ax.set_visible(False)
            continue
        if sns is not None:
            sns.boxplot(data=sub, x="RNA_BAM_Group", y=metric, order=order, ax=ax, color="#f1e0cc")
            sns.stripplot(data=sub, x="RNA_BAM_Group", y=metric, order=order, ax=ax, color="#7f4f24", size=3, alpha=0.7)
        else:
            labels: list[str] = []
            groups: list[np.ndarray] = []
            for label in order:
                values = pd.to_numeric(sub.loc[sub["RNA_BAM_Group"].astype(str) == label, metric], errors="coerce").dropna().values
                if len(values):
                    labels.append(label)
                    groups.append(values)
            if groups:
                ax.boxplot(groups, labels=labels, patch_artist=True)
                for i, values in enumerate(groups, start=1):
                    ax.scatter(np.full(len(values), i), values, color="#7f4f24", s=8, alpha=0.7)
        ax.set_title(title)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=25)
    fig.suptitle("Objective 2: RNA panel depth and per-gene read-count QC", fontsize=14, fontweight="bold")
    fig.savefig(Path(out_dir) / "20_RNA_Coverage_QC.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_rna_machine_vs_panel_qc(manifest: pd.DataFrame, out_path: str | Path) -> None:
    if manifest.empty or "RNA_Machine_Reads" not in manifest.columns:
        return
    plot_df = manifest.copy()
    if "DNA_Baseline_Present" in plot_df.columns:
        plot_df = plot_df[plot_df["DNA_Baseline_Present"].fillna(False)].copy()
    plot_df = plot_df[plot_df["RNA_Machine_Reads"].notna() & plot_df["RNA_Machine_Q20_Pct"].notna()].copy()
    if plot_df.empty:
        return
    status_colors = {
        "Analysis-ready": "#2e8b57",
        "Exploratory": "#e9a03b",
        "Excluded": "#d95d39",
    }
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), constrained_layout=True)
    for status in ["Analysis-ready", "Exploratory", "Excluded"]:
        sub = plot_df[plot_df["RNA_QC_Final_Status"].astype(str) == status].copy()
        if sub.empty:
            continue
        axes[0].scatter(
            pd.to_numeric(sub["RNA_Machine_Reads"], errors="coerce"),
            pd.to_numeric(sub["RNA_Panel_Targets_GE100_Reads_Pct"], errors="coerce") * 100.0,
            s=30, alpha=0.8, label=status, color=status_colors[status], edgecolor="white", linewidth=0.3,
        )
        axes[1].scatter(
            pd.to_numeric(sub["RNA_Machine_Q20_Pct"], errors="coerce"),
            pd.to_numeric(sub["RNA_Panel_Targets_GE20_Pct"], errors="coerce") * 100.0,
            s=30, alpha=0.8, label=status, color=status_colors[status], edgecolor="white", linewidth=0.3,
        )
    axes[0].axvline(RNA_MACHINE_HARD_FAIL_READS, color="#c0392b", linestyle="--", linewidth=1.1)
    axes[0].axvline(RNA_MACHINE_WARN_READS, color="#d68910", linestyle=":", linewidth=1.2)
    axes[0].axhline(RNA_GENE_STRICT_TARGET_PCT * 100.0, color="#2e86c1", linestyle="--", linewidth=1.1)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Machine reads (log scale)")
    axes[0].set_ylabel(f"Panel targets >={RNA_GENE_STRICT_MIN_READS} reads (%)")
    axes[0].set_title("Reads vs panel-gene uniformity")
    axes[0].grid(True, linestyle=":", alpha=0.3)

    axes[1].axvline(RNA_MACHINE_HARD_FAIL_Q20_PCT, color="#c0392b", linestyle="--", linewidth=1.1)
    axes[1].axvline(RNA_MACHINE_WARN_Q20_PCT, color="#d68910", linestyle=":", linewidth=1.2)
    axes[1].axhline(RNA_PANEL_STRICT_TARGET_PCT * 100.0, color="#2e86c1", linestyle="--", linewidth=1.1)
    axes[1].set_xlabel("Machine Q20 bases (%)")
    axes[1].set_ylabel(f"Panel targets >={RNA_PANEL_LOW_DEPTH:.0f}x (%)")
    axes[1].set_title("Q20 vs panel depth coverage")
    axes[1].grid(True, linestyle=":", alpha=0.3)

    handles, labels = axes[1].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, frameon=False, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.03))
    matched_n = int(plot_df["snp_code"].astype(str).nunique()) if "snp_code" in plot_df.columns else len(plot_df)
    fig.suptitle(f"RNA machine QC vs panel QC (matched baseline subset: n={matched_n})", fontsize=14, fontweight="bold")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
