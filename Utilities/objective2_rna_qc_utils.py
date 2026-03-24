from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    import seaborn as sns
except ImportError:
    sns = None

RNA_PANEL_LOW_DEPTH = 20.0
RNA_GSDMB_MODERATE_DEPTH = 100.0
RNA_PANEL_COVERAGE_PCT = 0.90
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
    "RNA_Panel_Low_Depth_Targets",
    "RNA_Panel_Mean_Depth",
    "RNA_Panel_Median_Depth",
    "RNA_GSDMB_Target_Count",
    "RNA_GSDMB_Targets_Covered",
    "RNA_GSDMB_Mean_Depth",
    "RNA_GSDMB_Min_Depth",
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
        parts = line.split("\t")
        coverage_values.append(float(parts[-1]))
    if len(coverage_values) != len(panel_design):
        raise RuntimeError(f"bedcov returned {len(coverage_values)} rows for {len(panel_design)} design targets")
    out = panel_design.copy()
    out["covered_base_sum"] = coverage_values
    out["mean_depth"] = np.where(out["target_length"] > 0, out["covered_base_sum"] / out["target_length"], np.nan)
    return out

def classify_rna_coverage(panel_pct: float, panel_mean_depth: float, gsdmb_target_count: int, gsdmb_targets_covered: int, gsdmb_mean_depth: float) -> tuple[bool, str, str]:
    if gsdmb_target_count == 0:
        return True, "Unavailable", "The RNA panel design did not contain a GSDMB target row to validate."
    if gsdmb_targets_covered == 0 or pd.isna(gsdmb_mean_depth) or gsdmb_mean_depth <= 0:
        return True, "High", "No measurable coverage was detected across the GSDMB RNA panel target; interpret isoform quantification very cautiously."
    if gsdmb_mean_depth < RNA_PANEL_LOW_DEPTH:
        return True, "High", f"GSDMB RNA panel coverage was very low (mean depth {gsdmb_mean_depth:.1f}); interpret isoform quantification very cautiously."
    if panel_pct < RNA_PANEL_COVERAGE_PCT or panel_mean_depth < RNA_PANEL_LOW_DEPTH or gsdmb_mean_depth < RNA_GSDMB_MODERATE_DEPTH:
        return True, "Moderate", f"RNA panel coverage was usable but suboptimal (panel target coverage {panel_pct:.1%}, GSDMB mean depth {gsdmb_mean_depth:.1f}); interpret with caution."
    return False, "None", "RNA panel and GSDMB target coverage were adequate for routine interpretation."


def summarise_bam_coverage(target_cov: pd.DataFrame) -> dict[str, object]:
    panel_target_count = len(target_cov)
    panel_targets_covered = int((target_cov["mean_depth"] > 0).sum())
    panel_pct = panel_targets_covered / panel_target_count if panel_target_count else np.nan
    panel_mean_depth = pd.to_numeric(target_cov["mean_depth"], errors="coerce").mean()
    panel_median_depth = pd.to_numeric(target_cov["mean_depth"], errors="coerce").median()
    panel_low_depth_targets = int((pd.to_numeric(target_cov["mean_depth"], errors="coerce") < RNA_PANEL_LOW_DEPTH).sum())
    gsdmb_cov = target_cov[target_cov["gene_symbol"].astype(str).str.upper() == "GSDMB"].copy()
    gsdmb_target_count = len(gsdmb_cov)
    gsdmb_targets_covered = int((pd.to_numeric(gsdmb_cov.get("mean_depth"), errors="coerce") > 0).sum()) if gsdmb_target_count else 0
    gsdmb_mean_depth = pd.to_numeric(gsdmb_cov.get("mean_depth"), errors="coerce").mean() if gsdmb_target_count else np.nan
    gsdmb_min_depth = pd.to_numeric(gsdmb_cov.get("mean_depth"), errors="coerce").min() if gsdmb_target_count else np.nan
    risk_flag, risk_level, risk_note = classify_rna_coverage(panel_pct, panel_mean_depth, gsdmb_target_count, gsdmb_targets_covered, gsdmb_mean_depth)
    return {
        "RNA_Coverage_QC_Status": "Complete",
        "RNA_Panel_Target_Count": panel_target_count,
        "RNA_Panel_Targets_Covered": panel_targets_covered,
        "RNA_Panel_Target_Coverage_Pct": panel_pct,
        "RNA_Panel_Low_Depth_Targets": panel_low_depth_targets,
        "RNA_Panel_Mean_Depth": panel_mean_depth,
        "RNA_Panel_Median_Depth": panel_median_depth,
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
        rank_map = {"None": 0, "Moderate": 1, "High": 2, "Unavailable": 3}
        ranked["_risk_rank"] = ranked["RNA_Coverage_Risk_Level"].map(rank_map).fillna(-1)
        ranked["_depth_rank"] = pd.to_numeric(ranked.get("RNA_GSDMB_Mean_Depth"), errors="coerce").fillna(-1)
        ranked = ranked.sort_values(["_risk_rank", "_depth_rank"], ascending=[False, False])
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
    for col in ["RNA_BAM_Group", "RNA_BAM_Path", "RNA_BAM_Name", "RNA_Coverage_QC_Status", "RNA_Coverage_Risk_Level", "RNA_Coverage_Risk_Note"]:
        out[col] = out[col].fillna("")
    return out


RNA_MANIFEST_DECISION_COLUMNS = [
    "RNA_QC_Analysis_Ready",
    "RNA_QC_Final_Status",
    "RNA_QC_Exclusion_Reason",
    "RNA_QC_Eligibility_Note",
]


def classify_rna_analysis_row(row: pd.Series) -> tuple[bool, str, str, str]:
    if not bool(row.get("master_join_success", False)) or pd.isna(row.get("snp_code")):
        return False, "Excluded", "Unmatched to harmonised master", "The expression row could not be reconciled to a harmonised RNA sample." 

    qc_status = str(row.get("RNA_Coverage_QC_Status", "") or "").strip()
    risk_level = str(row.get("RNA_Coverage_Risk_Level", "") or "").strip()
    risk_note = str(row.get("RNA_Coverage_Risk_Note", "") or "").strip()
    bam_path = str(row.get("RNA_BAM_Path", "") or "").strip()
    bai_present = bool(row.get("RNA_BAI_Present", False))
    gsdmb_target_count = pd.to_numeric(pd.Series([row.get("RNA_GSDMB_Target_Count")]), errors="coerce").iloc[0]
    gsdmb_targets_covered = pd.to_numeric(pd.Series([row.get("RNA_GSDMB_Targets_Covered")]), errors="coerce").iloc[0]
    gsdmb_mean_depth = pd.to_numeric(pd.Series([row.get("RNA_GSDMB_Mean_Depth")]), errors="coerce").iloc[0]

    if not bam_path or qc_status == "Missing_BAM":
        return False, "Excluded", "Missing BAM", "No matched RNA BAM was available for this sample."
    if not bai_present:
        return False, "Excluded", "BAM/index failure", risk_note or "RNA BAM indexing was not available for downstream QC."
    if qc_status != "Complete":
        lowered = risk_note.lower()
        if "bed" in lowered or "target" in lowered:
            return False, "Excluded", "Missing/invalid RNA BED QC context", risk_note or "RNA target design context was not available."
        if "index" in lowered or "samtools" in lowered:
            return False, "Excluded", "BAM/index failure", risk_note or "RNA QC could not be completed because BAM indexing failed."
        return False, "Excluded", "RNA QC incomplete", risk_note or "RNA coverage QC did not complete."
    if pd.isna(gsdmb_target_count) or gsdmb_target_count <= 0:
        return False, "Excluded", "Missing/invalid RNA BED QC context", risk_note or "The RNA panel design did not yield a measurable GSDMB target set."
    if pd.isna(gsdmb_targets_covered) or gsdmb_targets_covered <= 0 or pd.isna(gsdmb_mean_depth) or gsdmb_mean_depth <= 0:
        return False, "Excluded", "No measurable GSDMB target coverage", risk_note or "No measurable coverage was detected across the GSDMB RNA targets."
    if risk_level == "High":
        return False, "Excluded", "Severe low-depth/high-risk coverage failure", risk_note or "RNA target coverage was too weak for QC-gated analysis."
    if risk_level in {"Moderate", "Unavailable"}:
        return False, "Excluded", "Suboptimal RNA target coverage", risk_note or "RNA target coverage was below the analysis-ready threshold."
    return True, "Analysis-ready", "", risk_note or "RNA BAM and GSDMB target coverage passed the QC gate."


def build_analysis_manifest(expr_with_qc: pd.DataFrame) -> pd.DataFrame:
    manifest = expr_with_qc.copy()
    decisions = manifest.apply(classify_rna_analysis_row, axis=1, result_type="expand")
    decisions.columns = RNA_MANIFEST_DECISION_COLUMNS
    manifest = pd.concat([manifest, decisions], axis=1)
    manifest["RNA_QC_Analysis_Ready"] = manifest["RNA_QC_Analysis_Ready"].fillna(False).astype(bool)
    manifest["RNA_QC_Final_Status"] = manifest["RNA_QC_Final_Status"].fillna("Excluded")
    manifest["RNA_QC_Exclusion_Reason"] = manifest["RNA_QC_Exclusion_Reason"].fillna("")
    manifest["RNA_QC_Eligibility_Note"] = manifest["RNA_QC_Eligibility_Note"].fillna("")
    return manifest


def summarise_analysis_manifest(manifest: pd.DataFrame) -> pd.DataFrame:
    if manifest.empty:
        return pd.DataFrame(columns=["group", "total_rows", "matched_rows", "qc_passed_rows", "qc_excluded_rows", "pass_percentage"])
    rows = []
    groups = [g for g in ["Breast_Tumour", "Endometrial_Tumour", "Breast_Normal", "Endometrial_Normal"] if g in manifest["analysis_group"].astype(str).unique()]
    for group in groups + ["Overall"]:
        sub = manifest if group == "Overall" else manifest[manifest["analysis_group"].astype(str) == group]
        total = len(sub)
        matched = int(sub["master_join_success"].fillna(False).sum()) if total else 0
        passed = int(sub["RNA_QC_Analysis_Ready"].fillna(False).sum()) if total else 0
        excluded = total - passed
        rows.append({
            "group": group,
            "total_rows": total,
            "matched_rows": matched,
            "qc_passed_rows": passed,
            "qc_excluded_rows": excluded,
            "pass_percentage": (passed / total * 100.0) if total else 0.0,
        })
    return pd.DataFrame(rows)


def summarise_exclusion_reasons(manifest: pd.DataFrame) -> pd.DataFrame:
    if manifest.empty:
        return pd.DataFrame(columns=["group", "exclusion_reason", "sample_count", "excluded_total", "excluded_percentage"])
    rows = []
    groups = [g for g in ["Breast_Tumour", "Endometrial_Tumour", "Breast_Normal", "Endometrial_Normal"] if g in manifest["analysis_group"].astype(str).unique()]
    excluded = manifest[~manifest["RNA_QC_Analysis_Ready"].fillna(False)].copy()
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


def plot_rna_gate_summary(summary_df: pd.DataFrame, out_path: str | Path) -> None:
    if summary_df.empty:
        return
    plot_df = summary_df.copy()
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    xpos = np.arange(len(plot_df))
    ax.bar(xpos, plot_df["total_rows"], color="#ececec", edgecolor="#b0b0b0", width=0.7, label="Total")
    ax.bar(xpos, plot_df["qc_passed_rows"], color="#2e8b57", width=0.48, label="QC-passed")
    ax.bar(xpos, plot_df["qc_excluded_rows"], bottom=plot_df["qc_passed_rows"], color="#d95d39", width=0.48, label="Excluded")
    for idx, row in plot_df.iterrows():
        ax.text(idx, row["total_rows"] + max(plot_df["total_rows"].max() * 0.02, 0.2), f"{row['pass_percentage']:.1f}%", ha="center", va="bottom", fontsize=8.5)
    ax.set_xticks(xpos)
    ax.set_xticklabels([str(x).replace("_", "\n") for x in plot_df["group"]], fontsize=8.5)
    ax.set_ylabel("Rows")
    ax.set_title("Objective 2 RNA QC pass summary", fontweight="bold")
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
        ("RNA_Panel_Target_Coverage_Pct", "Panel target coverage proportion"),
        ("RNA_GSDMB_Mean_Depth", "GSDMB target mean depth"),
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
    fig.suptitle("Objective 2: RNA panel coverage QC", fontsize=14, fontweight="bold")
    fig.savefig(Path(out_dir) / "20_RNA_Coverage_QC.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
