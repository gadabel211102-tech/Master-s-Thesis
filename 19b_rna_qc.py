#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import pandas as pd

from association_runtime import script19b_defaults
from objective2_rna_qc_utils import (
    DNA_BASELINE_COLUMNS,
    RNA_QC_COLUMNS,
    attach_dna_baseline,
    attach_rna_qc,
    bam_validation,
    build_analysis_manifest,
    load_rna_panel_design,
    plot_rna_exclusion_summary,
    plot_rna_gate_summary,
    plot_rna_machine_vs_panel_qc,
    plot_rna_qc,
    summarise_analysis_manifest,
    summarise_bam_coverage,
    summarise_exclusion_reasons,
    summarise_machine_qc,
)

D = script19b_defaults()
TMAP = {
    "breast tumor": ("Breast", "Tumour", "Breast_Tumour", "primary_tumour"),
    "breast tumour": ("Breast", "Tumour", "Breast_Tumour", "primary_tumour"),
    "breast normal": ("Breast", "Healthy", "Breast_Normal", "context_control"),
    "endometrial cancer": ("Endometrial", "Tumour", "Endometrial_Tumour", "primary_tumour"),
    "endometrial normal": ("Endometrial", "Healthy", "Endometrial_Normal", "context_control"),
    "endometrial_normal": ("Endometrial", "Healthy", "Endometrial_Normal", "context_control"),
    "ovarian": ("Ovarian", "Tumour", "Ovarian", "excluded"),
    "ovarian organoids": ("Ovarian", "Tumour", "Ovarian_Organoids", "excluded"),
    "cell line": ("CellLine", "Unknown", "Cell_Line", "excluded"),
}


def s(v):
    return "" if pd.isna(v) else str(v).strip()


def nk(v):
    return re.sub(r"[^A-Z0-9]+", "", s(v).upper())


def first(values):
    for value in values:
        if pd.notna(value) and str(value).strip() not in {"", "nan", "None"}:
            return value
    return pd.NA


def canonical_snp_code(v):
    text = s(v).upper()
    if not text:
        return None
    text = text.replace("SNP_MT-N_", "SNP_MN_")
    text = re.sub(r"_REP$", "", text)
    return text


def sx(v):
    text = s(v)
    patterns = [
        (r"^(SNP_(?:AT|EN|MN|MT-T|MT-N)_\d+)", lambda m: m.group(1).upper()),
        (r"^SNP_DNA_(EN|MN|AT)_(\d+)", lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^SNP_DNA_MT[-_]T_(\d+)", lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^SNP_DNA_MT[-_]N_(\d+)", lambda m: f"SNP_MN_{m.group(1)}"),
        (r"^DNA_SNP_(EN|MN|AT)_(\d+)", lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^DNA_(AT|EN|MN)_(\d+)", lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^DNA_(?:SNP_)?MT[-_]T_(\d+)", lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^MAMAH2_MT[-_]T_(\d+)", lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^DNA_(?:SNP_)?MT[-_]N_(\d+)", lambda m: f"SNP_MN_{m.group(1)}"),
        (r"^MAMAH2_MT[-_]N_(\d+)", lambda m: f"SNP_MN_{m.group(1)}"),
        (r"^RNA_SNP_(AT|EN|MN|MT-T|MT-N)[-_]?(\d+)", lambda m: f"SNP_{'MN' if m.group(1).upper() == 'MT-N' else m.group(1).upper()}_{m.group(2)}"),
        (r"^SNP_RNA_(AT|EN|MN|MT-T|MT-N)_(\d+)", lambda m: f"SNP_{'MN' if m.group(1).upper() == 'MT-N' else m.group(1).upper()}_{m.group(2)}"),
        (r"^SNP_(AT|EN|MN|MT-T|MT-N)_RNA_(\d+)", lambda m: f"SNP_{'MN' if m.group(1).upper() == 'MT-N' else m.group(1).upper()}_{m.group(2)}"),
        (r"^RNA_(AT|EN|MN)_(\d+)", lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^RNA_MT[-_]T_(\d+)", lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^RNA_MT[-_]N_(\d+)", lambda m: f"SNP_MN_{m.group(1)}"),
    ]
    for pat, fmt in patterns:
        match = re.match(pat, text, re.I)
        if match:
            return canonical_snp_code(fmt(match))
    return None


def aliases(v):
    text = s(v)
    out = {nk(text)} if text else set()
    code = sx(text)
    if code:
        out.add(nk(code))
    match = re.match(r"^ENDOMDA[_ -]?(\d+)$", nk(text))
    if match:
        out.add(nk(f"MDA{int(match.group(1)):03d}"))
    for pat in [r"(?:ENDO)?(BL|MDA)[_ -]?0*(\d+)$", r"(?:^|_)(BL|MDA)[_ -]?0*(\d+)"]:
        match = re.search(pat, text, re.I) or re.search(pat, nk(text), re.I)
        if match:
            out.add(nk(f"{match.group(1).upper()}{int(match.group(2))}"))
            out.add(nk(f"{match.group(1).upper()}{int(match.group(2)):03d}"))
    return {x for x in out if x}


def load_expr(path: Path) -> pd.DataFrame:
    xls = pd.ExcelFile(path)
    if "Sheet3" not in xls.sheet_names or "eva" not in xls.sheet_names:
        raise ValueError("Workbook needs Sheet3 and eva")
    df = pd.read_excel(path, sheet_name="Sheet3")
    eva = pd.read_excel(path, sheet_name="eva")
    meta = df["Tissue"].astype(str).str.strip().str.lower().map(TMAP)
    df["cohort"] = meta.map(lambda x: x[0] if isinstance(x, tuple) else "Unknown")
    df["tumour_normal"] = meta.map(lambda x: x[1] if isinstance(x, tuple) else "Unknown")
    df["analysis_group"] = meta.map(lambda x: x[2] if isinstance(x, tuple) else "Unknown")
    df["analysis_role"] = meta.map(lambda x: x[3] if isinstance(x, tuple) else "excluded")
    key = [c for c in ["NOMBRE DE LA MUESTRA", "CODIGO JC"] if c in eva.columns]
    sup = [c for c in ["observaciones EVA", "TIENEN RNA", "FALTA MUESTRA"] if c in eva.columns]
    if key:
        df = df.merge(eva[key + sup].drop_duplicates(), on=key, how="left")
    df["raw1"] = df["CODIGO JC"].fillna(df["NOMBRE DE LA MUESTRA"])
    df["raw2"] = df["NOMBRE DE LA MUESTRA"].fillna(df["CODIGO JC"])
    df["sx1"] = df["raw1"].map(sx)
    df["sx2"] = df["raw2"].map(sx)
    return df


def load_master(path: Path) -> pd.DataFrame:
    master = pd.read_excel(path, sheet_name="harmonised_plus_canon")
    master["snp_code"] = master["snp_code"].astype(str).str.strip()
    master["nucleic_acid"] = master["nucleic_acid"].astype(str).str.upper()
    return master


def parse_machine_integer(value):
    text = s(value)
    digits = re.sub(r"[^0-9]", "", text)
    return int(digits) if digits else pd.NA


def machine_sample_to_snp(value):
    text = s(value)
    if not text:
        return None
    code = sx(text)
    if code:
        return canonical_snp_code(code)
    compact = re.sub(r"\s+", "", text.upper())
    patterns = [
        (r"^SNP_(AT|EN|MN|MT-T|MT-N)_RNA_(\d+)", lambda m: f"SNP_{'MN' if m.group(1) == 'MT-N' else m.group(1)}_{m.group(2)}"),
        (r"^RNA_SNP_(AT|EN|MN|MT-T|MT-N)[-_]?(\d+)", lambda m: f"SNP_{'MN' if m.group(1) == 'MT-N' else m.group(1)}_{m.group(2)}"),
        (r"^SNP_RNA_(AT|EN|MN|MT-T|MT-N)_(\d+)", lambda m: f"SNP_{'MN' if m.group(1) == 'MT-N' else m.group(1)}_{m.group(2)}"),
        (r"^RNA_(AT|EN|MN)_(\d+)", lambda m: f"SNP_{m.group(1)}_{m.group(2)}"),
        (r"^RNA_MT[-_]T_(\d+)", lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^RNA_MT[-_]N_(\d+)", lambda m: f"SNP_MN_{m.group(1)}"),
    ]
    for pattern, formatter in patterns:
        match = re.match(pattern, compact)
        if match:
            return canonical_snp_code(formatter(match))
    return None


def classify_machine_qc_flag(reads, q20_pct):
    reads_val = pd.to_numeric(pd.Series([reads]), errors="coerce").iloc[0]
    q20_val = pd.to_numeric(pd.Series([q20_pct]), errors="coerce").iloc[0]
    if pd.isna(reads_val) and pd.isna(q20_val):
        return "Unmatched"
    if (pd.notna(reads_val) and reads_val < D.get("rna_machine_hard_fail_reads", 10000)) or (pd.notna(q20_val) and q20_val < D.get("rna_machine_hard_fail_q20_pct", 85.0)):
        return "Hard_fail"
    if (pd.notna(reads_val) and reads_val < D.get("rna_machine_warn_reads", 50000)) or (pd.notna(q20_val) and q20_val < D.get("rna_machine_warn_q20_pct", 92.0)):
        return "Warning"
    return "Pass"


def load_machine_qc(path: Path | None) -> pd.DataFrame:
    columns = [
        "machine_match_code",
        "RNA_Machine_Reads",
        "RNA_Machine_Q20_Pct",
        "RNA_Machine_QC_Flag",
        "RNA_Machine_Source_Sample",
        "RNA_Machine_Barcode",
    ]
    if path is None or not path.exists():
        return pd.DataFrame(columns=columns)
    xls = pd.ExcelFile(path)
    sheet_name = "Hoja1" if "Hoja1" in xls.sheet_names else xls.sheet_names[0]
    df = pd.read_excel(path, sheet_name=sheet_name)
    if df.empty:
        return pd.DataFrame(columns=columns)
    df = df.copy()
    df["RNA_Machine_Source_Sample"] = df.get("SAMPLE", pd.Series(index=df.index, dtype=object)).map(s)
    df["RNA_Machine_Barcode"] = df.get("BARCODE NAME", pd.Series(index=df.index, dtype=object)).map(s)
    df["machine_match_code"] = df["RNA_Machine_Source_Sample"].map(machine_sample_to_snp).map(canonical_snp_code)
    df["_bases"] = df.get("BASES", pd.Series(index=df.index, dtype=object)).map(parse_machine_integer)
    df["_q20_bases"] = df.get(">=Q20 BASES", pd.Series(index=df.index, dtype=object)).map(parse_machine_integer)
    df["RNA_Machine_Reads"] = df.get("READS", pd.Series(index=df.index, dtype=object)).map(parse_machine_integer)
    df["RNA_Machine_Q20_Pct"] = pd.to_numeric(df["_q20_bases"], errors="coerce") / pd.to_numeric(df["_bases"], errors="coerce") * 100.0
    df["RNA_Machine_QC_Flag"] = [classify_machine_qc_flag(r, q) for r, q in zip(df["RNA_Machine_Reads"], df["RNA_Machine_Q20_Pct"])]
    df = df[df["machine_match_code"].notna()].copy()
    if df.empty:
        return pd.DataFrame(columns=columns)
    df = df.sort_values(["RNA_Machine_Reads", "RNA_Machine_Q20_Pct"], ascending=[False, False]).drop_duplicates("machine_match_code")
    return df[columns].reset_index(drop=True)


def attach_machine_qc(df: pd.DataFrame, machine_qc: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "snp_code" not in out.columns:
        out["snp_code"] = pd.NA
    stale_cols = [
        "RNA_Machine_Reads",
        "RNA_Machine_Q20_Pct",
        "RNA_Machine_QC_Flag",
        "RNA_Machine_Source_Sample",
        "RNA_Machine_Barcode",
    ]
    out = out.drop(columns=[col for col in stale_cols if col in out.columns], errors="ignore")
    out["_machine_match_code"] = out["snp_code"].map(canonical_snp_code)
    if not machine_qc.empty:
        keep = [
            "machine_match_code",
            "RNA_Machine_Reads",
            "RNA_Machine_Q20_Pct",
            "RNA_Machine_QC_Flag",
            "RNA_Machine_Source_Sample",
            "RNA_Machine_Barcode",
        ]
        out = out.merge(machine_qc[keep], left_on="_machine_match_code", right_on="machine_match_code", how="left")
        out = out.drop(columns=["machine_match_code"], errors="ignore")
    if "RNA_Machine_Reads" not in out.columns:
        out["RNA_Machine_Reads"] = pd.NA
    if "RNA_Machine_Q20_Pct" not in out.columns:
        out["RNA_Machine_Q20_Pct"] = pd.NA
    if "RNA_Machine_QC_Flag" not in out.columns:
        out["RNA_Machine_QC_Flag"] = "Unmatched"
    else:
        out["RNA_Machine_QC_Flag"] = out["RNA_Machine_QC_Flag"].fillna("Unmatched")
    return out.drop(columns=["_machine_match_code"], errors="ignore")


def match_expr(df: pd.DataFrame, master: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    mr = (master[master["nucleic_acid"] == "RNA"]
          .sort_values(["snp_code", "sample_id", "case_id"])
          .groupby("snp_code", dropna=False)
          .agg(first)
          .reset_index(drop=False))
    alias_map = {}
    for _, row in mr.iterrows():
        for alias in aliases(row.get("snp_code")) | aliases(row.get("sample_id")) | aliases(row.get("case_id")):
            alias_map.setdefault(alias, set()).add(row["snp_code"])
    rows = []
    audit = []
    codes = set(mr["snp_code"])
    for _, row in df.iterrows():
        candidates = {x for x in [row.get("sx1"), row.get("sx2")] if x in codes}
        method = "snp_code_regex" if len(candidates) == 1 else ""
        if not candidates:
            alias_candidates = aliases(row.get("raw1")) | aliases(row.get("raw2"))
            candidates = {code for alias in alias_candidates for code in alias_map.get(alias, set())}
            method = "sample_alias" if len(candidates) == 1 else ("ambiguous_alias" if len(candidates) > 1 else "unmatched")
        elif len(candidates) > 1:
            method = "ambiguous_regex"
        out = row.to_dict()
        out["snp_code"] = next(iter(candidates)) if len(candidates) == 1 else pd.NA
        out["match_status"] = "matched" if len(candidates) == 1 else ("ambiguous" if len(candidates) > 1 else "unmatched")
        out["match_method"] = method
        out["match_candidates"] = "; ".join(sorted(candidates)) if len(candidates) > 1 else ""
        rows.append(out)
        if out["match_status"] != "matched":
            audit.append({
                "NOMBRE DE LA MUESTRA": row.get("NOMBRE DE LA MUESTRA"),
                "CODIGO JC": row.get("CODIGO JC"),
                "Tissue": row.get("Tissue"),
                "analysis_group": row.get("analysis_group"),
                "match_status": out["match_status"],
                "match_method": method,
                "match_candidates": out["match_candidates"],
            })
    merged = pd.DataFrame(rows).merge(mr, on="snp_code", how="left", suffixes=("", "_master"))
    merged["master_join_success"] = merged["sample_id"].notna()
    merged["analysis_include_primary"] = (merged["analysis_role"] == "primary_tumour") & merged["master_join_success"]
    return merged, pd.DataFrame(audit)


def rebuild_cached_bam_validation(out_dir: Path, rna_bed: str | Path | None, expr: pd.DataFrame, sx_func):
    inventory_path = out_dir / "RNA_BAM_Inventory.tsv"
    target_cov_path = out_dir / "RNA_Target_Coverage.tsv.gz"
    if not inventory_path.exists() or not target_cov_path.exists() or not rna_bed or not Path(rna_bed).exists():
        return None
    inventory = pd.read_csv(inventory_path, sep="	")
    target_cov = pd.read_csv(target_cov_path, sep="	")
    panel_design = load_rna_panel_design(rna_bed)
    if inventory.empty:
        return None

    # Refresh parser-dependent identity fields from the current code so cache
    # reuse does not preserve stale snp_code assignments after naming fixes.
    expr_codes = set(expr["snp_code"].dropna())
    primary_codes = set(expr.loc[expr["analysis_include_primary"], "snp_code"].dropna())
    if "RNA_BAM_Name" in inventory.columns:
        bam_stems = inventory["RNA_BAM_Name"].astype(str).map(lambda x: Path(x).stem)
    else:
        bam_stems = inventory["RNA_BAM_Path"].astype(str).map(lambda x: Path(x).stem)
    inventory["snp_code"] = bam_stems.map(sx_func)
    inventory["matched_expression_sample"] = inventory["snp_code"].isin(expr_codes)
    inventory["matched_primary_sample"] = inventory["snp_code"].isin(primary_codes)

    metrics_rows = []
    if not target_cov.empty and "RNA_BAM_Path" in target_cov.columns:
        for bam_path, cov_df in target_cov.groupby("RNA_BAM_Path", dropna=False):
            metrics_rows.append({"RNA_BAM_Path": bam_path, **summarise_bam_coverage(cov_df)})
    metrics_df = pd.DataFrame(metrics_rows)
    if not metrics_df.empty:
        metric_cols = [col for col in metrics_df.columns if col != "RNA_BAM_Path"]
        base = inventory.copy()
        for col in metric_cols:
            if col in base.columns:
                base = base.drop(columns=[col])
        inventory = base.merge(metrics_df, on="RNA_BAM_Path", how="left")
    sample_qc = pd.DataFrame()
    ranked = inventory[inventory["snp_code"].notna()].copy() if not inventory.empty else pd.DataFrame()
    if not ranked.empty:
        rank_map = {"Low": 0, "Moderate": 1, "High": 2, "Unavailable": 3}
        ranked["_risk_rank"] = ranked["RNA_Coverage_Risk_Level"].map(rank_map).fillna(-1)
        ranked["_read_rank"] = pd.to_numeric(ranked.get("RNA_Panel_Min_Target_Reads"), errors="coerce").fillna(-1)
        ranked = ranked.sort_values(["_risk_rank", "_read_rank"], ascending=[True, False])
        sample_qc = ranked.drop_duplicates("snp_code").copy()
        sample_qc = sample_qc[["snp_code"] + [col for col in RNA_QC_COLUMNS if col in sample_qc.columns]]
    bam_summary = inventory.groupby("RNA_BAM_Group", dropna=False).agg(
        n_bams=("RNA_BAM_Name", "count"),
        matched_expression=("matched_expression_sample", "sum"),
        matched_primary=("matched_primary_sample", "sum"),
        indexed_bams=("RNA_BAI_Present", "sum"),
        new_indices=("RNA_BAI_Created", "sum"),
        completed_qc=("RNA_Coverage_QC_Status", lambda s: int((pd.Series(s) == "Complete").sum())),
        unavailable_qc=("RNA_Coverage_QC_Status", lambda s: int((pd.Series(s) != "Complete").sum())),
        high_risk_bams=("RNA_Coverage_Risk_Level", lambda s: int((pd.Series(s) == "High").sum())),
        moderate_risk_bams=("RNA_Coverage_Risk_Level", lambda s: int((pd.Series(s) == "Moderate").sum())),
        unavailable_risk_bams=("RNA_Coverage_Risk_Level", lambda s: int((pd.Series(s) == "Unavailable").sum())),
        mean_panel_depth=("RNA_Panel_Mean_Depth", "mean"),
        mean_min_target_reads=("RNA_Panel_Min_Target_Reads", "mean"),
        mean_gsdmb_depth=("RNA_GSDMB_Mean_Depth", "mean"),
    ).reset_index()
    bam_summary["rna_bed"] = str(rna_bed)
    bam_summary["validation_note"] = "RNA BAM metrics were rebuilt from cached target coverage and refreshed with the current sample parser."
    return inventory, target_cov, bam_summary, sample_qc, panel_design


def build_dna_baseline_inventory(root: Path, master: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    dna_master = (master[master["nucleic_acid"] == "DNA"]
                  .sort_values(["snp_code", "sample_id", "case_id"])
                  .groupby("snp_code", dropna=False)
                  .agg(first)
                  .reset_index(drop=False))
    alias_map: dict[str, set[str]] = {}
    for _, row in dna_master.iterrows():
        for alias in aliases(row.get("snp_code")) | aliases(row.get("sample_id")) | aliases(row.get("case_id")):
            alias_map.setdefault(alias, set()).add(row["snp_code"])

    manifest_map = [
        ("Breast_Tumour", root / "manifests" / "breast-tumour-pass_manifest.txt"),
        ("Breast_Normal", root / "manifests" / "breast-normal-pass_manifest.txt"),
        ("Endometrial_Tumour", root / "manifests" / "endometrium-tumour-pass_manifest.txt"),
        ("Endometrial_Normal", root / "manifests" / "endometrium-normal-pass_manifest.txt"),
    ]
    rows: list[dict[str, object]] = []
    for group_label, manifest_path in manifest_map:
        if not manifest_path.exists():
            continue
        for line in manifest_path.read_text().splitlines():
            bam_path = Path(line.strip())
            if not str(bam_path).strip():
                continue
            stem = bam_path.stem
            candidates = set()
            code = sx(stem)
            if code and code in set(dna_master["snp_code"]):
                candidates.add(code)
            if not candidates:
                for alias in aliases(stem):
                    candidates.update(alias_map.get(alias, set()))
            rows.append({
                "DNA_BAM_Group": group_label,
                "DNA_BAM_Path": str(bam_path),
                "DNA_BAM_Name": bam_path.name,
                "DNA_BAM_Size_MB": round(bam_path.stat().st_size / (1024 * 1024), 2) if bam_path.exists() else pd.NA,
                "DNA_Manifest_Source": manifest_path.name,
                "DNA_BAM_On_Disk": bam_path.exists(),
                "snp_code": next(iter(candidates)) if len(candidates) == 1 else pd.NA,
                "DNA_Match_Status": "matched" if len(candidates) == 1 else ("ambiguous" if len(candidates) > 1 else "unmatched"),
                "DNA_Match_Candidates": "; ".join(sorted(candidates)) if len(candidates) > 1 else "",
            })
    inventory = pd.DataFrame(rows)
    if inventory.empty:
        note = pd.DataFrame({"Note": ["No DNA pass-manifest BAMs were found for baseline anchoring."]})
        return inventory, note
    inventory["DNA_Baseline_Present"] = inventory["DNA_Match_Status"].eq("matched") & inventory["DNA_BAM_On_Disk"].eq(True)
    summary = inventory.groupby("DNA_BAM_Group", dropna=False).agg(
        manifest_rows=("DNA_BAM_Name", "count"),
        matched_baseline_bams=("DNA_Baseline_Present", "sum"),
        matched_unique_codes=("snp_code", lambda s: int(pd.Series(s).dropna().nunique())),
        unmatched_bams=("DNA_Match_Status", lambda s: int((pd.Series(s) == "unmatched").sum())),
        ambiguous_bams=("DNA_Match_Status", lambda s: int((pd.Series(s) == "ambiguous").sum())),
    ).reset_index()
    summary["baseline_note"] = "DNA baseline is anchored to current DNA QC pass manifests, then paired against matched RNA BAMs."
    return inventory, summary


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Standalone RNA QC stage for objective 2")
    ap.add_argument("--expr-xlsx", default=str(D["expr_xlsx"]))
    ap.add_argument("--master", default=str(D["master"]))
    ap.add_argument("--rna-bed", default=str(D.get("rna_bed", "")))
    ap.add_argument("--machine-qc-xlsx", default=str(D.get("machine_qc_xlsx", "")))
    ap.add_argument("--recompute-coverage", action="store_true", help="Ignore cached RNA target coverage and rerun BAM-level coverage QC.")
    ap.add_argument("--out-dir", default=str(D.get("out_dir", Path(__file__).resolve().parent / "analysis_results/19b_objective2_rna_qc")))
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    expr_path = Path(args.expr_xlsx)
    master_path = Path(args.master)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bed_path = Path(args.rna_bed) if args.rna_bed else None
    machine_qc_path = Path(args.machine_qc_xlsx) if args.machine_qc_xlsx else None

    for path, label in [(expr_path, "Expression workbook"), (master_path, "Harmonised master")]:
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    expr = load_expr(expr_path)
    master = load_master(master_path)
    machine_qc = load_machine_qc(machine_qc_path)
    merged, audit = match_expr(expr, master)
    dna_inventory, dna_summary = build_dna_baseline_inventory(Path(__file__).resolve().parent, master)
    cached = None if args.recompute_coverage else rebuild_cached_bam_validation(out_dir, bed_path, merged, sx)
    if cached is None:
        inventory, target_cov, bam_summary, sample_qc, panel_design = bam_validation(Path(__file__).resolve().parent, merged, bed_path, sx)
        cache_mode = "recomputed"
    else:
        inventory, target_cov, bam_summary, sample_qc, panel_design = cached
        cache_mode = "cached"
    merged_qc = attach_dna_baseline(merged, dna_inventory)
    merged_qc = attach_rna_qc(merged_qc, sample_qc)
    merged_qc = attach_machine_qc(merged_qc, machine_qc)
    sample_qc = attach_machine_qc(sample_qc, machine_qc)
    manifest = build_analysis_manifest(merged_qc)
    qc_summary = summarise_analysis_manifest(manifest)
    exclusion_summary = summarise_exclusion_reasons(manifest)
    machine_summary = summarise_machine_qc(manifest)

    manifest_path = out_dir / "RNA_QC_Manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    qc_summary.to_csv(out_dir / "RNA_QC_Summary.tsv", sep="\t", index=False)
    exclusion_summary.to_csv(out_dir / "RNA_QC_Exclusion_Reasons.tsv", sep="\t", index=False)
    machine_summary.to_csv(out_dir / "RNA_Machine_QC_Summary.tsv", sep="\t", index=False)
    machine_qc.to_csv(out_dir / "RNA_Machine_QC.tsv", sep="\t", index=False)
    inventory.to_csv(out_dir / "RNA_BAM_Inventory.tsv", sep="\t", index=False)
    dna_inventory.to_csv(out_dir / "DNA_Baseline_Inventory.tsv", sep="\t", index=False)
    dna_summary.to_csv(out_dir / "DNA_Baseline_Summary.tsv", sep="\t", index=False)
    audit.to_csv(out_dir / "RNA_Match_Audit.tsv", sep="\t", index=False)
    bam_summary.to_csv(out_dir / "RNA_BAM_Validation_Summary.tsv", sep="\t", index=False)
    sample_qc.to_csv(out_dir / "RNA_Sample_QC.tsv", sep="\t", index=False)
    if not target_cov.empty:
        target_cov.to_csv(out_dir / "RNA_Target_Coverage.tsv.gz", sep="\t", index=False, compression="gzip")
    if not panel_design.empty:
        panel_design.to_csv(out_dir / "RNA_Panel_Design.tsv", sep="\t", index=False)

    plot_rna_qc(sample_qc, out_dir)
    legacy_plot = out_dir / "20_RNA_Coverage_QC.png"
    if legacy_plot.exists():
        shutil.move(str(legacy_plot), out_dir / "19b_RNA_Coverage_QC.png")
    plot_rna_gate_summary(qc_summary, out_dir / "19b_RNA_QC_Pass_Summary.png")
    plot_rna_exclusion_summary(exclusion_summary, out_dir / "19b_RNA_QC_Exclusion_Reasons.png")
    plot_rna_machine_vs_panel_qc(manifest, out_dir / "19b_RNA_Machine_vs_Panel_QC.png")

    workbook_path = out_dir / "19b_Objective2_RNA_QC.xlsx"
    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        manifest.to_excel(writer, sheet_name="rna_qc_manifest", index=False)
        qc_summary.to_excel(writer, sheet_name="rna_qc_summary", index=False)
        exclusion_summary.to_excel(writer, sheet_name="rna_qc_exclusions", index=False)
        (machine_summary if not machine_summary.empty else pd.DataFrame({"Note": ["No matched machine-level RNA QC rows were available."]})).to_excel(writer, sheet_name="machine_qc_summary", index=False)
        (machine_qc if not machine_qc.empty else pd.DataFrame({"Note": ["No machine-level RNA QC workbook rows were matched to study samples."]})).to_excel(writer, sheet_name="machine_qc", index=False)
        inventory.to_excel(writer, sheet_name="bam_inventory", index=False)
        dna_inventory.to_excel(writer, sheet_name="dna_baseline_inventory", index=False)
        dna_summary.to_excel(writer, sheet_name="dna_baseline_summary", index=False)
        bam_summary.to_excel(writer, sheet_name="bam_validation_summary", index=False)
        audit.to_excel(writer, sheet_name="match_audit", index=False)
        (sample_qc if not sample_qc.empty else pd.DataFrame({"Note": ["No sample-level RNA QC rows were available."]})).to_excel(writer, sheet_name="sample_qc", index=False)
        (target_cov if not target_cov.empty else pd.DataFrame({"Note": ["No target-level RNA coverage rows were available."]})).to_excel(writer, sheet_name="target_coverage", index=False)
        (panel_design if not panel_design.empty else pd.DataFrame({"Note": ["No RNA panel design rows were parsed."]})).to_excel(writer, sheet_name="panel_design", index=False)

    baseline_manifest = manifest[manifest["DNA_Baseline_Present"].fillna(False)].copy() if not manifest.empty and "DNA_Baseline_Present" in manifest.columns else manifest.copy()
    total_rows = len(baseline_manifest)
    strict_rows = int(baseline_manifest["RNA_QC_Analysis_Ready"].sum()) if not baseline_manifest.empty else 0
    exploratory_rows = int(baseline_manifest["RNA_QC_Exploratory_Ready"].sum()) if not baseline_manifest.empty else 0
    exploratory_only_rows = max(exploratory_rows - strict_rows, 0)
    primary_strict = int(baseline_manifest.loc[baseline_manifest["analysis_include_primary"].fillna(False), "RNA_QC_Analysis_Ready"].sum()) if not baseline_manifest.empty else 0
    primary_exploratory = int(baseline_manifest.loc[baseline_manifest["analysis_include_primary"].fillna(False), "RNA_QC_Exploratory_Ready"].sum()) if not baseline_manifest.empty else 0
    excluded_rows = total_rows - exploratory_rows
    print("=== RNA QC stage summary ===")
    print(f"DNA-baseline rows reviewed   : {total_rows}")
    print(f"Matched master rows          : {int(baseline_manifest['master_join_success'].fillna(False).sum()) if not baseline_manifest.empty else 0}")
    print(f"Strict-ready rows            : {strict_rows}")
    print(f"Exploratory-ready rows       : {exploratory_rows}")
    print(f"Exploratory-only rows        : {exploratory_only_rows}")
    print(f"Excluded rows                : {excluded_rows}")
    print(f"Primary tumour rows strict   : {primary_strict}")
    print(f"Primary tumour rows usable   : {primary_exploratory}")
    matched_machine_rows = int((baseline_manifest["RNA_Machine_QC_Flag"].astype(str) != "Unmatched").sum()) if not baseline_manifest.empty and "RNA_Machine_QC_Flag" in baseline_manifest.columns else 0
    print(f"RNA BAM inventory rows       : {len(inventory)}")
    print(f"DNA baseline inventory rows  : {len(dna_inventory)}")
    print(f"Machine QC matched rows      : {matched_machine_rows}")
    print(f"Manifest                     : {manifest_path}")
    print(f"Workbook                     : {workbook_path}")
    print(f"Coverage source             : {cache_mode}")
    print(f"Output directory             : {out_dir}")


if __name__ == "__main__":
    main()

