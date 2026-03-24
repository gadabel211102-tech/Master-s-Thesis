#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import pandas as pd

from association_runtime import script20_defaults
from objective2_rna_qc import (
    RNA_QC_COLUMNS,
    attach_rna_qc,
    bam_validation,
    build_analysis_manifest,
    plot_rna_exclusion_summary,
    plot_rna_gate_summary,
    plot_rna_qc,
    summarise_analysis_manifest,
    summarise_exclusion_reasons,
)

D = script20_defaults()
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


def sx(v):
    text = s(v)
    patterns = [
        (r"^(SNP_(?:AT|EN|MN|MT-T)_\d+)$", lambda m: m.group(1).upper()),
        (r"^SNP_DNA_(EN|MN|AT)_(\d+)$", lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^SNP_DNA_MT[-_]T_(\d+)$", lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^DNA_SNP_(EN|MN)_(\d+)", lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^DNA_EN_(\d+)", lambda m: f"SNP_EN_{m.group(1)}"),
        (r"^DNA_MN_(\d+)", lambda m: f"SNP_MN_{m.group(1)}"),
        (r"^DNA_SNP_AT_(\d+)", lambda m: f"SNP_AT_{m.group(1)}"),
        (r"^DNA_(?:SNP_)?MT[-_]T_(\d+)", lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^RNA_SNP_(AT|EN|MN|MT-T)_(\d+)", lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^SNP_RNA_(AT|EN|MN|MT-T)_(\d+)$", lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^RNA_(AT|EN|MN)_(\d+)", lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^RNA_MT[-_]T_(\d+)", lambda m: f"SNP_MT-T_{m.group(1)}"),
    ]
    for pat, fmt in patterns:
        match = re.match(pat, text, re.I)
        if match:
            return fmt(match)
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


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Standalone RNA QC stage for objective 2")
    ap.add_argument("--expr-xlsx", default=str(D["expr_xlsx"]))
    ap.add_argument("--master", default=str(D["master"]))
    ap.add_argument("--rna-bed", default=str(D.get("rna_bed", "")))
    ap.add_argument("--out-dir", default="/home/gadeaalonsoj/tfm/analysis_results/19b_objective2_rna_qc")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    expr_path = Path(args.expr_xlsx)
    master_path = Path(args.master)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bed_path = Path(args.rna_bed) if args.rna_bed else None

    for path, label in [(expr_path, "Expression workbook"), (master_path, "Harmonised master")]:
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    expr = load_expr(expr_path)
    master = load_master(master_path)
    merged, audit = match_expr(expr, master)
    inventory, target_cov, bam_summary, sample_qc, panel_design = bam_validation(Path(__file__).resolve().parent, merged, bed_path, sx)
    merged_qc = attach_rna_qc(merged, sample_qc)
    manifest = build_analysis_manifest(merged_qc)
    qc_summary = summarise_analysis_manifest(manifest)
    exclusion_summary = summarise_exclusion_reasons(manifest)

    manifest_path = out_dir / "RNA_QC_Manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    qc_summary.to_csv(out_dir / "RNA_QC_Summary.tsv", sep="\t", index=False)
    exclusion_summary.to_csv(out_dir / "RNA_QC_Exclusion_Reasons.tsv", sep="\t", index=False)
    inventory.to_csv(out_dir / "RNA_BAM_Inventory.tsv", sep="\t", index=False)
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

    workbook_path = out_dir / "19b_Objective2_RNA_QC.xlsx"
    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        manifest.to_excel(writer, sheet_name="rna_qc_manifest", index=False)
        qc_summary.to_excel(writer, sheet_name="rna_qc_summary", index=False)
        exclusion_summary.to_excel(writer, sheet_name="rna_qc_exclusions", index=False)
        inventory.to_excel(writer, sheet_name="bam_inventory", index=False)
        bam_summary.to_excel(writer, sheet_name="bam_validation_summary", index=False)
        audit.to_excel(writer, sheet_name="match_audit", index=False)
        (sample_qc if not sample_qc.empty else pd.DataFrame({"Note": ["No sample-level RNA QC rows were available."]})).to_excel(writer, sheet_name="sample_qc", index=False)
        (target_cov if not target_cov.empty else pd.DataFrame({"Note": ["No target-level RNA coverage rows were available."]})).to_excel(writer, sheet_name="target_coverage", index=False)
        (panel_design if not panel_design.empty else pd.DataFrame({"Note": ["No RNA panel design rows were parsed."]})).to_excel(writer, sheet_name="panel_design", index=False)

    total_rows = len(manifest)
    passed_rows = int(manifest["RNA_QC_Analysis_Ready"].sum()) if not manifest.empty else 0
    primary_passed = int(manifest.loc[manifest["analysis_include_primary"].fillna(False), "RNA_QC_Analysis_Ready"].sum()) if not manifest.empty else 0
    excluded_rows = total_rows - passed_rows
    print("=== RNA QC stage summary ===")
    print(f"Expression rows reviewed      : {total_rows}")
    print(f"Matched master rows          : {int(manifest['master_join_success'].fillna(False).sum()) if not manifest.empty else 0}")
    print(f"QC-passed rows               : {passed_rows}")
    print(f"QC-excluded rows             : {excluded_rows}")
    print(f"Primary tumour rows passed   : {primary_passed}")
    print(f"BAM inventory rows           : {len(inventory)}")
    print(f"Manifest                     : {manifest_path}")
    print(f"Workbook                     : {workbook_path}")
    print(f"Output directory             : {out_dir}")


if __name__ == "__main__":
    main()
