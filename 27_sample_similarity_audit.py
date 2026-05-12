#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from itertools import combinations
from pathlib import Path
from typing import Iterable

import pandas as pd

from pipeline_utils import get_paths
from sample_identity_utils import build_analysis_sample_map

CALLABLE_CLASSES = {"WT", "Het_ALT", "Hom_ALT", "Other_ALT_Model", "Other_ALT_model"}
NONREF_CLASSES = {"Het_ALT", "Hom_ALT", "Other_ALT_Model", "Other_ALT_model"}
STATUS_RANK = {
    "Hom_ALT": 6,
    "Het_ALT": 5,
    "Other_ALT_Model": 4,
    "Other_ALT_model": 4,
    "WT": 3,
    "Below_DP_threshold": 2,
    "No_call": 1,
}


def _extract_snp_code(sample_name: str) -> str | None:
    s = str(sample_name)
    m = re.match(r"^(SNP_MT-T_\d+_REP)", s, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    m = re.match(r"^(SNP_(?:AT|EN|MN|MT-T)_\d+)", s, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    m = re.match(r"^DNA_MT[-_]T_(\d+)_Repeticion_", s, re.IGNORECASE)
    if m:
        return f"SNP_MT-T_{m.group(1)}_REP"
    m = re.match(r"^DNA_SNP_(EN|MN)_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_{m.group(1).upper()}_{m.group(2)}"
    m = re.match(r"^DNA_SNP_AT_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_AT_{int(m.group(1))}"
    m = re.match(r"^DNA_AT_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_AT_{int(m.group(1))}"
    m = re.match(r"^DNA_SNP_(?:CK|RSB)_ECLAI_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_AT_{int(m.group(1))}"
    m = re.match(r"^SNP_DNA_AT_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_AT_{m.group(1)}"
    m = re.match(r"^SNP_DNA_(MN|EN)_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_{m.group(1).upper()}_{m.group(2)}"
    m = re.match(r"^DNA_SNP_MT-T_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_MT-T_{m.group(1)}"
    m = re.match(r"^DNA_MT[-_]T_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_MT-T_{m.group(1)}"
    return None


def _find_col(columns: Iterable[str], *targets: str) -> str:
    normalized = {str(col).strip().upper(): str(col) for col in columns}
    for target in targets:
        hit = normalized.get(str(target).strip().upper())
        if hit is not None:
            return hit
    raise KeyError(f"Missing required column among: {targets}")


def _first_nonempty(series: pd.Series):
    for value in series:
        if pd.notna(value) and str(value).strip() and str(value).strip().lower() != "nan":
            return value
    return pd.NA


def load_master_metadata(master_path: Path) -> pd.DataFrame:
    master = pd.read_excel(master_path, sheet_name="harmonised_plus_canon")
    keep = [col for col in ["snp_code", "sample_id", "case_id", "sheet", "tissue", "tumour_normal"] if col in master.columns]
    meta = master[keep].copy()
    meta["snp_code"] = meta["snp_code"].astype(str).str.strip().replace({"": pd.NA, "nan": pd.NA, "None": pd.NA})
    meta = meta[meta["snp_code"].notna()].copy()
    return meta.groupby("snp_code", as_index=False).agg({col: _first_nonempty for col in keep if col != "snp_code"})


def load_identity_map(raw_sample_names: list[str], phased_path: Path, master_meta: pd.DataFrame) -> pd.DataFrame:
    identity = build_analysis_sample_map(raw_sample_names, phased_path, _extract_snp_code)
    raw_name = identity["raw_sample_name"].astype(str)
    identity["name_repetition_note"] = raw_name.str.contains(r"repeticion|_rep", case=False, regex=True).map(
        {True: "Name_contains_repetition_text_only", False: "No_repetition_text_in_name"}
    )
    identity["within_code_exact_note"] = identity["collapsed_exact_alias_group"].map(
        {True: "Exact_within_code_genotype_alias", False: "No_exact_within_code_alias"}
    )
    identity = identity.merge(master_meta, on="snp_code", how="left")
    return identity


def load_locus_frequency_classes(annot_path: Path) -> pd.DataFrame:
    variant = pd.read_excel(annot_path, sheet_name="Variant_Annotations")
    locus_col = _find_col(variant.columns, "Locus_Key")
    af_col = _find_col(variant.columns, "gnomAD_NFE_AF_combined", "gnomADe_NFE_AF", "gnomAD_NFE_AF")
    work = variant[[locus_col, af_col]].copy()
    work.columns = ["Locus_Key", "AF"]
    work["AF"] = pd.to_numeric(work["AF"], errors="coerce")
    grouped = work.groupby("Locus_Key", as_index=False).agg(
        all_common=("AF", lambda s: bool(len(s)) and s.notna().all() and (s > 0.01).all()),
        all_rare=("AF", lambda s: bool(len(s)) and ((s.isna()) | (s <= 0.01)).all()),
    )
    grouped["locus_frequency_class"] = "mixed"
    grouped.loc[grouped["all_common"], "locus_frequency_class"] = "common"
    grouped.loc[grouped["all_rare"], "locus_frequency_class"] = "rare"
    return grouped[["Locus_Key", "locus_frequency_class"]]


def load_genotype_matrix(genotype_path: Path) -> pd.DataFrame:
    geno = pd.read_csv(genotype_path, sep="\t", dtype=str)
    sample_col = _find_col(geno.columns, "Sample")
    locus_col = _find_col(geno.columns, "Locus_Key")
    class_col = _find_col(geno.columns, "Raw_Genotype_Class", "Genotype_Status")
    work = geno[[sample_col, locus_col, class_col]].copy()
    work.columns = ["Sample", "Locus_Key", "Genotype_Class"]
    work["rank"] = work["Genotype_Class"].map(STATUS_RANK).fillna(0)
    work = work.sort_values(["Sample", "Locus_Key", "rank"], ascending=[True, True, False])
    work = work.drop_duplicates(subset=["Sample", "Locus_Key"], keep="first")
    return work.pivot(index="Locus_Key", columns="Sample", values="Genotype_Class")


def _pair_metrics(sample_a: str, sample_b: str, matrix: pd.DataFrame, locus_classes: pd.Series) -> dict[str, object]:
    a = matrix[sample_a]
    b = matrix[sample_b]
    callable_mask = a.isin(CALLABLE_CLASSES) & b.isin(CALLABLE_CLASSES)

    def concordance_for(mask: pd.Series) -> tuple[int, int, int, float | pd.NA]:
        comparable = int(mask.sum())
        if comparable == 0:
            return 0, 0, 0, pd.NA
        exact = int((a[mask] == b[mask]).sum())
        mismatch = comparable - exact
        return comparable, exact, mismatch, exact / comparable

    common_mask = callable_mask & locus_classes.eq("common")
    rare_mask = callable_mask & locus_classes.eq("rare")

    comparable_all, exact_all, mismatch_all, concordance_all = concordance_for(callable_mask)
    comparable_common, exact_common, mismatch_common, concordance_common = concordance_for(common_mask)
    comparable_rare, exact_rare, mismatch_rare, concordance_rare = concordance_for(rare_mask)

    nonref_a = a.isin(NONREF_CLASSES)
    nonref_b = b.isin(NONREF_CLASSES)
    common_locus = locus_classes.eq("common")
    rare_locus = locus_classes.eq("rare")

    shared_nonref_common = int((nonref_a & nonref_b & common_locus).sum())
    union_nonref_common = int(((nonref_a | nonref_b) & common_locus).sum())
    shared_nonref_rare = int((nonref_a & nonref_b & rare_locus).sum())
    union_nonref_rare = int(((nonref_a | nonref_b) & rare_locus).sum())

    return {
        "sample_1": sample_a,
        "sample_2": sample_b,
        "comparable_all": comparable_all,
        "concordance_all": concordance_all,
        "exact_all": exact_all,
        "mismatch_all": mismatch_all,
        "comparable_common": comparable_common,
        "concordance_common": concordance_common,
        "exact_common": exact_common,
        "mismatch_common": mismatch_common,
        "comparable_rare": comparable_rare,
        "concordance_rare": concordance_rare,
        "exact_rare": exact_rare,
        "mismatch_rare": mismatch_rare,
        "shared_nonref_common": shared_nonref_common,
        "union_nonref_common": union_nonref_common,
        "common_nonref_jaccard": (shared_nonref_common / union_nonref_common) if union_nonref_common else pd.NA,
        "shared_nonref_rare": shared_nonref_rare,
        "union_nonref_rare": union_nonref_rare,
        "rare_nonref_jaccard": (shared_nonref_rare / union_nonref_rare) if union_nonref_rare else pd.NA,
    }


def classify_pair(row: pd.Series) -> str:
    if bool(row.get("same_snp_code", False)):
        if bool(row.get("same_code_but_distinct_genotype", False)):
            return "SAME_CODE_DISTINCT_GENOTYPE_KEEP_SEPARATE"
        if bool(row.get("same_phased_digest", False)):
            return "SAME_CODE_EXACT_ALIAS_WITHIN_CURRENT_RULE"
    concordance = pd.to_numeric(row.get("concordance_all"), errors="coerce")
    if pd.notna(concordance) and concordance == 1.0 and not bool(row.get("same_snp_code", False)):
        return "DIFFERENT_CODE_IDENTICAL_REVIEW_ONLY_KEEP_SEPARATE"
    if pd.notna(concordance) and concordance >= 0.99 and not bool(row.get("same_snp_code", False)):
        return "DIFFERENT_CODE_NEAR_IDENTICAL_REVIEW_ONLY_KEEP_SEPARATE"
    return "DISTINCT_SAMPLES"


def build_outputs(out_dir: Path, genotype_path: Path, annot_path: Path, master_path: Path, phased_path: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    matrix = load_genotype_matrix(genotype_path)
    raw_samples = [str(col) for col in matrix.columns]
    master_meta = load_master_metadata(master_path)
    identity = load_identity_map(raw_samples, phased_path, master_meta)
    locus_df = load_locus_frequency_classes(annot_path)
    locus_classes = locus_df.set_index("Locus_Key")["locus_frequency_class"].reindex(matrix.index).fillna("mixed")

    pair_rows = [_pair_metrics(a, b, matrix, locus_classes) for a, b in combinations(raw_samples, 2)]
    pairs = pd.DataFrame(pair_rows)

    meta_cols = [
        "raw_sample_name", "snp_code", "analysis_sample_id", "identity_root", "phased_digest",
        "group_order", "exact_group_size", "n_exact_groups_for_code", "collapsed_exact_alias_group",
        "in_phased_backbone", "name_repetition_note", "within_code_exact_note", "sample_id",
        "case_id", "sheet", "tissue", "tumour_normal",
    ]
    sample_meta = identity[meta_cols].copy()

    pairs = pairs.merge(sample_meta.add_suffix("_1"), left_on="sample_1", right_on="raw_sample_name_1", how="left")
    pairs = pairs.merge(sample_meta.add_suffix("_2"), left_on="sample_2", right_on="raw_sample_name_2", how="left")

    pairs["same_snp_code"] = pairs["snp_code_1"].fillna("").eq(pairs["snp_code_2"].fillna("")) & pairs["snp_code_1"].notna()
    pairs["same_case_id"] = pairs["case_id_1"].fillna("").eq(pairs["case_id_2"].fillna("")) & pairs["case_id_1"].notna()
    pairs["same_phased_digest"] = pairs["phased_digest_1"].fillna("").eq(pairs["phased_digest_2"].fillna("")) & pairs["phased_digest_1"].notna()
    pairs["same_code_but_distinct_genotype"] = pairs["same_snp_code"] & (~pairs["same_phased_digest"])
    pairs["classification"] = pairs.apply(classify_pair, axis=1)

    review_classes = {
        "SAME_CODE_DISTINCT_GENOTYPE_KEEP_SEPARATE",
        "DIFFERENT_CODE_IDENTICAL_REVIEW_ONLY_KEEP_SEPARATE",
        "DIFFERENT_CODE_NEAR_IDENTICAL_REVIEW_ONLY_KEEP_SEPARATE",
    }
    review = pairs[pairs["classification"].isin(review_classes)].copy()
    review = review.sort_values(["classification", "concordance_all", "shared_nonref_rare", "comparable_all"], ascending=[True, False, False, False])

    identity.to_csv(out_dir / "28_sample_identity_map.tsv", sep="\t", index=False)
    pairs.to_csv(out_dir / "28_pairwise_sample_similarity.tsv", sep="\t", index=False)
    review.to_csv(out_dir / "28_pairs_for_manual_review.tsv", sep="\t", index=False)

    summary_rows = [
        {"metric": "raw_samples", "value": int(identity["raw_sample_name"].nunique())},
        {"metric": "unique_snp_code", "value": int(identity["snp_code"].nunique(dropna=True))},
        {"metric": "analysis_sample_id", "value": int(identity["analysis_sample_id"].nunique())},
        {"metric": "pairwise_comparisons", "value": int(len(pairs))},
    ]
    for label, count in pairs["classification"].value_counts().items():
        summary_rows.append({"metric": label, "value": int(count)})
    pd.DataFrame(summary_rows).to_csv(out_dir / "28_similarity_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    paths = get_paths()
    parser = argparse.ArgumentParser(description="Audit pairwise DNA sample similarity without collapsing different formal sample identities.")
    parser.add_argument("--genotypes", default=str(paths["forced_genotypes_dir"] / "05b_per_sample_force_genotypes.tsv"))
    parser.add_argument("--annot", default=str(paths["annotated_report"]))
    parser.add_argument("--master", default=str(paths["harmonised_master"]))
    parser.add_argument("--phased", default=str(paths["haplotype_phased"]))
    parser.add_argument("--out-dir", default=str(paths["results_dir"] / "28_sample_similarity_audit"))
    args = parser.parse_args()
    build_outputs(Path(args.out_dir), Path(args.genotypes), Path(args.annot), Path(args.master), Path(args.phased))
