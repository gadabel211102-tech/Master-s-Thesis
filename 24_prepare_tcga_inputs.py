#!/usr/bin/env python3
"""Prepare TCGA BRCA/UCEC metadata and GDC manifests for stage-24 cross-validation.

This helper does the public part of the TCGA setup automatically:

- queries the GDC API for controlled-access TCGA genotyping-array germline files
- applies the same female and age-window logic used by stage 24
- chooses one representative germline file per case
- writes stage-24-compatible metadata TSVs plus GDC manifest TSVs

The remaining controlled-access download still requires an authorised GDC/dbGaP
user token. Public GDC metadata can be prepared without that token, which makes
it much easier to see what is missing locally and what needs to be downloaded.
"""

from __future__ import annotations

import argparse
import json
import math
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from pipeline_utils import ensure_directory

GDC_API = "https://api.gdc.cancer.gov"
PAGE_SIZE = 5000

PROJECTS = {
    "TCGA-BRCA": {
        "comparison": "Breast",
        "short_label": "BRCA",
    },
    "TCGA-UCEC": {
        "comparison": "Endometrium",
        "short_label": "UCEC",
    },
}


@dataclass(frozen=True)
class PreparedProject:
    project_id: str
    comparison: str
    short_label: str
    all_files_path: Path
    metadata_path: Path
    manifest_path: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare public TCGA metadata tables and controlled-access GDC manifests "
            "for the stage-24 external cross-validation branch."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--projects",
        nargs="+",
        default=["TCGA-BRCA", "TCGA-UCEC"],
        choices=sorted(PROJECTS),
        help="TCGA projects to prepare.",
    )
    parser.add_argument("--age-min", type=int, default=45, help="Minimum diagnosis age in years.")
    parser.add_argument("--age-max", type=int, default=65, help="Maximum diagnosis age in years.")
    parser.add_argument(
        "--out-dir",
        default="/home/gadeaalonsoj/tfm/analysis_results/24_external_cross_validation/tcga_prep",
        help="Directory for audit tables and the workbook summary.",
    )
    parser.add_argument(
        "--manifest-dir",
        default="/home/gadeaalonsoj/tfm/manifests",
        help="Directory where GDC manifest TSVs will be written.",
    )
    return parser.parse_args()


def gdc_get(endpoint: str, params: dict[str, object]) -> dict[str, object]:
    query = urllib.parse.urlencode(params)
    url = f"{GDC_API}/{endpoint}?{query}"
    with urllib.request.urlopen(url, timeout=120) as response:
        return json.load(response)


def build_filters(project_id: str) -> dict[str, object]:
    return {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project_id]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Simple Germline Variation"]}},
            {"op": "in", "content": {"field": "experimental_strategy", "value": ["Genotyping Array"]}},
        ],
    }


def fetch_project_hits(project_id: str) -> list[dict[str, object]]:
    fields = ",".join(
        [
            "file_id",
            "file_name",
            "md5sum",
            "file_size",
            "state",
            "access",
            "data_category",
            "data_type",
            "experimental_strategy",
            "cases.submitter_id",
            "cases.demographic.gender",
            "cases.diagnoses.age_at_diagnosis",
            "cases.samples.submitter_id",
            "cases.samples.sample_type",
            "cases.samples.tissue_type",
        ]
    )
    data = gdc_get(
        "files",
        {
            "filters": json.dumps(build_filters(project_id)),
            "fields": fields,
            "format": "JSON",
            "size": str(PAGE_SIZE),
        },
    )
    return list(data["data"]["hits"])


def age_years_from_case(case: dict[str, object]) -> float:
    diagnosis_rows = case.get("diagnoses") or []
    values = []
    for row in diagnosis_rows:
        age_days = row.get("age_at_diagnosis")
        if age_days is None:
            continue
        try:
            values.append(float(age_days) / 365.25)
        except (TypeError, ValueError):
            continue
    return min(values) if values else math.nan


def flatten_hits(project_id: str, hits: list[dict[str, object]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for hit in hits:
        cases = hit.get("cases") or []
        if not cases:
            continue
        case = cases[0]
        demographic = case.get("demographic") or {}
        samples = case.get("samples") or []
        sample_ids = sorted({str(sample.get("submitter_id", "")).strip() for sample in samples if str(sample.get("submitter_id", "")).strip()})
        sample_types = sorted({str(sample.get("sample_type", "")).strip() for sample in samples if str(sample.get("sample_type", "")).strip()})
        tissue_types = sorted({str(sample.get("tissue_type", "")).strip() for sample in samples if str(sample.get("tissue_type", "")).strip()})
        case_id = str(case.get("submitter_id", "")).strip()
        if not case_id:
            continue
        rows.append(
            {
                "Project": project_id,
                "Case_Submitter_ID": case_id,
                "Gender": str(demographic.get("gender", "")).strip().lower(),
                "Age_Years": age_years_from_case(case),
                "File_ID": hit.get("file_id") or hit.get("id"),
                "File_Name": hit.get("file_name"),
                "MD5": hit.get("md5sum"),
                "File_Size": pd.to_numeric(hit.get("file_size"), errors="coerce"),
                "State": hit.get("state"),
                "Access": hit.get("access"),
                "Data_Category": hit.get("data_category"),
                "Data_Type": hit.get("data_type"),
                "Experimental_Strategy": hit.get("experimental_strategy"),
                "Associated_Sample_IDs": "; ".join(sample_ids),
                "Associated_Sample_Types": "; ".join(sample_types),
                "Associated_Tissue_Types": "; ".join(tissue_types),
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values(["Case_Submitter_ID", "File_Name", "File_ID"]).reset_index(drop=True)


def choose_representative_files(df: pd.DataFrame, age_min: int, age_max: int) -> pd.DataFrame:
    if df.empty:
        return df
    filtered = df[df["Gender"].eq("female")].copy()
    filtered = filtered[filtered["Age_Years"].between(age_min, age_max, inclusive="both")].copy()
    if filtered.empty:
        return filtered

    filtered["Duplicate_File_Count"] = filtered.groupby("Case_Submitter_ID")["File_ID"].transform("count")
    filtered["All_File_Names_For_Case"] = filtered.groupby("Case_Submitter_ID")["File_Name"].transform(lambda s: "; ".join(sorted(set(s.astype(str)))))
    filtered["All_File_IDs_For_Case"] = filtered.groupby("Case_Submitter_ID")["File_ID"].transform(lambda s: "; ".join(sorted(set(s.astype(str)))))
    selected = filtered.sort_values(["Case_Submitter_ID", "File_Name", "File_ID"]).groupby("Case_Submitter_ID", as_index=False).head(1).copy()
    selected["Sample"] = selected["Case_Submitter_ID"]
    selected["gender.demographic"] = "female"
    selected["age"] = selected["Age_Years"].round(2)
    return selected.reset_index(drop=True)


def build_stage24_metadata(selected: pd.DataFrame) -> pd.DataFrame:
    if selected.empty:
        return pd.DataFrame(columns=["sample", "gender.demographic", "age", "project", "case_submitter_id", "file_id", "file_name", "duplicate_germline_file_count", "all_candidate_file_names"])
    out = selected[
        [
            "Sample",
            "gender.demographic",
            "age",
            "Project",
            "Case_Submitter_ID",
            "File_ID",
            "File_Name",
            "Duplicate_File_Count",
            "All_File_Names_For_Case",
        ]
    ].copy()
    out.columns = [
        "sample",
        "gender.demographic",
        "age",
        "project",
        "case_submitter_id",
        "file_id",
        "file_name",
        "duplicate_germline_file_count",
        "all_candidate_file_names",
    ]
    return out.sort_values(["sample", "file_name"]).reset_index(drop=True)


def build_gdc_manifest(selected: pd.DataFrame) -> pd.DataFrame:
    if selected.empty:
        return pd.DataFrame(columns=["id", "filename", "md5", "size", "state"])
    out = selected[["File_ID", "File_Name", "MD5", "File_Size", "State"]].copy()
    out.columns = ["id", "filename", "md5", "size", "state"]
    out["size"] = out["size"].fillna(0).astype(int)
    return out.sort_values(["filename", "id"]).reset_index(drop=True)


def build_summary_row(project_id: str, all_files: pd.DataFrame, selected: pd.DataFrame, age_min: int, age_max: int) -> dict[str, object]:
    total_cases = all_files["Case_Submitter_ID"].nunique() if not all_files.empty else 0
    female_cases = all_files.loc[all_files["Gender"].eq("female"), "Case_Submitter_ID"].nunique() if not all_files.empty else 0
    selected_cases = selected["Case_Submitter_ID"].nunique() if not selected.empty else 0
    selected_size_gib = float(selected["File_Size"].fillna(0).sum() / 1024 / 1024 / 1024) if not selected.empty else 0.0
    return {
        "Project": project_id,
        "Comparison": PROJECTS[project_id]["comparison"],
        "Controlled_File_Data_Type": "Simple Germline Variation",
        "Controlled_File_Strategy": "Genotyping Array",
        "Total_Germline_Files": int(all_files.shape[0]),
        "Total_Cases_With_Germline_Files": int(total_cases),
        "Female_Cases": int(female_cases),
        "Female_Cases_Age_%d_%d" % (age_min, age_max): int(selected_cases),
        "Representative_Files_Selected": int(selected.shape[0]),
        "Selected_File_Size_GiB": round(selected_size_gib, 2),
        "Token_Required_For_Download": True,
        "Stage24_Metadata_Ready": True,
        "Stage24_VCF_Ready": False,
        "Next_Blocker": "Controlled-access germline-array files must still be downloaded and converted into a stage-24-compatible VCF/BCF.",
    }


def write_readme(out_dir: Path, prepared: list[PreparedProject], age_min: int, age_max: int) -> None:
    lines = [
        "TCGA stage-24 preparation outputs",
        "",
        "What this helper prepared:",
        "- public metadata tables compatible with the stage-24 metadata loader",
        "- GDC manifests listing one representative controlled germline-array file per female case",
        f"- the same age window used by stage 24: {age_min}-{age_max} years",
        "",
        "What still blocks full TCGA cross-validation:",
        "- the TCGA genotype files are controlled-access Simple Germline Variation files",
        "- a GDC/dbGaP token is required to download them",
        "- after download, those files still need conversion into a multi-sample VCF/BCF whose sample names match the `sample` column in the metadata TSV",
        "",
        "Suggested download commands after you obtain a token:",
    ]
    for item in prepared:
        lines.append(
            f"- {item.project_id}: gdc-client download -t /path/to/gdc-user-token.txt -m {item.manifest_path} -d /path/to/{item.short_label.lower()}_germline_download"
        )
        lines.append(
            f"  stage-24 metadata table: {item.metadata_path}"
        )
    lines.extend(
        [
            "",
            "After conversion to VCF/BCF, rerun stage 24 with:",
            "- --tcga-breast-vcf /path/to/TCGA_BRCA_germline.vcf.gz --tcga-breast-metadata /path/to/24_TCGA_BRCA_Metadata.tsv",
            "- --tcga-endometrial-vcf /path/to/TCGA_UCEC_germline.vcf.gz --tcga-endometrial-metadata /path/to/24_TCGA_UCEC_Metadata.tsv",
        ]
    )
    (out_dir / "24_TCGA_PREP_README.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    out_dir = ensure_directory(args.out_dir)
    manifest_dir = ensure_directory(args.manifest_dir)

    summary_rows: list[dict[str, object]] = []
    workbook_tables: dict[str, pd.DataFrame] = {}
    prepared: list[PreparedProject] = []

    for project_id in args.projects:
        cfg = PROJECTS[project_id]
        short_label = cfg["short_label"]
        all_files_path = out_dir / f"24_{short_label}_All_Germline_Array_Files.tsv"
        metadata_path = out_dir / f"24_{short_label}_Metadata.tsv"
        manifest_path = manifest_dir / f"24_{short_label}_Germline_Array_GDC_Manifest.tsv"

        hits = fetch_project_hits(project_id)
        all_files = flatten_hits(project_id, hits)
        selected = choose_representative_files(all_files, args.age_min, args.age_max)
        metadata = build_stage24_metadata(selected)
        manifest = build_gdc_manifest(selected)

        all_files.to_csv(all_files_path, sep="\t", index=False)
        metadata.to_csv(metadata_path, sep="\t", index=False)
        manifest.to_csv(manifest_path, sep="\t", index=False)

        workbook_tables[f"{short_label}_All_Files"] = all_files
        workbook_tables[f"{short_label}_Selected_Cases"] = selected
        workbook_tables[f"{short_label}_Stage24_Metadata"] = metadata
        workbook_tables[f"{short_label}_GDC_Manifest"] = manifest
        summary_rows.append(build_summary_row(project_id, all_files, selected, args.age_min, args.age_max))
        prepared.append(
            PreparedProject(
                project_id=project_id,
                comparison=cfg["comparison"],
                short_label=short_label,
                all_files_path=all_files_path,
                metadata_path=metadata_path,
                manifest_path=manifest_path,
            )
        )

    summary_df = pd.DataFrame(summary_rows)
    workbook_path = out_dir / "24_TCGA_Input_Preparation.xlsx"
    readme_df = pd.DataFrame(
        {
            "Notes": [
                "Purpose: prepare the public TCGA metadata and controlled-download manifests needed for stage-24 cross-validation.",
                "The relevant TCGA file type is controlled-access Simple Germline Variation generated by Genotyping Array.",
                "The helper keeps one representative germline-array file per case after filtering to female cases in the requested age range.",
                "Metadata TSVs are stage-24-compatible now, but the genotype VCF/BCF still has to be produced from the downloaded controlled files.",
                "Exact stage-24 sample matching will work if the future TCGA VCF/BCF uses the `sample` column from the metadata TSV as its header sample IDs.",
            ]
        }
    )
    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        readme_df.to_excel(writer, sheet_name="README", index=False)
        summary_df.to_excel(writer, sheet_name="Summary", index=False)
        for sheet_name, table in workbook_tables.items():
            table.to_excel(writer, sheet_name=sheet_name[:31], index=False)

    write_readme(out_dir, prepared, args.age_min, args.age_max)

    print("Prepared TCGA stage-24 inputs")
    print(f"Workbook: {workbook_path}")
    for row in summary_rows:
        project_id = row["Project"]
        short_label = PROJECTS[project_id]["short_label"]
        print(
            f"{project_id}: {row['Representative_Files_Selected']} representative female cases in the {args.age_min}-{args.age_max} year window; manifest -> {manifest_dir / ('24_' + short_label + '_Germline_Array_GDC_Manifest.tsv')}"
        )


if __name__ == "__main__":
    main()
