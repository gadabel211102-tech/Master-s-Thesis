#!/usr/bin/env python3
"""Collect significant KM/survival curve pages into thesis-friendly PDFs.

This script works from the active survival PDFs themselves so the output bundle
always reflects the curves that were actually rendered, rather than larger
workbook tables that may contain rows not exported to the current PDFs.

Recommended run command in this project:
  /home/gadeaalonsoj/tfm_env/bin/python 29_collect_significant_survival_curves.py
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd
from pypdf import PdfReader, PdfWriter


ROOT = Path(__file__).resolve().parent
RES = ROOT / "analysis_results"
OUT_DIR = RES / "survival_curve_shortlists"

STAGE17_PDFS = [
    RES / "17_snp_clinical_associations" / "supplementary_figures" / "17_KM_Curves_Survival.pdf",
    RES / "17_snp_clinical_associations" / "supplementary_figures" / "17_KM_Curves_Survival_age.pdf",
    RES / "17_snp_clinical_associations" / "supplementary_figures" / "17_KM_Curves_Survival_age_bmi.pdf",
    RES / "17_snp_clinical_associations" / "supplementary_figures" / "17_KM_Curves_Survival_bmi.pdf",
    RES / "17_snp_clinical_associations" / "supplementary_figures" / "17_KM_Curves_Survival_carrier_age.pdf",
    RES / "17_snp_clinical_associations" / "supplementary_figures" / "17_KM_Curves_Survival_carrier_age_bmi.pdf",
    RES / "17_snp_clinical_associations" / "supplementary_figures" / "17_KM_Curves_Survival_carrier_bmi.pdf",
    RES / "17_snp_clinical_associations" / "supplementary_figures" / "17_KM_Curves_Survival_carrier_only.pdf",
]

STAGE18_PDFS = [
    RES / "18_haplotype_clinical_associations" / "supplementary_figures" / "GSDMB_Haplotype_KM_Curves_Endometrial.pdf",
    RES / "18_haplotype_clinical_associations" / "supplementary_figures" / "GSDMB_Haplotype_KM_Curves_All_Cohorts.pdf",
]

STAGE17_PRIORITY = {
    "carrier_age_bmi": 0,
    "carrier_age": 1,
    "carrier_bmi": 2,
    "carrier_only": 3,
    "genotype_age_bmi": 4,
    "genotype_age": 5,
    "genotype_bmi": 6,
    "genotype_only": 7,
}

STAGE18_PRIORITY = {
    "GSDMB_Haplotype_KM_Curves_Endometrial.pdf": 0,
    "GSDMB_Haplotype_KM_Curves_All_Cohorts.pdf": 1,
}


def s(value) -> str:
    return "" if value is None else str(value).strip()


def normalise_cohort(text: str) -> str:
    text = s(text).lower()
    if "endomet" in text:
        return "endometrium"
    if "breast" in text:
        return "breast"
    return text


def normalise_endpoint(text: str) -> str:
    text = s(text).lower()
    if "progression-free survival" in text or re.search(r"\bpfs\b", text):
        return "pfs"
    if "overall survival" in text or re.search(r"\bos\b", text):
        return "os"
    return text


def parse_float(pattern: str, text: str) -> float:
    match = re.search(pattern, text, flags=re.I)
    if not match:
        return math.inf
    value = match.group(1).strip().lower()
    value = value.replace("***", "").replace("**", "").replace("*", "").replace("(ns)", "").strip()
    try:
        return float(value)
    except Exception:
        return math.inf


def stage17_model_slug(pdf_name: str) -> str:
    base = pdf_name.removeprefix("17_KM_Curves_Survival").removesuffix(".pdf")
    mapping = {
        "": "genotype_only",
        "_age": "genotype_age",
        "_age_bmi": "genotype_age_bmi",
        "_bmi": "genotype_bmi",
        "_carrier_age": "carrier_age",
        "_carrier_age_bmi": "carrier_age_bmi",
        "_carrier_bmi": "carrier_bmi",
        "_carrier_only": "carrier_only",
        "_genotype_age": "genotype_age",
        "_genotype_age_bmi": "genotype_age_bmi",
        "_genotype_bmi": "genotype_bmi",
    }
    if base not in mapping:
        raise ValueError(f"Unknown stage-17 PDF name: {pdf_name}")
    return mapping[base]


def extract_page_rows(pdf_path: Path) -> list[dict]:
    reader = PdfReader(str(pdf_path))
    rows = []
    for page_index, page in enumerate(reader.pages):
        text = page.extract_text() or ""
        text_lower = text.lower()
        rsid_match = re.search(r"(rs\d+)", text_lower)
        hap_match = re.search(r"according to ([a-z0-9_]+) genotype grouping", text_lower)
        p_logrank = parse_float(r"log-rank p=([0-9.eE+-]+(?:\s*\*+)?)", text)
        p_cox = parse_float(r"cox additive p=([0-9.eE+-]+(?:\s*\*+)?)", text)
        q_value = parse_float(r"\bq=([0-9.eE+-]+)", text)
        cohort = normalise_cohort(text_lower)
        endpoint = normalise_endpoint(text_lower)

        if pdf_path in STAGE17_PDFS:
            source_stage = "stage17_snp"
            model_slug = stage17_model_slug(pdf_path.name)
            family = "carrier" if model_slug.startswith("carrier") else "genotype"
            entity_id = rsid_match.group(1) if rsid_match else ""
            label = f"Stage 17 SNP | {entity_id} | {cohort} | {endpoint} | {model_slug}"
            shortlist_key = f"{source_stage}|{entity_id}|{cohort}|{endpoint}|{family}"
            shortlist_priority = STAGE17_PRIORITY[model_slug]
        else:
            source_stage = "stage18_haplotype"
            model_slug = pdf_path.name
            family = "haplotype"
            entity_id = hap_match.group(1) if hap_match else ""
            label = f"Stage 18 haplotype | {entity_id} | {cohort} | {endpoint}"
            shortlist_key = f"{source_stage}|{entity_id}|{cohort}|{endpoint}"
            shortlist_priority = STAGE18_PRIORITY.get(pdf_path.name, 999)

        significant = (p_logrank < 0.05) or (p_cox < 0.05) or (q_value < 0.10)
        rows.append(
            {
                "stage": source_stage,
                "pdf_path": pdf_path,
                "pdf_name": pdf_path.name,
                "page_index": page_index,
                "page_number": page_index + 1,
                "text": text,
                "label": label,
                "entity_id": entity_id,
                "cohort": cohort,
                "endpoint": endpoint,
                "model_slug": model_slug,
                "family": family,
                "p_logrank": p_logrank,
                "p_cox": p_cox,
                "q_value": q_value,
                "significant": significant,
                "shortlist_key": shortlist_key,
                "shortlist_priority": shortlist_priority,
            }
        )
    return rows


def write_bundle(pages_df: pd.DataFrame, output_pdf: Path) -> None:
    writer = PdfWriter()
    cache: dict[Path, PdfReader] = {}
    for row in pages_df.itertuples(index=False):
        reader = cache.setdefault(row.pdf_path, PdfReader(str(row.pdf_path)))
        writer.add_page(reader.pages[row.page_index])
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    with output_pdf.open("wb") as handle:
        writer.write(handle)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    active_pdfs = [path for path in STAGE17_PDFS + STAGE18_PDFS if path.exists()]
    rows: list[dict] = []
    for pdf_path in active_pdfs:
        rows.extend(extract_page_rows(pdf_path))

    pages = pd.DataFrame(rows)
    sig = pages[pages["significant"]].copy()
    if sig.empty:
        raise RuntimeError("No significant survival/KM pages were detected in the active PDFs.")

    sig = sig.sort_values(["stage", "cohort", "endpoint", "entity_id", "shortlist_priority", "p_logrank", "p_cox", "pdf_name", "page_number"]).reset_index(drop=True)

    shortlist = (
        sig.sort_values(["shortlist_key", "shortlist_priority", "q_value", "p_logrank", "p_cox", "pdf_name", "page_number"])
        .drop_duplicates(["shortlist_key"])
        .reset_index(drop=True)
    )

    all_pdf = OUT_DIR / "Significant_Survival_Curves_All.pdf"
    shortlist_pdf = OUT_DIR / "Significant_Survival_Curves_Thesis_Shortlist.pdf"
    index_tsv = OUT_DIR / "Significant_Survival_Curves_Index.tsv"

    write_bundle(sig, all_pdf)
    write_bundle(shortlist, shortlist_pdf)

    index_df = pd.concat(
        [
            sig.assign(bundle="all"),
            shortlist.assign(bundle="shortlist"),
        ],
        ignore_index=True,
    )[
        [
            "bundle",
            "stage",
            "label",
            "entity_id",
            "cohort",
            "endpoint",
            "model_slug",
            "pdf_name",
            "page_number",
            "p_logrank",
            "p_cox",
            "q_value",
        ]
    ]
    index_df.to_csv(index_tsv, sep="\t", index=False)

    print("Created:")
    print(f"  {all_pdf}")
    print(f"  {shortlist_pdf}")
    print(f"  {index_tsv}")
    print("")
    print(f"All-pages bundle pages   : {len(sig)}")
    print(f"Shortlist bundle pages   : {len(shortlist)}")
    print(f"Stage-17 pages included  : {int(sig['stage'].eq('stage17_snp').sum())}")
    print(f"Stage-18 pages included  : {int(sig['stage'].eq('stage18_haplotype').sum())}")


if __name__ == "__main__":
    main()
