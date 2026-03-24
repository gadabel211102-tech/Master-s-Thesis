"""Shared runtime defaults for the association scripts.

This module centralises the file paths, thresholds and grouping choices used by
scripts 17 and 18. Keeping those defaults in one place makes the analytical
assumptions easier to audit and reduces the risk that the SNP and haplotype
association pipelines drift out of sync.
"""

from __future__ import annotations

from pipeline_utils import get_grouping, get_paths, get_thresholds


def script17_defaults() -> dict[str, object]:
    """Return the canonical configuration bundle for SNP-clinical analysis."""
    paths = get_paths()
    thresholds = get_thresholds()
    manifests_dir = paths["manifests_dir"]
    return {
        "gsdmb": paths["annotated_report"],
        "master": paths["harmonised_master"],
        "out_dir": paths["snp_clinical_dir"],
        "manifests": {
            "endometrium-tumour": manifests_dir / "endometrium-tumour-pass_manifest.txt",
            "endometrium-normal": manifests_dir / "endometrium-normal-pass_manifest.txt",
            "breast-tumour": manifests_dir / "breast-tumour-pass_manifest.txt",
            "breast-normal": manifests_dir / "breast-normal-pass_manifest.txt",
        },
        "fdr_threshold": thresholds["fdr_threshold"],
        "min_carriers": thresholds["min_carriers"],
        "min_nfe_af": thresholds["min_nfe_af"],
        "grouping": get_grouping(),
    }


def script18_defaults() -> dict[str, object]:
    """Return the canonical configuration bundle for haplotype analysis."""
    paths = get_paths()
    thresholds = get_thresholds()
    return {
        "phased": paths["haplotype_phased"],
        "master": paths["harmonised_master"],
        "annot": paths["annotated_report"],
        "out_dir": paths["haplotype_clinical_dir"],
        "haplotype_results": paths["haplotype_results"],
        "min_hap_freq": thresholds["min_hap_freq"],
        "min_carriers": thresholds["min_carriers"],
        "min_events_logistic": thresholds["min_events_logistic"],
        "fdr_threshold": thresholds["fdr_threshold"],
        "grouping": get_grouping(),
    }


def script19_defaults() -> dict[str, object]:
    """Return the canonical configuration bundle for 1000 Genomes haplotype comparison."""
    paths = get_paths()
    out_dir = paths.get("thousand_genomes_results", paths["results_dir"] / "19_1000g_haplotype_comparison")
    return {
        "annot": paths["annotated_report"],
        "haplotype_results": paths["haplotype_results"],
        "out_dir": out_dir,
        "populations": ["EUR", "ALL", "IBS"],
    }



def script19b_defaults() -> dict[str, object]:
    """Return the canonical configuration bundle for the standalone RNA QC stage."""
    paths = get_paths()
    return {
        "expr_xlsx": paths["objective2_expression_workbook"],
        "master": paths["harmonised_master"],
        "rna_bed": paths["objective2_rna_bed"],
        "rna_qc_manifest": paths["objective2_rna_qc_manifest"],
        "out_dir": paths["objective2_rna_qc_dir"],
        "rna_qc_manifest": paths["objective2_rna_qc_manifest"],
    }


def script20_defaults() -> dict[str, object]:
    """Return the canonical configuration bundle for objective-2 isoform analysis."""
    paths = get_paths()
    thresholds = get_thresholds()
    return {
        "expr_xlsx": paths["objective2_expression_workbook"],
        "master": paths["harmonised_master"],
        "variant_workbook": paths["annotated_report"],
        "haplotype_input": paths["haplotype_phased"],
        "out_dir": paths["isoform_expression_dir"],
        "rna_bed": paths["objective2_rna_bed"],
        "rna_qc_manifest": paths["objective2_rna_qc_manifest"],
        "min_nfe_af": thresholds["min_nfe_af"],
        "min_hap_freq": thresholds["min_hap_freq"],
        "min_carriers": thresholds["min_carriers"],
        "fdr_threshold": thresholds["fdr_threshold"],
        "grouping": get_grouping(),
    }


