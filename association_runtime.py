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
        "out_dir": paths["results_dir"],
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
        "out_dir": paths["results_dir"],
        "haplotype_results": paths["haplotype_results"],
        "min_hap_freq": thresholds["min_hap_freq"],
        "min_carriers": thresholds["min_carriers"],
        "min_events_logistic": thresholds["min_events_logistic"],
        "fdr_threshold": thresholds["fdr_threshold"],
        "grouping": get_grouping(),
    }
