"""
Shared visual style definitions for the GSDMB thesis pipeline.

This module centralises semantic colour choices so that the same biological
concepts are rendered consistently across descriptive, comparative, static, and
interactive outputs.

The palettes below prioritise colour-blind-friendly combinations based on the
Okabe-Ito family of colours, supplemented with neutral greys for non-carrier or
non-significant states.
"""

from __future__ import annotations

from typing import Dict, List


OKABE_ITO: Dict[str, str] = {
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "black": "#000000",
    "grey": "#7F7F7F",
    "light_grey": "#D9D9D9",
}

TUMOUR_CONTROL_COLORS: Dict[str, str] = {
    "Tumour": OKABE_ITO["vermillion"],
    "Control": OKABE_ITO["blue"],
    "Healthy": OKABE_ITO["blue"],
    "All tumours": OKABE_ITO["vermillion"],
    "All controls": OKABE_ITO["blue"],
    "Pooled control": OKABE_ITO["blue"],
    "Breast tumour": OKABE_ITO["vermillion"],
    "Endometrium tumour": OKABE_ITO["vermillion"],
}

COHORT_COLORS: Dict[str, str] = {
    "Breast": OKABE_ITO["reddish_purple"],
    "Endometrium": OKABE_ITO["bluish_green"],
    "Endometrial": OKABE_ITO["bluish_green"],
    "Breast_Tumour": OKABE_ITO["reddish_purple"],
    "Breast_Healthy": "#F4CAE4",
    "Endometrium_Tumour": OKABE_ITO["bluish_green"],
    "Endometrial_Tumour": OKABE_ITO["bluish_green"],
    "Endometrium_Healthy": "#B3E2CD",
    "Endometrial_Healthy": "#B3E2CD",
}

IMPACT_COLORS: Dict[str, str] = {
    "HIGH": OKABE_ITO["vermillion"],
    "MODERATE": OKABE_ITO["orange"],
    "LOW": OKABE_ITO["bluish_green"],
    "MODIFIER": OKABE_ITO["grey"],
}

GENOTYPE_COLORS: Dict[str, str] = {
    "WT": "#BDBDBD",
    "Het": OKABE_ITO["orange"],
    "Hom": OKABE_ITO["vermillion"],
}

HAPLOTYPE_STATUS_COLORS: Dict[str, str] = {
    "Wild-type": "#BDBDBD",
    "Heterozygous": OKABE_ITO["orange"],
    "Homozygous": OKABE_ITO["vermillion"],
}

QUALITATIVE_COLORBLIND_SEQUENCE: List[str] = [
    OKABE_ITO["blue"],
    OKABE_ITO["vermillion"],
    OKABE_ITO["bluish_green"],
    OKABE_ITO["reddish_purple"],
    OKABE_ITO["orange"],
    OKABE_ITO["sky_blue"],
    OKABE_ITO["yellow"],
    OKABE_ITO["black"],
    OKABE_ITO["grey"],
]

SEQUENTIAL_COLORBLIND_SCALE = "Cividis"

DASHBOARD_THEME: Dict[str, str] = {
    "header_start": OKABE_ITO["blue"],
    "header_end": OKABE_ITO["bluish_green"],
    "section_start": OKABE_ITO["sky_blue"],
    "section_end": OKABE_ITO["bluish_green"],
    "section_accent": OKABE_ITO["orange"],
    "card_number": OKABE_ITO["blue"],
    "table_header": OKABE_ITO["blue"],
    "table_alt": "#EEF4F8",
}

DESCRIPTIVE_TAG = "Descriptive only"
COMPARATIVE_TAG = ""


def tagged_title(title: str, tag: str) -> str:
    """Append a bracketed tag only when the tag is non-empty."""
    clean_tag = str(tag).strip()
    return title if not clean_tag else f"{title} [{clean_tag}]"


def cohort_color(label: str, default: str = "#666666") -> str:
    """Return a stable cohort colour, with substring fallback for renamed cohorts."""
    if label in COHORT_COLORS:
        return COHORT_COLORS[label]
    if "Breast" in str(label):
        return COHORT_COLORS["Breast"]
    if "Endometri" in str(label):
        return COHORT_COLORS["Endometrial"]
    return default


def arm_color(label: str, default: str = "#666666") -> str:
    """Return a stable tumour/control colour, using simple semantic fallbacks."""
    if label in TUMOUR_CONTROL_COLORS:
        return TUMOUR_CONTROL_COLORS[label]
    text = str(label).lower()
    # Keep the American-spelt fallback for compatibility with raw source labels.
    if "tumour" in text or "tumor" in text:
        return TUMOUR_CONTROL_COLORS["Tumour"]
    if "control" in text or "healthy" in text or "normal" in text:
        return TUMOUR_CONTROL_COLORS["Control"]
    return default

def cohort_colour(label: str, default: str = "#666666") -> str:
    """British-English alias for cohort_color()."""
    return cohort_color(label, default=default)


def arm_colour(label: str, default: str = "#666666") -> str:
    """British-English alias for arm_color()."""
    return arm_color(label, default=default)
