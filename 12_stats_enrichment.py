"""
Script 12: SNP Enrichment Testing
=================================

Purpose
-------
Test whether common SNPs are enriched in tumour samples relative to the control
arm, both globally and within each tumour cohort.

Statistical approach
--------------------
For each SNP, the script constructs a 2x2 carrier table and applies Fisher's
exact test. Odds ratios are reported, with a Haldane-Anscombe correction when
zero cells would otherwise make the estimate unstable. False-discovery-rate
correction is then applied within each analysis group.

Control definition
------------------
Cohort-specific tumour analyses use the pooled control group agreed for this
project: breast normal plus endometrium normal.
"""

import pandas as pd
import textwrap
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
import seaborn as sns
from statsmodels.stats.multitest import multipletests

from figure_style import COMPARATIVE_TAG, IMPACT_COLORS, arm_color, tagged_title
from pipeline_utils import (
    attach_amplicon_warning_columns,
    build_amplicon_warning_lookup,
    build_variant_id_series,
    combine_gnomad_nfe,
    ensure_directory,
    find_col,
    get_paths,
    get_thresholds,
    pooled_control_frame,
    standardize_cohort_labels,
    standardize_tissue_labels,
    load_common_snp_ld_reference,
)
from pipeline_validation import print_validation_summary, validate_file_exists, validate_nonempty, validate_percentage_columns, validate_required_columns, validate_tissue_values

# --- 1. PATH CONFIGURATION ---
PATHS = get_paths()
THRESHOLDS = get_thresholds()
output_dir = ensure_directory(PATHS["snp_enrichment_dir"])
input_file = str(PATHS["annotated_report"])
output_xlsx = str(output_dir / "GSDMB_SNP_Enrichment_Results.xlsx")

# Separate volcano outputs
output_plot_raw = str(output_dir / "GSDMB_SNP_Enrichment_Volcano_RAW.png")
output_plot_fdr = str(output_dir / "GSDMB_SNP_Enrichment_Volcano_FDR.png")
output_ld_summary_xlsx = str(output_dir / "GSDMB_SNP_Enrichment_LD_Block_Summary.xlsx")
LD_R2_THRESHOLD = 0.80
RAW_P_THRESHOLD = 0.05
FDR_THRESHOLD = float(THRESHOLDS.get("fdr_threshold", 0.10))


def _short_label(value, max_len=28):
    text = " ".join(str(value).split())
    return text if len(text) <= max_len else text[: max_len - 1].rstrip() + "?"


def _annotate_ranked_hits(ax, labelled, label_col, x_col, y_col, max_labels=5):
    if labelled.empty:
        return
    ranked = labelled.sort_values([y_col, x_col], ascending=[False, True]).head(max_labels)
    y_min, y_max = ax.get_ylim()
    pad = (y_max - y_min) * 0.08 if y_max != y_min else 0.5
    if len(ranked) == 1:
        y_positions = [min(y_max - pad, ranked.iloc[0][y_col] + pad * 0.25)]
    else:
        upper = max(y_max - pad, y_min + pad)
        lower = max(y_min + pad, y_min + (y_max - y_min) * 0.18)
        y_positions = np.linspace(upper, lower, len(ranked))

    text_transform = blended_transform_factory(ax.transAxes, ax.transData)
    x_frac = 0.74

    for (_, row), y_text in zip(ranked.iterrows(), y_positions):
        ax.annotate(
            _short_label(row[label_col], max_len=22),
            xy=(row[x_col], row[y_col]),
            xycoords="data",
            xytext=(x_frac, y_text),
            textcoords=text_transform,
            ha="left",
            va="center",
            fontsize=9.2,
            fontweight="bold",
            clip_on=False,
            bbox=dict(boxstyle="round,pad=0.26", fc="white", ec="#9e9e9e", alpha=0.96, linewidth=0.8),
            arrowprops=dict(arrowstyle="-", color="#777777", lw=0.8, shrinkA=0, shrinkB=4),
            zorder=5,
        )



def _singleton_block_id(snp_id):
    return f"Singleton::{snp_id}"


def _initialise_ld_block_columns(res_df, ld_reference):
    block_df = ld_reference["block_membership"].copy()
    if not block_df.empty:
        block_df["LD_Block_ID"] = block_df["Block"].astype(str)
        block_map = dict(zip(block_df["rsID"], block_df["LD_Block_ID"]))
        block_sizes = block_df.groupby("LD_Block_ID")["rsID"].size().to_dict()
        block_members = (
            block_df.groupby("LD_Block_ID")["rsID"]
            .apply(lambda series: ", ".join(series.astype(str)))
            .to_dict()
        )
    else:
        block_map = {}
        block_sizes = {}
        block_members = {}

    out = res_df.copy()
    out["LD_Threshold"] = ld_reference["threshold"]
    out["Threshold_Label"] = ld_reference["threshold_label"]
    out["LD_Block_ID"] = out["SNP_ID"].map(block_map)
    out["LD_Block_Assignment"] = np.where(out["LD_Block_ID"].notna(), "LD_Block", "Singleton")
    out.loc[out["LD_Block_ID"].isna(), "LD_Block_ID"] = out.loc[out["LD_Block_ID"].isna(), "SNP_ID"].map(_singleton_block_id)
    out["LD_Block_Size"] = out["LD_Block_ID"].map(block_sizes).fillna(1).astype(int)
    out["LD_Block_Members"] = out["LD_Block_ID"].map(block_members).fillna(out["SNP_ID"])
    return out, block_df


def _apply_block_fdr(res_df, group_col, p_col, block_col, primary_fdr_col, legacy_fdr_col):
    final_parts = []
    for group_value, group_df in res_df.groupby(group_col, sort=False):
        group_df = group_df.copy()
        _, group_df[legacy_fdr_col], _, _ = multipletests(group_df[p_col], method="fdr_bh")
        block_parts = []
        for _, block_df in group_df.groupby(block_col, sort=False):
            block_df = block_df.copy()
            if len(block_df) == 1:
                block_df[primary_fdr_col] = block_df[p_col]
            else:
                _, block_df[primary_fdr_col], _, _ = multipletests(block_df[p_col], method="fdr_bh")
            block_df["LD_Block_Test_Count"] = int(len(block_df))
            block_parts.append(block_df)
        group_df = pd.concat(block_parts, ignore_index=True)
        threshold_label = group_df["Threshold_Label"].iloc[0] if "Threshold_Label" in group_df.columns else "LD"
        group_df["FDR_Scope"] = f"Within_{group_value}_{threshold_label}_blocks"
        final_parts.append(group_df)
    return pd.concat(final_parts, ignore_index=True)


def _build_ld_summary_tables(final_res, ld_reference, block_membership_df):
    global_df = final_res[final_res["Analysis_Group"] == "Global"].copy()
    if global_df.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    global_df["Legacy_Global_FDR_Significant"] = global_df["FDR_P_Value_Global_Legacy"] < FDR_THRESHOLD
    global_df["Block_FDR_Significant"] = global_df["FDR_P_Value"] < FDR_THRESHOLD

    block_context = (
        global_df.groupby("LD_Block_ID", as_index=False)
        .agg(
            Block_Test_Count=("LD_Block_Test_Count", "first"),
            Block_Min_P_Value=("P_Value", "min"),
            Block_Min_FDR_P_Value=("FDR_P_Value", "min"),
            Block_Min_Global_Legacy_FDR=("FDR_P_Value_Global_Legacy", "min"),
            Block_FDR_Significant_SNPs=("Block_FDR_Significant", "sum"),
            Legacy_Global_FDR_Significant_SNPs=("Legacy_Global_FDR_Significant", "sum"),
            Block_Tumour_Freq_Min_pct=("Tumour_Freq_%", "min"),
            Block_Tumour_Freq_Max_pct=("Tumour_Freq_%", "max"),
            Block_Control_Freq_Min_pct=("Control_Freq_%", "min"),
            Block_Control_Freq_Max_pct=("Control_Freq_%", "max"),
        )
    )

    pairwise = ld_reference["ld_long"].copy()
    sig_snps = sorted(global_df.loc[global_df["Legacy_Global_FDR_Significant"], "SNP_ID"].astype(str).unique())
    if not sig_snps:
        sig_snps = sorted(global_df.loc[global_df["Block_FDR_Significant"], "SNP_ID"].astype(str).unique())
    pairwise_sig = pairwise[
        pairwise["SNP1"].isin(sig_snps)
        & pairwise["SNP2"].isin(sig_snps)
        & (pairwise["SNP1"] < pairwise["SNP2"])
    ].copy()
    if not pairwise_sig.empty:
        block_lookup = dict(zip(global_df["SNP_ID"], global_df["LD_Block_ID"]))
        pairwise_sig["LD_Block_ID_1"] = pairwise_sig["SNP1"].map(block_lookup)
        pairwise_sig["LD_Block_ID_2"] = pairwise_sig["SNP2"].map(block_lookup)
        pairwise_sig["Same_LD_Block"] = pairwise_sig["LD_Block_ID_1"] == pairwise_sig["LD_Block_ID_2"]
        pairwise_sig["Pass_r2_0_80"] = pairwise_sig["R2"] >= LD_R2_THRESHOLD

    pairwise_strings = {}
    if not pairwise_sig.empty:
        for snp in sig_snps:
            rows = pairwise_sig[(pairwise_sig["SNP1"] == snp) | (pairwise_sig["SNP2"] == snp)].copy()
            if rows.empty:
                pairwise_strings[snp] = ""
                continue
            bits = []
            for _, row in rows.iterrows():
                other = row["SNP2"] if row["SNP1"] == snp else row["SNP1"]
                same = "same block" if row["Same_LD_Block"] else "different block"
                bits.append(f"{other} (r2={row['R2']:.3f}, {same})")
            pairwise_strings[snp] = "; ".join(bits)

    global_summary = global_df.merge(block_context, on="LD_Block_ID", how="left")
    global_summary["Pairwise_LD_to_Global_Significant_SNPs"] = global_summary["SNP_ID"].map(pairwise_strings).fillna("")
    ordered_cols = [
        "SNP_ID", "LD_Block_ID", "LD_Block_Assignment", "LD_Block_Size", "LD_Block_Members",
        "P_Value", "FDR_P_Value", "FDR_P_Value_Global_Legacy",
        "Raw_Significant", "FDR_Significant", "Legacy_Global_FDR_Significant",
        "Tumour_Freq_%", "Control_Freq_%", "Odds_Ratio",
        "Block_Test_Count", "Block_Min_P_Value", "Block_Min_FDR_P_Value", "Block_Min_Global_Legacy_FDR",
        "Block_FDR_Significant_SNPs", "Legacy_Global_FDR_Significant_SNPs",
        "Block_Tumour_Freq_Min_pct", "Block_Tumour_Freq_Max_pct",
        "Block_Control_Freq_Min_pct", "Block_Control_Freq_Max_pct",
        "Pairwise_LD_to_Global_Significant_SNPs",
    ]
    global_summary = global_summary[ordered_cols + [c for c in global_summary.columns if c not in ordered_cols]]
    return global_summary, pairwise_sig, block_membership_df

def run_snp_association_analysis():
    """Run tumour-versus-control SNP enrichment tests and export the results."""
    print(">>> Starting Script 12: SNP statistical enrichment with separate raw and FDR volcano plots...")

    validate_file_exists(input_file, "Script 12 input")

    # --- 2. DATA LOADING & CLEANING ---
    df = pd.read_excel(input_file, sheet_name="Biological_Annotations")

    sym_c = find_col(df, "SYMBOL")
    imp_c = find_col(df, "IMPACT")
    con_c = find_col(df, "Consequence")
    tis_c = find_col(df, "TISSUE")
    sam_c = find_col(df, "SAMPLE")
    var_c = find_col(df, "Existing_variation")
    hgv_c = find_col(df, "HGVSp")
    coh_c = find_col(df, "COHORT")
    exome_nfe_c = find_col(df, "gnomADe_NFE_AF")
    genome_nfe_c = find_col(df, "gnomADg_NFE_AF")

    validate_required_columns(df, [sym_c, imp_c, con_c, tis_c, sam_c, var_c, hgv_c, coh_c, exome_nfe_c, genome_nfe_c], "Script 12 Biological_Annotations")

    # --- 3. BUILD COMBINED NFE AF COLUMN ---
    df = combine_gnomad_nfe(df, exome_nfe_c, genome_nfe_c)

    # --- 4. SNP FILTER: established common SNPs ---
    df = df[df["gnomAD_NFE_AF_combined"] > THRESHOLDS["min_nfe_af"]].copy()
    validate_nonempty(df, "Script 12 common SNP subset")

    df[tis_c] = standardize_tissue_labels(df[tis_c])
    df[coh_c] = standardize_cohort_labels(df[coh_c])

    # Keep only expected tissue groups
    df = df[df[tis_c].isin(["Tumour", "Healthy"])].copy()
    validate_tissue_values(df, tis_c, "Script 12 filtered tissues")

    # --- 5. ROBUST SNP ID CREATION ---
    df["Variant_ID"] = build_variant_id_series(df[var_c], df[sym_c], df[hgv_c])
    warning_lookup = build_amplicon_warning_lookup(df, variant_col="Variant_ID", chrom_col="CHROM", pos_col="POS")
    print_validation_summary(df, sam_c, "Script 12 filtered SNPs", [coh_c, tis_c])

    cohort_list = ["Global", "Breast", "Endometrium"]
    cohort_filters = {"Breast": "Breast", "Endometrium": "Endometri"}
    all_results = []

    # --- 6. FISHER'S EXACT TEST WITH HALDANE-ANSCOMBE CORRECTION FOR OR ---
    for cohort_name in cohort_list:
        if cohort_name == "Global":
            cohort_df = pooled_control_frame(df, coh_c, tis_c)
        else:
            cohort_df = pooled_control_frame(df, coh_c, tis_c, cohort_name, cohort_filters[cohort_name])

        if cohort_df.empty:
            continue

        # Sample totals are calculated on unique carriers/non-carriers within
        # the tumour arm and the pooled control arm defined above.
        n_tumour = cohort_df[cohort_df[tis_c] == "Tumour"][sam_c].nunique()
        n_healthy = cohort_df[cohort_df[tis_c] == "Healthy"][sam_c].nunique()

        if n_tumour == 0 or n_healthy == 0:
            continue

        unique_snps = cohort_df["Variant_ID"].dropna().unique()

        for var in unique_snps:
            var_data = cohort_df[cohort_df["Variant_ID"] == var]

            count_tumour = var_data[var_data[tis_c] == "Tumour"][sam_c].nunique()
            count_healthy = var_data[var_data[tis_c] == "Healthy"][sam_c].nunique()

            a = int(count_tumour)
            b = int(max(0, n_tumour - count_tumour))
            c = int(count_healthy)
            d = int(max(0, n_healthy - count_healthy))

            # Raw Fisher test
            _, p_value = stats.fisher_exact([[a, b], [c, d]])

            # Odds ratio, corrected only if any cell is zero
            if a == 0 or b == 0 or c == 0 or d == 0:
                or_val = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
                note = "Corrected (+0.5)"
            else:
                or_val = (a * d) / (b * c) if (b * c) != 0 else np.inf
                note = "Standard"

            all_results.append({
                "Analysis_Group": cohort_name,
                "SNP_ID": var,
                "Symbol": var_data[sym_c].iloc[0],
                "Consequence": var_data[con_c].iloc[0],
                "Impact": var_data[imp_c].iloc[0],
                "gnomAD_NFE_AF": var_data["gnomAD_NFE_AF_combined"].iloc[0],
                "gnomAD_NFE_Source": var_data["gnomAD_NFE_Source"].iloc[0],
                "Tumour_Carriers": count_tumour,
                "Control_Carriers": count_healthy,
                "Tumour_Total": n_tumour,
                "Control_Total": n_healthy,
                "Tumour_Freq_%": (count_tumour / n_tumour) * 100,
                "Control_Freq_%": (count_healthy / n_healthy) * 100,
                "Odds_Ratio": or_val,
                "P_Value": p_value,
                "Calculation_Method": note
            })

    if not all_results:
        print("No SNP association results generated.")
        return

    # --- 7. LD STRUCTURE + MULTIPLE TESTING CORRECTION ---
    res_df = pd.DataFrame(all_results)
    ld_reference = load_common_snp_ld_reference(
        annotated_report=PATHS["annotated_report"],
        phased_genotypes=PATHS["haplotype_phased"],
        min_nfe_af=THRESHOLDS["min_nfe_af"],
        common_snp_whitelist=PATHS["common_snps_dir"] / "GSDMB_Common_SNP_Frequency_Summary.xlsx",
        threshold=LD_R2_THRESHOLD,
    )
    res_df, ld_block_membership = _initialise_ld_block_columns(res_df, ld_reference)
    final_res = _apply_block_fdr(
        res_df,
        group_col="Analysis_Group",
        p_col="P_Value",
        block_col="LD_Block_ID",
        primary_fdr_col="FDR_P_Value",
        legacy_fdr_col="FDR_P_Value_Global_Legacy",
    ).sort_values(["Analysis_Group", "P_Value"])

    # --- 8. DERIVED COLUMNS FOR PLOTTING ---
    final_res["P_Value_Safe"] = final_res["P_Value"].replace(0, 1e-300)
    final_res["FDR_P_Value_Safe"] = final_res["FDR_P_Value"].replace(0, 1e-300)

    final_res["neglog10_raw_p"] = -np.log10(final_res["P_Value_Safe"])
    final_res["neglog10_fdr_p"] = -np.log10(final_res["FDR_P_Value_Safe"])

    # Safer OR for plotting on log scale
    final_res["Odds_Ratio_Plot"] = (
        final_res["Odds_Ratio"]
        .replace(0, 1e-6)
        .replace(np.inf, 1e6)
    )

    # Separate significance flags
    final_res["Raw_Significant"] = final_res["P_Value"] < RAW_P_THRESHOLD
    final_res["FDR_Significant"] = final_res["FDR_P_Value"] < FDR_THRESHOLD
    final_res["FDR_Significant_Global_Legacy"] = final_res["FDR_P_Value_Global_Legacy"] < FDR_THRESHOLD
    final_res["Raw_P_Threshold"] = RAW_P_THRESHOLD
    final_res["FDR_Threshold"] = FDR_THRESHOLD
    validate_percentage_columns(final_res, ["Tumour_Freq_%", "Control_Freq_%"], "Script 12 results")
    final_res = attach_amplicon_warning_columns(
        final_res,
        warning_lookup.rename(columns={"Variant_ID": "SNP_ID"}),
        variant_col="SNP_ID",
    )

    # Separate labels for each plot
    final_res["Label_RAW"] = ""
    final_res["Label_FDR"] = ""

    for group in final_res["Analysis_Group"].unique():
        group_mask = final_res["Analysis_Group"] == group

        raw_sig = final_res[group_mask & final_res["Raw_Significant"]]
        if not raw_sig.empty:
            top_raw = raw_sig.nsmallest(5, "P_Value").index
            final_res.loc[top_raw, "Label_RAW"] = final_res.loc[top_raw, "SNP_ID"]

        fdr_sig = final_res[group_mask & final_res["FDR_Significant"]]
        if not fdr_sig.empty:
            top_fdr = fdr_sig.nsmallest(5, "FDR_P_Value").index
            final_res.loc[top_fdr, "Label_FDR"] = final_res.loc[top_fdr, "SNP_ID"]

    # Order columns more cleanly for Excel
    preferred_cols = [
        "Analysis_Group", "SNP_ID", "LD_Threshold", "Threshold_Label", "LD_Block_ID", "LD_Block_Assignment",
        "LD_Block_Size", "LD_Block_Test_Count", "LD_Block_Members", "FDR_Scope",
        "Coverage_Risk_Flag", "Coverage_Risk_Amplicon", "Coverage_Risk_Region", "Coverage_Risk_Note",
        "Symbol", "Consequence", "Impact", "gnomAD_NFE_AF", "gnomAD_NFE_Source",
        "Tumour_Carriers", "Control_Carriers", "Tumour_Total", "Control_Total",
        "Tumour_Freq_%", "Control_Freq_%", "Odds_Ratio",
        "P_Value", "FDR_P_Value", "FDR_P_Value_Global_Legacy",
        "Raw_Significant", "FDR_Significant", "FDR_Significant_Global_Legacy",
        "Calculation_Method"
    ]
    other_cols = [c for c in final_res.columns if c not in preferred_cols]
    final_res = final_res[preferred_cols + other_cols]

    # --- 9. SAVE EXCEL RESULTS ---
    global_ld_summary, global_pairwise_ld, ld_membership_export = _build_ld_summary_tables(final_res, ld_reference, ld_block_membership)

    with pd.ExcelWriter(output_xlsx) as writer:
        for group in cohort_list:
            if group in final_res["Analysis_Group"].unique():
                final_res[final_res["Analysis_Group"] == group].to_excel(
                    writer, sheet_name=group, index=False
                )
        if not ld_membership_export.empty:
            ld_membership_export.to_excel(writer, sheet_name="LD_Block_Membership", index=False)
        if not global_pairwise_ld.empty:
            global_pairwise_ld.to_excel(writer, sheet_name="Global_Significant_LD", index=False)
        if not global_ld_summary.empty:
            global_ld_summary.to_excel(writer, sheet_name="Global_Block_Context", index=False)

    with pd.ExcelWriter(output_ld_summary_xlsx) as writer:
        if not global_ld_summary.empty:
            global_ld_summary.to_excel(writer, sheet_name="Global_Block_Context", index=False)
        if not global_pairwise_ld.empty:
            global_pairwise_ld.to_excel(writer, sheet_name="Global_Significant_LD", index=False)
        if not ld_membership_export.empty:
            ld_membership_export.to_excel(writer, sheet_name="LD_Block_Membership", index=False)

    # --- 10. SHARED VOLCANO HELPERS ---
    def add_volcano_decorations(ax, group_name, data, y_col, label_col, threshold_p, y_label):
        """Add threshold lines, shading and labels to one volcano axis."""
        y_thresh = -np.log10(threshold_p)

        ax.axhline(
            y_thresh,
            color=arm_color('Tumour'),
            linestyle="--",
            alpha=0.7,
            linewidth=2,
            label=f"threshold = {threshold_p}"
        )

        ax.axvline(
            1,
            color=arm_color('Control'),
            linestyle="--",
            alpha=0.7,
            linewidth=2,
            label="OR = 1"
        )

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        ax.axvspan(xlim[0], 1, alpha=0.05, color=arm_color('Control'), label="Protective (OR<1)")
        ax.axvspan(1, xlim[1], alpha=0.05, color=arm_color('Tumour'), label="Risk (OR>1)")
        ax.axhspan(y_thresh, ylim[1], alpha=0.08, color=arm_color('Tumour'), label="Above threshold")

        group_data = data[data["Analysis_Group"] == group_name]
        labelled = group_data[group_data[label_col] != ""]
        if not group_data.empty:
            n_tumour = int(group_data['Tumour_Total'].iloc[0])
            n_control = int(group_data['Control_Total'].iloc[0])
            ax.set_title(f"{group_name}\nTumour samples: n={n_tumour}; controls: n={n_control}", fontweight='bold')

        _annotate_ranked_hits(ax, labelled, label_col, "Odds_Ratio_Plot", y_col, max_labels=6)

        ax.set_ylabel(y_label, fontweight="bold")
        ax.set_xlabel("Odds Ratio (log scale)", fontweight="bold")
        ax.grid(True, linestyle=":", alpha=0.4)
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        ax.legend(unique.values(), unique.keys(), loc="lower left", fontsize=7.5, frameon=True)

    def make_volcano_plot(data, y_col, label_col, threshold_p, y_label, title, output_path):
        group_order = ["Breast", "Endometrium", "Global"]
        plot_data = data[data["Analysis_Group"].isin(group_order)].copy()
        plot_data["Analysis_Group"] = pd.Categorical(
            plot_data["Analysis_Group"],
            categories=group_order,
            ordered=True
        )
        plot_data = plot_data.sort_values("Analysis_Group")

        g = sns.FacetGrid(
            plot_data,
            col="Analysis_Group",
            hue="Impact",
            col_order=group_order,
            palette=IMPACT_COLORS,
            height=6.4,
            aspect=1.52,
            despine=False,
            sharex=True,
            sharey=True
        )

        g.map_dataframe(
            sns.scatterplot,
            x="Odds_Ratio_Plot",
            y=y_col,
            s=120,
            edgecolor="black",
            alpha=0.8
        )

        for ax, group_name in zip(g.axes.flat, group_order):
            ax.set_xscale("log")
            add_volcano_decorations(
                ax=ax,
                group_name=group_name,
                data=plot_data,
                y_col=y_col,
                label_col=label_col,
                threshold_p=threshold_p,
                y_label=y_label
            )

        g.add_legend(title="VEP Impact")
        if g._legend is not None:
            g._legend.set_bbox_to_anchor((0.90, 0.56), transform=g.fig.transFigure)
            g._legend._loc = 6
            g._legend.get_title().set_fontsize(11)
            for text in g._legend.texts:
                text.set_fontsize(10)

        g.fig.subplots_adjust(top=0.82, right=0.86, left=0.055, bottom=0.12, wspace=0.10)
        g.fig.suptitle(tagged_title(f"SNP Enrichment Volcano: {title}", COMPARATIVE_TAG), fontsize=16, fontweight="bold", y=0.97)

        g.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(g.fig)

    # --- 11. PLOT 1: RAW P-VALUE VOLCANO ---
    make_volcano_plot(
        data=final_res,
        y_col="neglog10_raw_p",
        label_col="Label_RAW",
        threshold_p=RAW_P_THRESHOLD,
        y_label="-log10(raw P-value)",
        title="Raw p-values",
        output_path=output_plot_raw
    )

    # --- 12. PLOT 2: FDR-ADJUSTED VOLCANO ---
    make_volcano_plot(
        data=final_res,
        y_col="neglog10_fdr_p",
        label_col="Label_FDR",
        threshold_p=FDR_THRESHOLD,
        y_label="-log10(FDR-adjusted P-value)",
        title="FDR-adjusted",
        output_path=output_plot_fdr
    )

    print("\n>>> SUCCESS: All SNPs processed.")
    print(f">>> Results saved to: {output_xlsx}")
    print(f">>> LD-block summary saved to: {output_ld_summary_xlsx}")
    print(">>> Combined NFE AF logic used:")
    print("    - gnomADe_NFE_AF if available")
    print("    - otherwise gnomADg_NFE_AF")
    print(">>> Volcano plots saved:")
    print(f"    - RAW: {output_plot_raw}")
    print(f"    - FDR: {output_plot_fdr}")


if __name__ == "__main__":
    run_snp_association_analysis()