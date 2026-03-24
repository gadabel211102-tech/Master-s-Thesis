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
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
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

    # --- 7. MULTIPLE TESTING CORRECTION ---
    res_df = pd.DataFrame(all_results)
    final_dfs = []

    for group in res_df["Analysis_Group"].unique():
        sub = res_df[res_df["Analysis_Group"] == group].copy()
        _, sub["FDR_P_Value"], _, _ = multipletests(sub["P_Value"], method="fdr_bh")
        final_dfs.append(sub)

    final_res = pd.concat(final_dfs, ignore_index=True).sort_values(
        ["Analysis_Group", "P_Value"]
    )

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
    final_res["Raw_Significant"] = final_res["P_Value"] < 0.05
    final_res["FDR_Significant"] = final_res["FDR_P_Value"] < 0.05
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
        "Analysis_Group", "SNP_ID", "Coverage_Risk_Flag", "Coverage_Risk_Amplicon",
        "Coverage_Risk_Region", "Coverage_Risk_Note", "Symbol", "Consequence", "Impact",
        "gnomAD_NFE_AF", "gnomAD_NFE_Source",
        "Tumour_Carriers", "Control_Carriers",
        "Tumour_Total", "Control_Total",
        "Tumour_Freq_%", "Control_Freq_%",
        "Odds_Ratio", "P_Value", "FDR_P_Value",
        "Raw_Significant", "FDR_Significant",
        "Calculation_Method"
    ]
    other_cols = [c for c in final_res.columns if c not in preferred_cols]
    final_res = final_res[preferred_cols + other_cols]

    # --- 9. SAVE EXCEL RESULTS ---
    with pd.ExcelWriter(output_xlsx) as writer:
        for group in cohort_list:
            if group in final_res["Analysis_Group"].unique():
                final_res[final_res["Analysis_Group"] == group].to_excel(
                    writer, sheet_name=group, index=False
                )

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

        for _, row in labelled.iterrows():
            ax.annotate(
                row[label_col],
                xy=(row["Odds_Ratio_Plot"], row[y_col]),
                xytext=(8, 8),
                textcoords="offset points",
                fontsize=8,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="yellow", alpha=0.7),
                arrowprops=dict(arrowstyle="->", color="black", lw=0.8)
            )

        ax.set_ylabel(y_label, fontweight="bold")
        ax.set_xlabel("Odds Ratio (log scale)", fontweight="bold")
        ax.grid(True, linestyle=":", alpha=0.4)
        ax.legend(loc="lower left", fontsize=8, frameon=True)

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
            height=5.8,
            aspect=1.15,
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
            g._legend.set_bbox_to_anchor((1.02, 0.5))
            g._legend._loc = 6

        g.fig.subplots_adjust(top=0.82, right=0.86, wspace=0.08)
        g.fig.suptitle(tagged_title(f"SNP Enrichment Volcano: {title}", COMPARATIVE_TAG), fontsize=16, fontweight="bold", y=0.97)

        g.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(g.fig)

    # --- 11. PLOT 1: RAW P-VALUE VOLCANO ---
    make_volcano_plot(
        data=final_res,
        y_col="neglog10_raw_p",
        label_col="Label_RAW",
        threshold_p=0.05,
        y_label="-log10(raw P-value)",
        title="Raw p-values",
        output_path=output_plot_raw
    )

    # --- 12. PLOT 2: FDR-ADJUSTED VOLCANO ---
    make_volcano_plot(
        data=final_res,
        y_col="neglog10_fdr_p",
        label_col="Label_FDR",
        threshold_p=0.05,
        y_label="-log10(FDR-adjusted P-value)",
        title="FDR-adjusted",
        output_path=output_plot_fdr
    )

    print("\n>>> SUCCESS: All SNPs processed.")
    print(f">>> Results saved to: {output_xlsx}")
    print(">>> Combined NFE AF logic used:")
    print("    - gnomADe_NFE_AF if available")
    print("    - otherwise gnomADg_NFE_AF")
    print(">>> Volcano plots saved:")
    print(f"    - RAW: {output_plot_raw}")
    print(f"    - FDR: {output_plot_fdr}")


if __name__ == "__main__":
    run_snp_association_analysis()
