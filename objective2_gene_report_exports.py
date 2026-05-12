from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

try:
    import seaborn as sns
except ImportError:
    sns = None


def z(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    sd = numeric.std(ddof=0)
    return (numeric - numeric.mean()) / sd if pd.notna(sd) and sd != 0 else pd.Series(np.nan, index=numeric.index)


def bh(values) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    out = np.full(arr.shape, np.nan)
    mask = np.isfinite(arr)
    if not mask.any():
        return out
    pvals = arr[mask]
    order = np.argsort(pvals)
    ranked = pvals[order]
    n = len(ranked)
    adj = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, ranked[i] * n / (i + 1))
        adj[i] = min(prev, 1.0)
    back = np.empty(n)
    back[order] = adj
    out[mask] = back
    return out


def add_fdr(df: pd.DataFrame, p_col: str, group_cols: list[str], fdr_threshold: float) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    out["FDR"] = np.nan
    for _, idx in out.groupby(group_cols, dropna=False).groups.items():
        out.loc[idx, "FDR"] = bh(pd.to_numeric(out.loc[idx, p_col], errors="coerce"))
    out["Nominal_Sig"] = pd.to_numeric(out[p_col], errors="coerce") < 0.05
    out["FDR_Sig"] = pd.to_numeric(out["FDR"], errors="coerce") < fdr_threshold
    return out


def safe_gene_slug(gene: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", str(gene).strip())
    slug = re.sub(r"_+", "_", slug).strip("._")
    return slug or "gene"



def pretty_group(group: str) -> str:
    mapping = {
        "Breast_Tumour": "Breast\nTumour",
        "Endometrial_Tumour": "Endometrial\nTumour",
        "Pooled_Control": "Pooled\nControl",
        "Breast_Normal": "Breast\nControl",
        "Endometrial_Normal": "Endometrial\nControl",
    }
    return mapping.get(str(group), str(group).replace("_", "\n"))


def _top_rows(df: pd.DataFrame, p_col: str = "P_Value", top_n: int = 5) -> pd.DataFrame:
    if df is None or df.empty or p_col not in df.columns:
        return pd.DataFrame()
    out = df[df[p_col].notna()].copy()
    if out.empty:
        return pd.DataFrame()
    return out.sort_values(p_col).head(top_n).copy()


def _plot_top_hits(ax, df: pd.DataFrame, title: str, label_fn, colour: str) -> None:
    top = _top_rows(df)
    if top.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, f"No {title.lower()} hits", ha="center", va="center", fontsize=10, color="#666666")
        return
    labels = [label_fn(row) for _, row in top.iterrows()]
    scores = -np.log10(pd.to_numeric(top["P_Value"], errors="coerce").clip(lower=1e-300))
    ypos = np.arange(len(top))[::-1]
    edge = ["#222222" if bool(row.get("FDR_Sig", False)) else "white" for _, row in top.iterrows()]
    ax.barh(ypos, scores, color=colour, alpha=0.82, edgecolor=edge, linewidth=1.3)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("-log10(P value)", fontsize=8.5)
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.grid(axis="x", linestyle=":", alpha=0.25)
    for y, score, (_, row) in zip(ypos, scores, top.iterrows()):
        p = row.get("P_Value", np.nan)
        fdr_flag = " FDR" if bool(row.get("FDR_Sig", False)) else ""
        ax.text(score + 0.05, y, f"p={p:.3g}{fdr_flag}", va="center", fontsize=7.5, color="#444444")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


def make_panel_gene_figure(gene: str, merged: pd.DataFrame, summary_df: pd.DataFrame, snp_df: pd.DataFrame, hap_df: pd.DataFrame, clin_df: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.5), constrained_layout=True)
    ax_expr, ax_snp, ax_hap, ax_clin = axes.flat

    expr = merged[["analysis_group", gene]].copy()
    expr[gene] = gene_series(merged, gene)
    expr = expr[expr["analysis_group"].isin(["Breast_Tumour", "Endometrial_Tumour", "Breast_Normal", "Endometrial_Normal"])].copy()
    expr["plot_group"] = expr["analysis_group"].replace({"Breast_Normal": "Pooled_Control", "Endometrial_Normal": "Pooled_Control"})
    expr = expr.dropna(subset=[gene])
    order = [g for g in ["Breast_Tumour", "Endometrial_Tumour", "Pooled_Control"] if g in set(expr["plot_group"])]
    if not expr.empty and order:
        if sns is not None:
            sns.boxplot(data=expr, x="plot_group", y=gene, order=order, ax=ax_expr, color="#d9e8f5", fliersize=0, width=0.62)
            sns.stripplot(data=expr, x="plot_group", y=gene, order=order, ax=ax_expr, color="#2f5d80", size=2.8, alpha=0.6, jitter=0.22)
        else:
            groups = [pd.to_numeric(expr.loc[expr["plot_group"] == g, gene], errors="coerce").dropna().values for g in order]
            ax_expr.boxplot(groups, labels=order, patch_artist=True, boxprops=dict(facecolor="#d9e8f5"))
        ax_expr.set_xticklabels([pretty_group(g) for g in order], fontsize=8.5)
        ax_expr.set_title("Expression by tumour/control context", fontsize=10, fontweight="bold")
        ax_expr.set_xlabel("")
        ax_expr.set_ylabel(gene, fontsize=9)
        ax_expr.grid(axis="y", linestyle=":", alpha=0.25)
        if not summary_df.empty:
            parts = []
            for grp in ["Breast_Tumour", "Endometrial_Tumour", "Breast_Normal", "Endometrial_Normal"]:
                sub = summary_df[summary_df["analysis_group"] == grp]
                if not sub.empty and pd.notna(sub.iloc[0].get("nonmissing")):
                    parts.append(f"{grp}: n={int(sub.iloc[0]['nonmissing'])}")
            if parts:
                ax_expr.text(0.02, 0.98, " | ".join(parts), transform=ax_expr.transAxes, ha="left", va="top", fontsize=7.2, color="#555555")
    else:
        ax_expr.axis("off")
        ax_expr.text(0.5, 0.5, "No non-missing expression values", ha="center", va="center", fontsize=10, color="#666666")

    _plot_top_hits(ax_snp, snp_df, "Top SNP associations", lambda r: f"{r.get('Cohort','')} | {r.get('Variant_ID','')}", "#4C78A8")
    _plot_top_hits(ax_hap, hap_df, "Top haplotype associations", lambda r: f"{r.get('Cohort','')} | {r.get('Haplotype_ID','')}", "#F58518")
    _plot_top_hits(ax_clin, clin_df, "Top clinical associations", lambda r: f"{r.get('Cohort','')} | {r.get('Clin_Label','')}", "#54A24B")

    fig.suptitle(f"Panel-gene report: {gene}", fontsize=14, fontweight="bold")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_panel_gene_significant_summary(summary_df: pd.DataFrame, out_path: Path) -> None:
    if summary_df.empty:
        return
    sig = summary_df[summary_df[["SNP_FDR_Sig_Count", "Haplotype_FDR_Sig_Count", "Clinical_FDR_Sig_Count"]].sum(axis=1) > 0].copy()
    if sig.empty:
        return
    sig["Best_P_Value"] = pd.to_numeric(sig["Best_P_Value"], errors="coerce")
    sig = sig.sort_values("Best_P_Value")
    top = sig.head(18).copy()
    top["Best_Score"] = -np.log10(top["Best_P_Value"].clip(lower=1e-300))
    domain_colours = {"SNP": "#4C78A8", "Haplotype": "#F58518", "Clinical": "#54A24B"}
    colours = [domain_colours.get(v, "#9E9E9E") for v in top["Best_Association_Type"]]

    fig = plt.figure(figsize=(14, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.28, 1.0], width_ratios=[1.05, 1.8], hspace=0.28, wspace=0.25)
    ax_counts = fig.add_subplot(gs[0, 0])
    ax_note = fig.add_subplot(gs[0, 1])
    ax_bar = fig.add_subplot(gs[1, :])

    count_map = {
        "SNP FDR": int((sig["SNP_FDR_Sig_Count"] > 0).sum()),
        "Haplotype FDR": int((sig["Haplotype_FDR_Sig_Count"] > 0).sum()),
        "Clinical FDR": int((sig["Clinical_FDR_Sig_Count"] > 0).sum()),
    }
    count_labels = list(count_map.keys())
    count_vals = list(count_map.values())
    count_cols = ["#4C78A8", "#F58518", "#54A24B"]
    ax_counts.bar(count_labels, count_vals, color=count_cols, alpha=0.88)
    ax_counts.set_title("Genes with at least one FDR-significant hit", fontsize=10, fontweight="bold")
    ax_counts.set_ylabel("Genes", fontsize=9)
    ax_counts.tick_params(axis="x", labelrotation=10, labelsize=8)
    for i, v in enumerate(count_vals):
        ax_counts.text(i, v + 0.3, str(v), ha="center", va="bottom", fontsize=8)
    for spine in ["top", "right"]:
        ax_counts.spines[spine].set_visible(False)

    ax_note.axis("off")
    ax_note.text(0, 0.85, "Exploratory panel-gene summary", fontsize=12, fontweight="bold")
    ax_note.text(0, 0.58, f"{len(sig)} of {len(summary_df)} genes showed at least one FDR-significant association across SNP, haplotype, or clinical analyses.", fontsize=10)
    ax_note.text(0, 0.28, "Bars below rank genes by best p-value. Colour shows the evidence type of the strongest association; badges indicate additional FDR-positive evidence types.", fontsize=9.2, color="#444444")

    ypos = np.arange(len(top))[::-1]
    ax_bar.barh(ypos, top["Best_Score"], color=colours, alpha=0.88)
    labels = []
    for _, row in top.iterrows():
        target = row.get("Best_Target_Label", "") or ""
        labels.append(f"{row['Gene']} | {target}")
    ax_bar.set_yticks(ypos)
    ax_bar.set_yticklabels(labels, fontsize=8.5)
    ax_bar.set_xlabel("-log10(best P value)", fontsize=10)
    ax_bar.set_title("Top FDR-significant panel genes", fontsize=11, fontweight="bold")
    ax_bar.grid(axis="x", linestyle=":", alpha=0.25)
    for spine in ["top", "right"]:
        ax_bar.spines[spine].set_visible(False)
    for y, (_, row) in zip(ypos, top.iterrows()):
        badges = []
        if row["SNP_FDR_Sig_Count"] > 0:
            badges.append("S")
        if row["Haplotype_FDR_Sig_Count"] > 0:
            badges.append("H")
        if row["Clinical_FDR_Sig_Count"] > 0:
            badges.append("C")
        ax_bar.text(row["Best_Score"] + 0.08, y, f"p={row['Best_P_Value']:.3g} [{' '.join(badges)}]", va="center", fontsize=7.8, color="#444444")

    fig.suptitle("Objective 2 exploratory panel-gene signals", fontsize=14, fontweight="bold")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def gene_series(df: pd.DataFrame, gene: str) -> pd.Series:
    value = df[gene]
    if isinstance(value, pd.DataFrame):
        return pd.to_numeric(value.iloc[:, 0], errors="coerce")
    return pd.to_numeric(value, errors="coerce")
def gene_summary(expr: pd.DataFrame, gene: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for group, sub in expr.groupby("analysis_group", dropna=False):
        values = gene_series(sub, gene)
        rows.append(
            {
                "Gene": gene,
                "analysis_group": group,
                "analysis_role": sub["analysis_role"].dropna().astype(str).iloc[0] if sub["analysis_role"].notna().any() else "",
                "n_rows": len(sub),
                "n_samples": sub["snp_code"].nunique(dropna=True),
                "nonmissing": int(values.notna().sum()),
                "mean": values.mean(skipna=True),
                "median": values.median(skipna=True),
                "std": values.std(skipna=True),
            }
        )
    return pd.DataFrame(rows).sort_values("analysis_group")


def assoc_genetic_gene(expr: pd.DataFrame, long_df: pd.DataFrame, id_col: str, analysis_type: str, gene: str, min_carriers: int, fdr_threshold: float) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    cohorts = {
        "Breast": expr[expr["cohort"] == "Breast"].copy(),
        "Endometrial": expr[expr["cohort"] == "Endometrial"].copy(),
        "All_Primary_Tumours": expr.copy(),
    }
    for cohort_name, cohort_df in cohorts.items():
        codes = set(cohort_df["snp_code"].dropna())
        subset = long_df[long_df["snp_code"].isin(codes)].copy()
        for key, grouped in subset.groupby(id_col, dropna=False):
            meta_cols = [col for col in [id_col, "snp_code", "Dosage", "Carrier", "SNP_Label", "gnomAD_NFE_AF", "Haplotype", "Global_Freq"] if col in grouped.columns]
            merged = cohort_df.merge(grouped[meta_cols], on="snp_code", how="inner")
            base_cols = [col for col in ["Dosage", "Carrier", "SNP_Label", "gnomAD_NFE_AF", "Haplotype", "Global_Freq"] if col in merged.columns]
            test_df = merged[base_cols].copy()
            test_df[gene] = gene_series(merged, gene)
            test_df = test_df.dropna(subset=[gene, "Dosage"])
            if len(test_df) < 10 or test_df["Dosage"].nunique() < 2 or int((test_df["Dosage"] > 0).sum()) < min_carriers:
                continue
            slope, intercept, r_value, p_value, se = stats.linregress(test_df["Dosage"], test_df[gene])
            row = {
                "Analysis_Type": analysis_type,
                "Cohort": cohort_name,
                id_col: key,
                "Gene": gene,
                "N": len(test_df),
                "Carrier_Count": int((test_df["Dosage"] > 0).sum()),
                "N_Dosage_0": int((test_df["Dosage"] == 0).sum()),
                "N_Dosage_1": int((test_df["Dosage"] == 1).sum()),
                "N_Dosage_2": int((test_df["Dosage"] == 2).sum()),
                "Slope": slope,
                "Intercept": intercept,
                "R": r_value,
                "SE": se,
                "P_Value": p_value,
            }
            if "Carrier" in test_df.columns:
                carrier_vals = pd.to_numeric(test_df.loc[test_df["Carrier"] == 1, gene], errors="coerce").dropna()
                noncarrier_vals = pd.to_numeric(test_df.loc[test_df["Carrier"] == 0, gene], errors="coerce").dropna()
                row["P_Carrier"] = stats.mannwhitneyu(carrier_vals, noncarrier_vals, alternative="two-sided")[1] if len(carrier_vals) >= min_carriers and len(noncarrier_vals) >= min_carriers else np.nan
                row["Carrier_Median_Diff"] = carrier_vals.median() - noncarrier_vals.median() if len(carrier_vals) and len(noncarrier_vals) else np.nan
            for col in ["SNP_Label", "gnomAD_NFE_AF", "Haplotype", "Global_Freq"]:
                if col in test_df.columns:
                    row[col] = test_df[col].dropna().astype(str).iloc[0] if test_df[col].notna().any() else ""
            rows.append(row)
    return add_fdr(pd.DataFrame(rows), "P_Value", ["Analysis_Type", "Cohort"], fdr_threshold) if rows else pd.DataFrame()

def clin_bin_gene(df: pd.DataFrame, y_col: str, x_col: str, bmi_adjust: int) -> dict[str, object]:
    out = {
        "OR_Adj_Age": np.nan,
        "OR_Adj_Age_CI95": np.nan,
        "P_Adj_Age": np.nan,
        "N_Adj_Age": np.nan,
        "Adj_Age_Note": "statsmodels not available",
        "OR_Adj_AgeBMI": np.nan,
        "OR_Adj_AgeBMI_CI95": np.nan,
        "P_Adj_AgeBMI": np.nan,
        "N_Adj_AgeBMI": np.nan,
        "Adj_AgeBMI_Note": "not run",
    }
    try:
        import statsmodels.formula.api as smf
    except ImportError:
        return out
    age_df = df[[y_col, x_col, "canon__age"]].dropna().copy()
    if len(age_df) >= 10 and age_df[y_col].nunique() == 2:
        age_df["endpoint_z"] = z(age_df[x_col])
        age_df["age_z"] = z(age_df["canon__age"])
        age_df = age_df.dropna()
        if len(age_df) >= 10:
            try:
                fit = smf.logit(f"{y_col} ~ endpoint_z + age_z", data=age_df).fit(disp=0, maxiter=200)
                ci = fit.conf_int().loc["endpoint_z"]
                out.update({
                    "OR_Adj_Age": float(np.exp(fit.params["endpoint_z"])),
                    "OR_Adj_Age_CI95": f"[{float(np.exp(ci.iloc[0])):.3f}, {float(np.exp(ci.iloc[1])):.3f}]",
                    "P_Adj_Age": float(fit.pvalues["endpoint_z"]),
                    "N_Adj_Age": len(age_df),
                    "Adj_Age_Note": "converged",
                })
            except Exception as exc:
                out["Adj_Age_Note"] = str(exc)[:80]
    if bmi_adjust:
        bmi_df = df[[y_col, x_col, "canon__age", "canon__bmi"]].dropna().copy()
        if len(bmi_df) >= 10 and bmi_df[y_col].nunique() == 2:
            bmi_df["endpoint_z"] = z(bmi_df[x_col])
            bmi_df["age_z"] = z(bmi_df["canon__age"])
            bmi_df["bmi_z"] = z(bmi_df["canon__bmi"])
            bmi_df = bmi_df.dropna()
            if len(bmi_df) >= 10:
                try:
                    fit = smf.logit(f"{y_col} ~ endpoint_z + age_z + bmi_z", data=bmi_df).fit(disp=0, maxiter=200)
                    ci = fit.conf_int().loc["endpoint_z"]
                    out.update({
                        "OR_Adj_AgeBMI": float(np.exp(fit.params["endpoint_z"])),
                        "OR_Adj_AgeBMI_CI95": f"[{float(np.exp(ci.iloc[0])):.3f}, {float(np.exp(ci.iloc[1])):.3f}]",
                        "P_Adj_AgeBMI": float(fit.pvalues["endpoint_z"]),
                        "N_Adj_AgeBMI": len(bmi_df),
                        "Adj_AgeBMI_Note": "converged",
                    })
                except Exception as exc:
                    out["Adj_AgeBMI_Note"] = str(exc)[:80]
    return out


def assoc_clin_gene(expr: pd.DataFrame, gene: str, breast_cfg: dict, endo_cfg: dict, min_carriers: int, fdr_threshold: float) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cohort_name, cfg in [("Breast", breast_cfg), ("Endometrial", endo_cfg)]:
        cohort_df = expr[expr["cohort"] == cohort_name].copy()
        for clin_var, (clin_type, label, bmi_adjust) in cfg.items():
            if clin_var not in cohort_df.columns:
                continue
            test_df = cohort_df[[clin_var, "canon__age", "canon__bmi"]].copy()
            test_df[gene] = gene_series(cohort_df, gene)
            test_df = test_df.dropna(subset=[gene, clin_var])
            if len(test_df) < 10:
                continue
            row = {
                "Analysis_Type": "PanelGene_vs_Clinical",
                "Cohort": cohort_name,
                "Gene": gene,
                "Clin_Var": clin_var,
                "Clin_Label": label,
                "Clin_Type": clin_type,
                "N": len(test_df),
                "P_Value": np.nan,
                "Effect": np.nan,
                "Test": "",
                "OR_Adj_Age": np.nan,
                "OR_Adj_Age_CI95": np.nan,
                "P_Adj_Age": np.nan,
                "N_Adj_Age": np.nan,
                "Adj_Age_Note": "",
                "OR_Adj_AgeBMI": np.nan,
                "OR_Adj_AgeBMI_CI95": np.nan,
                "P_Adj_AgeBMI": np.nan,
                "N_Adj_AgeBMI": np.nan,
                "Adj_AgeBMI_Note": "",
            }
            if clin_type == "continuous":
                test_df[clin_var] = pd.to_numeric(test_df[clin_var], errors="coerce")
                test_df = test_df.dropna(subset=[clin_var])
                if len(test_df) < 10:
                    continue
                rho, p_value = stats.spearmanr(test_df[gene], test_df[clin_var], nan_policy="omit")
                row.update({"Test": "Spearman", "P_Value": p_value, "Effect": rho})
            elif clin_type == "binary":
                test_df[clin_var] = pd.to_numeric(test_df[clin_var], errors="coerce")
                test_df = test_df[test_df[clin_var].isin([0, 1])].copy()
                if test_df.empty or test_df[clin_var].nunique() < 2:
                    continue
                pos_vals = pd.to_numeric(test_df.loc[test_df[clin_var] == 1, gene], errors="coerce").dropna()
                neg_vals = pd.to_numeric(test_df.loc[test_df[clin_var] == 0, gene], errors="coerce").dropna()
                if len(pos_vals) < min_carriers or len(neg_vals) < min_carriers:
                    continue
                row.update({
                    "Test": "Mann-Whitney U",
                    "P_Value": stats.mannwhitneyu(pos_vals, neg_vals, alternative="two-sided")[1],
                    "Effect": pos_vals.median() - neg_vals.median(),
                })
                row.update(clin_bin_gene(test_df, clin_var, gene, bmi_adjust))
            else:
                test_df[clin_var] = test_df[clin_var].astype(str).str.strip()
                groups = [pd.to_numeric(group[gene], errors="coerce").dropna() for _, group in test_df.groupby(clin_var)]
                groups = [group for group in groups if len(group) >= 3]
                if len(groups) < 2:
                    continue
                statistic, p_value = stats.kruskal(*groups)
                row.update({"Test": "Kruskal-Wallis", "P_Value": p_value, "Effect": statistic})
            rows.append(row)
    return add_fdr(pd.DataFrame(rows), "P_Value", ["Analysis_Type", "Cohort"], fdr_threshold) if rows else pd.DataFrame()


def export_panel_gene_reports(merged: pd.DataFrame, primary: pd.DataFrame, genes: list[str], genetic_long: pd.DataFrame, hap_carrier: pd.DataFrame, out_dir: str | Path, breast_cfg: dict, endo_cfg: dict, min_carriers: int, fdr_threshold: float, manifest_cols: list[str]) -> pd.DataFrame:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    figure_dir = out_path / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    expected_names = {f"{safe_gene_slug(gene)}.xlsx" for gene in genes}
    expected_figs = {f"{safe_gene_slug(gene)}.png" for gene in genes}
    for old_workbook in out_path.glob("*.xlsx"):
        if old_workbook.name not in expected_names:
            old_workbook.unlink()
    for old_figure in figure_dir.glob("*.png"):
        if old_figure.name not in expected_figs and old_figure.name != "20_PanelGene_Significant_Summary.png":
            old_figure.unlink()
    manifest_rows: list[dict[str, object]] = []
    for gene in genes:
        slug = safe_gene_slug(gene)
        workbook_path = out_path / f"{slug}.xlsx"
        figure_path = figure_dir / f"{slug}.png"
        gene_values = merged[[col for col in manifest_cols if col in merged.columns]].copy()
        gene_values[gene] = gene_series(merged, gene)
        gene_values = gene_values.sort_values([col for col in ["analysis_group", "snp_code", "CODIGO JC"] if col in gene_values.columns])
        summary_df = gene_summary(merged, gene)
        snp_df = assoc_genetic_gene(primary, genetic_long, "Variant_ID", "SNP_vs_PanelGene", gene, min_carriers, fdr_threshold)
        hap_df = assoc_genetic_gene(primary, hap_carrier, "Haplotype_ID", "Haplotype_vs_PanelGene", gene, min_carriers, fdr_threshold)
        clin_df = assoc_clin_gene(primary, gene, breast_cfg, endo_cfg, min_carriers, fdr_threshold)
        with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
            gene_values.to_excel(writer, sheet_name="expression_values", index=False)
            summary_df.to_excel(writer, sheet_name="cohort_summary", index=False)
            (snp_df if not snp_df.empty else pd.DataFrame({"Note": ["No SNP associations met the minimum thresholds for this gene."]})).to_excel(writer, sheet_name="snp_assoc", index=False)
            (hap_df if not hap_df.empty else pd.DataFrame({"Note": ["No haplotype associations met the minimum thresholds for this gene."]})).to_excel(writer, sheet_name="haplotype_assoc", index=False)
            (clin_df if not clin_df.empty else pd.DataFrame({"Note": ["No clinical associations met the minimum thresholds for this gene."]})).to_excel(writer, sheet_name="clinical_assoc", index=False)
        make_panel_gene_figure(gene, merged, summary_df, snp_df, hap_df, clin_df, figure_path)

        best_rows = []
        if not snp_df.empty and snp_df["P_Value"].notna().any():
            r = snp_df.loc[snp_df["P_Value"].idxmin()]
            best_rows.append(("SNP", float(r["P_Value"]), f"{r.get('Cohort','')} | {r.get('Variant_ID','')}"))
        if not hap_df.empty and hap_df["P_Value"].notna().any():
            r = hap_df.loc[hap_df["P_Value"].idxmin()]
            best_rows.append(("Haplotype", float(r["P_Value"]), f"{r.get('Cohort','')} | {r.get('Haplotype_ID','')}"))
        if not clin_df.empty and clin_df["P_Value"].notna().any():
            r = clin_df.loc[clin_df["P_Value"].idxmin()]
            best_rows.append(("Clinical", float(r["P_Value"]), f"{r.get('Cohort','')} | {r.get('Clin_Label','')}"))
        best_type, best_p, best_label = (None, np.nan, "") if not best_rows else sorted(best_rows, key=lambda x: x[1])[0]

        manifest_rows.append(
            {
                "Gene": gene,
                "Workbook": str(workbook_path),
                "Figure": str(figure_path),
                "Rows": len(gene_values),
                "Primary_Tumour_N": int(gene_series(primary, gene).dropna().shape[0]) if gene in primary.columns else 0,
                "SNP_Association_Rows": len(snp_df),
                "Haplotype_Association_Rows": len(hap_df),
                "Clinical_Association_Rows": len(clin_df),
                "SNP_FDR_Sig_Count": int(snp_df.get("FDR_Sig", pd.Series(dtype=bool)).fillna(False).sum()) if not snp_df.empty else 0,
                "Haplotype_FDR_Sig_Count": int(hap_df.get("FDR_Sig", pd.Series(dtype=bool)).fillna(False).sum()) if not hap_df.empty else 0,
                "Clinical_FDR_Sig_Count": int(clin_df.get("FDR_Sig", pd.Series(dtype=bool)).fillna(False).sum()) if not clin_df.empty else 0,
                "Best_Association_Type": best_type,
                "Best_P_Value": best_p,
                "Best_Target_Label": best_label,
            }
        )
    manifest_df = pd.DataFrame(manifest_rows).sort_values("Gene")
    make_panel_gene_significant_summary(manifest_df, figure_dir / "20_PanelGene_Significant_Summary.png")
    return manifest_df


