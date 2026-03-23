from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


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
    expected_names = {f"{safe_gene_slug(gene)}.xlsx" for gene in genes}
    for old_workbook in out_path.glob("*.xlsx"):
        if old_workbook.name not in expected_names:
            old_workbook.unlink()
    manifest_rows: list[dict[str, object]] = []
    for gene in genes:
        slug = safe_gene_slug(gene)
        workbook_path = out_path / f"{slug}.xlsx"
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
        manifest_rows.append(
            {
                "Gene": gene,
                "Workbook": str(workbook_path),
                "Rows": len(gene_values),
                "Primary_Tumour_N": int(gene_series(primary, gene).dropna().shape[0]) if gene in primary.columns else 0,
                "SNP_Association_Rows": len(snp_df),
                "Haplotype_Association_Rows": len(hap_df),
                "Clinical_Association_Rows": len(clin_df),
            }
        )
    return pd.DataFrame(manifest_rows).sort_values("Gene")



