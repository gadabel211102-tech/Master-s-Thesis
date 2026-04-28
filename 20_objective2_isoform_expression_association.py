#!/usr/bin/env python3
"""Objective 2 workbook-first RNA isoform/expression analysis.

Pipeline-facing text is written in British English. Some American-spelt tissue
labels are still accepted below because they appear in upstream source files and
must remain supported for robust sample matching.
"""

from __future__ import annotations
import argparse,re,shutil,warnings
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
try:
    import seaborn as sns
except ImportError:
    sns = None
from scipy import stats
from association_runtime import script20_defaults
from figure_style import COHORT_COLORS, TUMOUR_CONTROL_COLORS
from pipeline_utils import build_genomic_variant_id_series, combine_gnomad_nfe, common_nfe_variant_mask
from objective2_rna_qc_utils import DNA_BASELINE_COLUMNS, RNA_QC_COLUMNS, attach_rna_qc, bam_validation, plot_rna_qc
from objective2_gene_report_exports import export_panel_gene_reports
from sample_identity_utils import build_analysis_sample_map
warnings.filterwarnings("ignore")
D=script20_defaults(); MIN_NFE=float(D["min_nfe_af"]); MIN_HAP=float(D["min_hap_freq"]); MIN_C=int(D["min_carriers"]); FDR=float(D["fdr_threshold"])
ISO=["G1","G2","G3","G3b","G4"]; ENDPTS=ISO+["GSDMB","pyroptotic_isoform_fraction","non_pyroptotic_isoform_fraction","G4_vs_total","G2_vs_total"]
TMAP={"breast tumor":("Breast","Tumour","Breast_Tumour","primary_tumour"),"breast tumour":("Breast","Tumour","Breast_Tumour","primary_tumour"),"breast normal":("Breast","Healthy","Breast_Normal","context_control"),"endometrial cancer":("Endometrial","Tumour","Endometrial_Tumour","primary_tumour"),"endometrial normal":("Endometrial","Healthy","Endometrial_Normal","context_control"),"endometrial_normal":("Endometrial","Healthy","Endometrial_Normal","context_control"),"ovarian":("Ovarian","Tumour","Ovarian","excluded"),"ovarian organoids":("Ovarian","Tumour","Ovarian_Organoids","excluded"),"cell line":("CellLine","Unknown","Cell_Line","excluded")}
CB_CORE={"BREAST_RECURRENCE_DERIVED":("binary","Recurrence / progression",1),"BREAST_METASTASIS_DERIVED":("binary","Distant metastasis",1),"BREAST_ANY_METASTASIS_BIN":("binary","Any metastasis",1),"BREAST_EXITUS_DERIVED":("binary","Death",1),"BREAST_OS_MONTHS_DERIVED":("continuous","Overall survival (months)",1)}
CB_SUPPLEMENTARY={"BREAST_HER2_SUBTYPE":("nominal","HER2 subtype",0),"BREAST_DX_TYPE":("nominal","Histological diagnosis",0),"BREAST_P53_NUMERIC":("continuous","p53",0)}
CE_CORE={"ENDO_FIGO_NUMERIC":("continuous","FIGO stage",1),"ENDO_GRADE_NUMERIC":("continuous","Tumour grade",1),"ENDO_RISK_ORDINAL":("continuous","Risk of recurrence",1),"ENDO_LVSI_BIN":("binary","LVSI",0),"ENDO_MYOINV_BIN":("binary","Myometrial invasion ?50%",1),"ENDO_PD_BIN":("binary","Disease progression",1),"ENDO_EXITUS_BIN":("binary","All-cause death",1),"ENDO_OS_MONTHS_DERIVED":("continuous","Overall survival (months)",1),"ENDO_PFS_MONTHS_DERIVED":("continuous","Progression-free survival (months)",1),"canon__distant_mets_flag":("binary","Distant metastasis",1)}
CE_SUPPLEMENTARY={"ENDO_MSI_BIN":("binary","MSI-H",0),"ENDO_TP53_ABN_BIN":("binary","TP53 IHC abnormal",0),"ENDO_NEEC_BIN":("binary","Non-endometrioid histology",0),"ENDO_N_STAGE_BIN":("binary","Lymph node involvement",1),"ENDO_ER_BIN":("binary","ER positive",0),"ENDO_PR_BIN":("binary","PR positive",0),"ENDO_PDL1_NUMERIC":("continuous","PD-L1",0),"ENDO_CD8_NUMERIC":("continuous","CD8+ TILs",0),"ENDO_PTEN_BIN":("binary","PTEN loss/reduced",0),"ENDO_MLH1_BIN":("binary","MLH1 loss/reduced",0),"ENDO_GENE_AMP_BIN":("binary","Gene amplification",0),"canon__molecular_class":("nominal","Molecular classification",0)}
BREAST_THERAPY={"canon__treatment_chemotherapy":("binary","Chemotherapy received",0),"canon__treatment_anti_her2":("binary","Anti-HER2 therapy received",0),"canon__treatment_endocrine":("binary","Endocrine therapy received",0),"canon__treatment_radiotherapy":("binary","Radiotherapy received",0)}
ENDO_THERAPY={"clin_au_endo__PD_TREATMENT":("nominal","Progression treatment",0)}
IMMUNE_SIGNATURE_GENES=("PTPRC","CD14","CD68","PDCD1","CD274","CD8A","CD8B","GZMB","CXCL9","CXCL10","LAG3","CTLA4","TIGIT","IFNG")
CB=CB_CORE
CE=CE_CORE

def clinical_var_sets(include_supplementary: bool = False):
    if include_supplementary:
        breast = {**CB_CORE, **CB_SUPPLEMENTARY}
        endo = {**CE_CORE, **CE_SUPPLEMENTARY}
        return breast, endo
    return CB_CORE, CE_CORE

def therapy_var_sets():
    return BREAST_THERAPY, ENDO_THERAPY

s=lambda v: "" if pd.isna(v) else str(v).strip()
nk=lambda v: re.sub(r"[^A-Z0-9]+","",s(v).upper())
def first(x):
    for v in x:
        if pd.notna(v) and str(v).strip() not in {"","nan","None"}: return v
    return np.nan
def z(x):
    x=pd.to_numeric(x,errors="coerce"); sd=x.std(ddof=0)
    return (x-x.mean())/sd if pd.notna(sd) and sd!=0 else pd.Series(np.nan,index=x.index)
def bh(vals):
    a=np.asarray(vals,dtype=float); out=np.full(a.shape,np.nan); m=np.isfinite(a)
    if not m.any(): return out
    p=a[m]; o=np.argsort(p); r=p[o]; n=len(r); adj=np.empty(n); prev=1.0
    for i in range(n-1,-1,-1): prev=min(prev,r[i]*n/(i+1)); adj[i]=min(prev,1.0)
    b=np.empty(n); b[o]=adj; out[m]=b; return out
def add_fdr(df,p,grp):
    if df.empty: return df
    df=df.copy(); df["FDR"]=np.nan
    for _,ix in df.groupby(list(grp),dropna=False).groups.items(): df.loc[ix,"FDR"]=bh(pd.to_numeric(df.loc[ix,p],errors="coerce"))
    df["Nominal_Sig"]=pd.to_numeric(df[p],errors="coerce")<0.05; df["FDR_Sig"]=pd.to_numeric(df["FDR"],errors="coerce")<FDR
    return df
def sx(v):
    t=s(v)
    for pat,fmt in [(r"^(SNP_(?:AT|EN|MN|MT-T)_\d+)",lambda m:m.group(1).upper()),(r"^SNP_DNA_(EN|MN|AT)_(\d+)",lambda m:f"SNP_{m.group(1).upper()}_{m.group(2)}"),(r"^SNP_DNA_MT[-_]T_(\d+)",lambda m:f"SNP_MT-T_{m.group(1)}"),(r"^DNA_SNP_(EN|MN|AT)_(\d+)",lambda m:f"SNP_{m.group(1).upper()}_{m.group(2)}"),(r"^DNA_(AT|EN|MN)_(\d+)",lambda m:f"SNP_{m.group(1).upper()}_{m.group(2)}"),(r"^DNA_(?:SNP_)?MT[-_]T_(\d+)",lambda m:f"SNP_MT-T_{m.group(1)}"),(r"^MAMAH2_MT[-_]T_(\d+)",lambda m:f"SNP_MT-T_{m.group(1)}"),(r"^MAMAH2_MT[-_]N_(\d+)",lambda m:f"SNP_MN_{m.group(1)}"),(r"^RNA_SNP_(AT|EN|MN|MT-T)_(\d+)",lambda m:f"SNP_{m.group(1).upper()}_{m.group(2)}"),(r"^SNP_RNA_(AT|EN|MN|MT-T)_(\d+)",lambda m:f"SNP_{m.group(1).upper()}_{m.group(2)}"),(r"^RNA_(AT|EN|MN)_(\d+)",lambda m:f"SNP_{m.group(1).upper()}_{m.group(2)}"),(r"^RNA_MT[-_]T_(\d+)",lambda m:f"SNP_MT-T_{m.group(1)}")]:
        m=re.match(pat,t,re.I)
        if m: return fmt(m)
    return None
def aliases(v):
    t=s(v); out={nk(t)} if t else set(); c=sx(t)
    if c: out.add(nk(c))
    m=re.match(r"^ENDOMDA[_ -]?(\d+)$",nk(t))
    if m: out.add(nk(f"MDA{int(m.group(1)):03d}"))
    for pat in [r"(?:ENDO)?(BL|MDA)[_ -]?0*(\d+)$",r"(?:^|_)(BL|MDA)[_ -]?0*(\d+)"]:
        m=re.search(pat,t,re.I) or re.search(pat,nk(t),re.I)
        if m:
            out.add(nk(f"{m.group(1).upper()}{int(m.group(2))}")); out.add(nk(f"{m.group(1).upper()}{int(m.group(2)):03d}"))
    return {x for x in out if x}
def gt_dose(gt):
    t=s(gt)
    if t in {"",".",".|.","./."}: return np.nan
    parts=t.split("|") if "|" in t else t.split("/") if "/" in t else None
    try: return float(sum(int(x) for x in parts)) if parts else np.nan
    except: return np.nan
def gt_haps(gt):
    t=s(gt)
    if t in {"",".",".|.","./."}: return "N","N"
    if "|" in t: a,b=t.split("|",1); return a,b
    if "/" in t: a,b=t.split("/",1); return a,b
    return "N","N"

def pretty_group(label):
    label = s(label).replace('_', ' ')
    parts = label.split()
    if len(parts) >= 2:
        return '\n'.join([' '.join(parts[:-1]), parts[-1]])
    return label

def short_label(label, max_chars=26):
    label = s(label)
    if len(label) <= max_chars:
        return label
    keep = max(6, (max_chars - 1) // 2)
    tail = max(6, max_chars - keep - 1)
    return f"{label[:keep]}...{label[-tail:]}"

def primary_rsid(*values):
    for value in values:
        if pd.isna(value):
            continue
        hits = re.findall(r"rs\d+", str(value), flags=re.I)
        if hits:
            return hits[0].lower()
    return np.nan

def variant_display_label(variant_id, rsid, snp_label=None):
    rid = s(rsid)
    if rid:
        return rid
    lab = s(snp_label)
    if lab:
        return lab
    return s(variant_id)

def summarise_group_counts(df, group_col):
    counts = df[group_col].value_counts(dropna=False).to_dict()
    return counts

def heat_annotation_symbol(p_value, fdr_value):
    p = pd.to_numeric(pd.Series([p_value]), errors="coerce").iloc[0]
    fdr = pd.to_numeric(pd.Series([fdr_value]), errors="coerce").iloc[0]
    if pd.notna(fdr) and fdr < FDR:
        return "**"
    if pd.notna(p) and p < 0.05:
        return "*"
    return ""

def format_p_value(p):
    if pd.isna(p):
        return "NA"
    p = float(p)
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.3f}"

def pairwise_mwu_pvals(df, endpoint, group_col, group_pairs):
    out = []
    for g1, g2, label in group_pairs:
        a = pd.to_numeric(df.loc[df[group_col] == g1, endpoint], errors="coerce").dropna()
        b = pd.to_numeric(df.loc[df[group_col] == g2, endpoint], errors="coerce").dropna()
        if len(a) < 2 or len(b) < 2:
            continue
        try:
            p = stats.mannwhitneyu(a, b, alternative="two-sided")[1]
        except Exception:
            p = np.nan
        out.append(f"{label}: p={format_p_value(p)}")
    return out

def load_expr(path):
    xls=pd.ExcelFile(path)
    if "Sheet3" not in xls.sheet_names or "eva" not in xls.sheet_names: raise ValueError("Workbook needs Sheet3 and eva")
    df=pd.read_excel(path,sheet_name="Sheet3"); eva=pd.read_excel(path,sheet_name="eva")
    meta=df["Tissue"].astype(str).str.strip().str.lower().map(TMAP)
    df["cohort"]=meta.map(lambda x:x[0] if isinstance(x,tuple) else "Unknown"); df["tumour_normal"]=meta.map(lambda x:x[1] if isinstance(x,tuple) else "Unknown")
    df["analysis_group"]=meta.map(lambda x:x[2] if isinstance(x,tuple) else "Unknown"); df["analysis_role"]=meta.map(lambda x:x[3] if isinstance(x,tuple) else "excluded")
    key=[c for c in ["NOMBRE DE LA MUESTRA","CODIGO JC"] if c in eva.columns]; sup=[c for c in ["observaciones EVA","TIENEN RNA","FALTA MUESTRA"] if c in eva.columns]
    if key: df=df.merge(eva[key+sup].drop_duplicates(),on=key,how="left")
    for c in set(ISO+["GSDMB"]): df[c]=pd.to_numeric(df.get(c),errors="coerce")
    df["pyroptotic_isoform_fraction"]=df[["G3","G3b","G4"]].sum(axis=1,min_count=1); df["non_pyroptotic_isoform_fraction"]=df[["G1","G2"]].sum(axis=1,min_count=1); df["isoform_total_sum"]=df[ISO].sum(axis=1,min_count=1)
    med=pd.to_numeric(df["isoform_total_sum"],errors="coerce").median(skipna=True); lo,hi=(95,105) if pd.notna(med) and med>2 else (0.95,1.05)
    df["isoform_total_scale"]="percentage" if pd.notna(med) and med>2 else "fraction"; df["isoform_total_flag"]=df["isoform_total_sum"].notna()&~df["isoform_total_sum"].between(lo,hi)
    df["G4_vs_total"]=np.where(df["isoform_total_sum"]>0,df["G4"]/df["isoform_total_sum"],np.nan); df["G2_vs_total"]=np.where(df["isoform_total_sum"]>0,df["G2"]/df["isoform_total_sum"],np.nan)
    df["endpoint_missing_count"]=df[["G1","G2","G3","G3b","G4","GSDMB","pyroptotic_isoform_fraction"]].isna().sum(axis=1)
    df["raw1"]=df["CODIGO JC"].fillna(df["NOMBRE DE LA MUESTRA"]); df["raw2"]=df["NOMBRE DE LA MUESTRA"].fillna(df["CODIGO JC"]); df["sx1"]=df["raw1"].map(sx); df["sx2"]=df["raw2"].map(sx)
    return df
def load_master(path):
    m=pd.read_excel(path,sheet_name="harmonised_plus_canon"); m["snp_code"]=m["snp_code"].astype(str).str.strip(); m["nucleic_acid"]=m["nucleic_acid"].astype(str).str.upper(); return m
def match_expr(df,m):
    mr=(m[m["nucleic_acid"]=="RNA"].sort_values(["snp_code","sample_id","case_id"]).groupby("snp_code",dropna=False).agg(first).reset_index(drop=False))
    amap={}
    for _,r in mr.iterrows():
        for a in aliases(r.get("snp_code"))|aliases(r.get("sample_id"))|aliases(r.get("case_id")): amap.setdefault(a,set()).add(r["snp_code"])
    rows=[]; aud=[]; codes=set(mr["snp_code"])
    for _,r in df.iterrows():
        cand={x for x in [r.get("sx1"),r.get("sx2")] if x in codes}; meth="snp_code_regex" if len(cand)==1 else ""
        if not cand:
            ac=aliases(r.get("raw1"))|aliases(r.get("raw2")); cand={c for a in ac for c in amap.get(a,set())}; meth="sample_alias" if len(cand)==1 else "ambiguous_alias" if len(cand)>1 else "unmatched"
        elif len(cand)>1: meth="ambiguous_regex"
        out=r.to_dict(); out["snp_code"]=next(iter(cand)) if len(cand)==1 else np.nan; out["match_status"]="matched" if len(cand)==1 else "ambiguous" if len(cand)>1 else "unmatched"; out["match_method"]=meth; out["match_candidates"]="; ".join(sorted(cand)) if len(cand)>1 else ""
        rows.append(out)
        if out["match_status"]!="matched": aud.append({"NOMBRE DE LA MUESTRA":r.get("NOMBRE DE LA MUESTRA"),"CODIGO JC":r.get("CODIGO JC"),"Tissue":r.get("Tissue"),"analysis_group":r.get("analysis_group"),"match_status":out["match_status"],"match_method":meth,"match_candidates":out["match_candidates"]})
    merged=pd.DataFrame(rows).merge(mr,on="snp_code",how="left",suffixes=("","_master")); merged["master_join_success"]=merged["sample_id"].notna(); merged["analysis_include_primary"]=(merged["analysis_role"]=="primary_tumour")&merged["master_join_success"]
    return merged,pd.DataFrame(aud)
def load_backbone(path):
    # Prefer the explicit variant-level sheets in the annotated workbook so the
    # backbone does not depend on whichever sheet happens to be first.
    workbook = pd.ExcelFile(path)
    preferred_sheets = ["Biological_Annotations", "Variant_Annotations", "Sample_Variant_Annotations"]
    chosen_sheet = next((sheet for sheet in preferred_sheets if sheet in workbook.sheet_names), None)
    if chosen_sheet is None:
        non_readme = [sheet for sheet in workbook.sheet_names if str(sheet).strip().lower() != "readme"]
        chosen_sheet = non_readme[0] if non_readme else workbook.sheet_names[0]
    a=pd.read_excel(path, sheet_name=chosen_sheet)
    if "gnomADe_NFE_AF" in a.columns or "gnomADg_NFE_AF" in a.columns: a=combine_gnomad_nfe(a); fc="gnomAD_NFE_AF_combined"
    else:
        c=[x for x in a.columns if re.search(r"gnomad.*nfe.*af|nfe.*af",str(x),re.I)]
        if not c: raise ValueError("No gnomAD NFE AF column found in variant workbook")
        fc=c[0]
    a[fc]=pd.to_numeric(a[fc],errors="coerce"); a["POS"]=pd.to_numeric(a["POS"],errors="coerce"); a["REF"]=a["REF"].astype(str).str.strip().str.upper(); a["ALT"]=a["ALT"].astype(str).str.strip().str.upper(); a["SNP_Label"]=a["POS"].astype("Int64").astype(str)+"_"+a["REF"]+">"+a["ALT"]
    a["Variant_ID"]=build_genomic_variant_id_series(a["Variant_Key"] if "Variant_Key" in a.columns else None, a["CHROM"] if "CHROM" in a.columns else None, a["POS"], a["REF"], a["ALT"])
    if "rsID" not in a.columns:
        a["rsID"] = np.nan
    if "Existing_variation" not in a.columns:
        a["Existing_variation"] = np.nan
    a["Primary_rsID"] = [primary_rsid(rsid, existing) for rsid, existing in zip(a["rsID"], a["Existing_variation"])]
    a["Variant_Display"] = [variant_display_label(vid, rsid, lab) for vid, rsid, lab in zip(a["Variant_ID"], a["Primary_rsID"], a["SNP_Label"])]
    b=(a[common_nfe_variant_mask(a, MIN_NFE)].dropna(subset=["POS"]).sort_values(["POS","ALT"]).drop_duplicates("SNP_Label").copy()); b["POS"]=b["POS"].astype(int)
    keep=["POS","REF","ALT","SNP_Label","Variant_ID","Primary_rsID","Variant_Display",fc]+(["rsID"] if "rsID" in b.columns else [])+(["Existing_variation"] if "Existing_variation" in b.columns else [])
    return b[keep].rename(columns={fc:"gnomAD_NFE_AF"})
def load_phased(path,back):
    g=pd.read_csv(path,sep="	",dtype=str); g["POS"]=pd.to_numeric(g["POS"],errors="coerce"); g["REF"]=g["REF"].astype(str).str.strip().str.upper(); g["ALT"]=g["ALT"].astype(str).str.strip().str.upper(); g["SNP_Label"]=g["POS"].astype("Int64").astype(str)+"_"+g["REF"]+">"+g["ALT"]
    g=g[g["SNP_Label"].isin(set(back["SNP_Label"]))].merge(back[["SNP_Label","Variant_ID","Primary_rsID","Variant_Display","gnomAD_NFE_AF"]],on="SNP_Label",how="left").sort_values("POS").reset_index(drop=True)
    if g.empty: raise ValueError("No overlap between backbone SNPs and phased table")
    sc=[c for c in g.columns if c not in {"CHROM","POS","ID","REF","ALT","SNP_Label","Variant_ID","gnomAD_NFE_AF"}]
    identity_map = build_analysis_sample_map(sc, path, sx)
    representative_cols = identity_map.sort_values(["analysis_sample_id","raw_sample_name"]).drop_duplicates("analysis_sample_id")[["raw_sample_name","analysis_sample_id","snp_code"]]
    rows=[]; h=[]
    for _,id_row in representative_cols.iterrows():
        c=id_row["raw_sample_name"]; code=id_row["snp_code"]; analysis_sample_id=id_row["analysis_sample_id"]
        if c not in g.columns or not code: continue
        h1=[]; h2=[]
        for _,r in g.iterrows():
            a,b=gt_haps(r[c]); h1.append(a); h2.append(b); rows.append({"analysis_sample_id":analysis_sample_id,"Sample_phased":c,"snp_code":code,"Variant_ID":r["Variant_ID"],"Primary_rsID":r.get("Primary_rsID"),"Variant_Display":r.get("Variant_Display"),"SNP_Label":r["SNP_Label"],"POS":r["POS"],"Dosage":gt_dose(r[c]),"GT":r[c],"gnomAD_NFE_AF":r["gnomAD_NFE_AF"]})
        h.append({"analysis_sample_id":analysis_sample_id,"Sample_phased":c,"snp_code":code,"hap1":"".join(h1),"hap2":"".join(h2),"callable_haplotype":int("N" not in h1 and "N" not in h2)})
    gl=pd.DataFrame(rows).drop_duplicates(["analysis_sample_id","Variant_ID"])
    gw=gl.pivot_table(index="analysis_sample_id",columns="Variant_ID",values="Dosage",aggfunc="first").reset_index()
    if not gw.empty:
        gw=gw.merge(representative_cols[["analysis_sample_id","snp_code"]].drop_duplicates("analysis_sample_id"),on="analysis_sample_id",how="left")
    hd=pd.DataFrame(h)
    allh={}
    for col in ["hap1","hap2"]:
        for hp in hd[col].dropna():
            if hp and "N" not in hp: allh[hp]=allh.get(hp,0)+1
    tot=sum(allh.values()); hf=pd.DataFrame([{"Haplotype":k,"Count":v,"Global_Freq":v/tot if tot else np.nan} for k,v in sorted(allh.items(),key=lambda kv:(-kv[1],kv[0]))])
    if hf.empty: return gl,gw,hd,hf,pd.DataFrame()
    hf=hf[hf["Global_Freq"]>=MIN_HAP].reset_index(drop=True); hf["Haplotype_ID"]=[f"H{i+1}" for i in range(len(hf))]
    hc=[]
    for _,fr in hf.iterrows():
        for _,r in hd.iterrows():
            ok=int("N" not in str(r["hap1"]) and "N" not in str(r["hap2"])); dose=np.nan if not ok else int(str(r["hap1"])==fr["Haplotype"])+int(str(r["hap2"])==fr["Haplotype"])
            hc.append({"analysis_sample_id":r["analysis_sample_id"],"snp_code":r["snp_code"],"Sample_phased":r["Sample_phased"],"Haplotype_ID":fr["Haplotype_ID"],"Haplotype":fr["Haplotype"],"Global_Freq":fr["Global_Freq"],"Dosage":dose,"Carrier":np.nan if not ok else int(dose>=1),"Callable":ok})
    return gl,gw,hd,hf,pd.DataFrame(hc)
def summarise(df):
    rows=[]
    for g,sub in df.groupby("analysis_group",dropna=False):
        row={"analysis_group":g,"analysis_role":first(sub["analysis_role"]),"n_samples":sub.get("analysis_sample_id", sub["snp_code"]).nunique(dropna=True),"n_rows":len(sub),"n_matched_master":int(sub["master_join_success"].sum()),"n_flagged_isoform_total":int(sub["isoform_total_flag"].sum())}
        for e in ENDPTS:
            v=pd.to_numeric(sub[e],errors="coerce"); row[f"{e}__nonmissing"]=int(v.notna().sum()); row[f"{e}__median"]=v.median(skipna=True); row[f"{e}__mean"]=v.mean(skipna=True)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("analysis_group")

def assoc_genetic(expr,longdf,idcol,pcol):
    rows=[]; cohorts={"Breast":expr[expr["cohort"]=="Breast"].copy(),"Endometrial":expr[expr["cohort"]=="Endometrial"].copy(),"All_Primary_Tumours":expr.copy()}
    for cname,cdf in cohorts.items():
        codes=set(cdf["snp_code"].dropna()); sdf=longdf[longdf["snp_code"].isin(codes)].copy()
        for key,g in sdf.groupby(idcol,dropna=False):
            meta=[c for c in [idcol,"analysis_sample_id","snp_code","Dosage","Carrier","SNP_Label","Primary_rsID","Variant_Display","gnomAD_NFE_AF","Haplotype","Global_Freq"] if c in g.columns]
            dedup_keys=[c for c in ["analysis_sample_id"] if c in meta]
            m=cdf.merge(g[meta].drop_duplicates(dedup_keys if dedup_keys else None),on="snp_code",how="inner")
            if m.empty: continue
            for e in ENDPTS:
                t=m[[c for c in [e,"Dosage","Carrier","SNP_Label","Primary_rsID","Variant_Display","gnomAD_NFE_AF","Haplotype","Global_Freq"] if c in m.columns]].dropna(subset=[e,"Dosage"]).copy()
                if len(t)<10 or t["Dosage"].nunique()<2 or int((t["Dosage"]>0).sum())<MIN_C: continue
                slope,inter,r,p,se=stats.linregress(t["Dosage"],t[e]); row={"Analysis_Type":pcol,"Cohort":cname,idcol:key,"Endpoint":e,"N":len(t),"Carrier_Count":int((t["Dosage"]>0).sum()),"N_Dosage_0":int((t["Dosage"]==0).sum()),"N_Dosage_1":int((t["Dosage"]==1).sum()),"N_Dosage_2":int((t["Dosage"]==2).sum()),"Slope":slope,"Intercept":inter,"R":r,"SE":se,"P_Value":p}
                if "Carrier" in t.columns:
                    a=pd.to_numeric(t.loc[t["Carrier"]==1,e],errors="coerce").dropna(); b=pd.to_numeric(t.loc[t["Carrier"]==0,e],errors="coerce").dropna(); row["P_Carrier"]=stats.mannwhitneyu(a,b,alternative="two-sided")[1] if len(a)>=MIN_C and len(b)>=MIN_C else np.nan; row["Carrier_Median_Diff"]=a.median()-b.median() if len(a) and len(b) else np.nan
                for c in ["SNP_Label","Primary_rsID","Variant_Display","gnomAD_NFE_AF","Haplotype","Global_Freq"]:
                    if c in t.columns: row[c]=first(t[c])
                rows.append(row)
    return add_fdr(pd.DataFrame(rows),"P_Value",["Analysis_Type","Cohort"]) if rows else pd.DataFrame()

def clin_bin(df,y,x,bmi_adj):
    out={"OR_Adj_Age":np.nan,"OR_Adj_Age_CI95":np.nan,"P_Adj_Age":np.nan,"N_Adj_Age":np.nan,"Adj_Age_Note":"statsmodels not available","OR_Adj_AgeBMI":np.nan,"OR_Adj_AgeBMI_CI95":np.nan,"P_Adj_AgeBMI":np.nan,"N_Adj_AgeBMI":np.nan,"Adj_AgeBMI_Note":"not run"}
    try:
        import statsmodels.formula.api as smf
    except ImportError:
        return out
    a=df[[y,x,"canon__age"]].dropna().copy();
    if len(a)>=10 and a[y].nunique()==2:
        a["endpoint_z"]=z(a[x]); a["age_z"]=z(a["canon__age"]); a=a.dropna()
        if len(a)>=10:
            try:
                fit=smf.logit(f"{y} ~ endpoint_z + age_z",data=a).fit(disp=0,maxiter=200); ci=fit.conf_int().loc["endpoint_z"]; out.update({"OR_Adj_Age":float(np.exp(fit.params['endpoint_z'])),"OR_Adj_Age_CI95":f"[{float(np.exp(ci.iloc[0])):.3f}, {float(np.exp(ci.iloc[1])):.3f}]","P_Adj_Age":float(fit.pvalues['endpoint_z']),"N_Adj_Age":len(a),"Adj_Age_Note":"converged"})
            except Exception as e: out["Adj_Age_Note"]=str(e)[:80]
    if bmi_adj:
        b=df[[y,x,"canon__age","canon__bmi"]].dropna().copy()
        if len(b)>=10 and b[y].nunique()==2:
            b["endpoint_z"]=z(b[x]); b["age_z"]=z(b["canon__age"]); b["bmi_z"]=z(b["canon__bmi"]); b=b.dropna()
            if len(b)>=10:
                try:
                    fit=smf.logit(f"{y} ~ endpoint_z + age_z + bmi_z",data=b).fit(disp=0,maxiter=200); ci=fit.conf_int().loc["endpoint_z"]; out.update({"OR_Adj_AgeBMI":float(np.exp(fit.params['endpoint_z'])),"OR_Adj_AgeBMI_CI95":f"[{float(np.exp(ci.iloc[0])):.3f}, {float(np.exp(ci.iloc[1])):.3f}]","P_Adj_AgeBMI":float(fit.pvalues['endpoint_z']),"N_Adj_AgeBMI":len(b),"Adj_AgeBMI_Note":"converged"})
                except Exception as e: out["Adj_AgeBMI_Note"]=str(e)[:80]
    return out

def assoc_clin(expr, breast_vars=CB, endo_vars=CE):
    rows=[]
    for cname,cfg in [("Breast",breast_vars),("Endometrial",endo_vars)]:
        cdf=expr[expr["cohort"]==cname].copy()
        for e in ENDPTS:
            for cv,(typ,label,bmi_adj) in cfg.items():
                if cv not in cdf.columns: continue
                t=cdf[[e,cv,"canon__age","canon__bmi"]].copy(); t[e]=pd.to_numeric(t[e],errors="coerce"); t=t.dropna(subset=[e,cv])
                if len(t)<10: continue
                row={"Analysis_Type":"Isoform_vs_Clinical","Cohort":cname,"Endpoint":e,"Clin_Var":cv,"Clin_Label":label,"Clin_Type":typ,"N":len(t),"P_Value":np.nan,"Effect":np.nan,"Test":"","OR_Adj_Age":np.nan,"OR_Adj_Age_CI95":np.nan,"P_Adj_Age":np.nan,"N_Adj_Age":np.nan,"Adj_Age_Note":"","OR_Adj_AgeBMI":np.nan,"OR_Adj_AgeBMI_CI95":np.nan,"P_Adj_AgeBMI":np.nan,"N_Adj_AgeBMI":np.nan,"Adj_AgeBMI_Note":""}
                if typ=="continuous":
                    t[cv]=pd.to_numeric(t[cv],errors="coerce"); t=t.dropna(subset=[cv]);
                    if len(t)<10: continue
                    rho,p=stats.spearmanr(t[e],t[cv],nan_policy="omit"); row.update({"Test":"Spearman","P_Value":p,"Effect":rho})
                elif typ=="binary":
                    t[cv]=pd.to_numeric(t[cv],errors="coerce"); t=t[t[cv].isin([0,1])].copy();
                    if t.empty or t[cv].nunique()<2: continue
                    a=pd.to_numeric(t.loc[t[cv]==1,e],errors="coerce").dropna(); b=pd.to_numeric(t.loc[t[cv]==0,e],errors="coerce").dropna();
                    if len(a)<MIN_C or len(b)<MIN_C: continue
                    row.update({"Test":"Mann-Whitney U","P_Value":stats.mannwhitneyu(a,b,alternative="two-sided")[1],"Effect":a.median()-b.median()}); row.update(clin_bin(t,cv,e,bmi_adj))
                else:
                    t[cv]=t[cv].astype(str).str.strip(); groups=[pd.to_numeric(g[e],errors="coerce").dropna() for _,g in t.groupby(cv)]; groups=[g for g in groups if len(g)>=3]
                    if len(groups)<2: continue
                    st,p=stats.kruskal(*groups); row.update({"Test":"Kruskal-Wallis","P_Value":p,"Effect":st})
                rows.append(row)
    return add_fdr(pd.DataFrame(rows),"P_Value",["Analysis_Type","Cohort"]) if rows else pd.DataFrame()

def bridge(snp,hap,clin):
    rows=[]
    for df,idc in [(snp,"Variant_ID"),(hap,"Haplotype_ID")]:
        if df.empty: continue
        for (c,e),g in df[df["Nominal_Sig"]].groupby(["Cohort","Endpoint"],dropna=False): rows.append({"Cohort":c,"Endpoint":e,"Genetic_Type":idc,"Genetic_Hits":"; ".join(sorted(g[idc].astype(str).unique())),"Num_Genetic_Hits":g[idc].nunique()})
    gdf=pd.DataFrame(rows)
    if not gdf.empty: gdf=(gdf.groupby(["Cohort","Endpoint"],dropna=False).agg({"Genetic_Type":lambda s:"; ".join(sorted(set(s))),"Genetic_Hits":lambda s:"; ".join(sorted({v for t in s for v in str(t).split('; ') if v})),"Num_Genetic_Hits":"sum"}).reset_index())
    else: gdf=pd.DataFrame(columns=["Cohort","Endpoint","Genetic_Type","Genetic_Hits","Num_Genetic_Hits"])
    if clin.empty: return gdf
    cdf=(clin[clin["Nominal_Sig"]].groupby(["Cohort","Endpoint"],dropna=False).agg({"Clin_Label":lambda s:"; ".join(sorted(pd.Series(s).dropna().astype(str).unique())),"Clin_Var":"nunique"}).rename(columns={"Clin_Label":"Clinical_Hits","Clin_Var":"Num_Clinical_Hits"}).reset_index())
    out=gdf.merge(cdf,on=["Cohort","Endpoint"],how="outer"); out["Bridge_Signal"]=(out["Num_Genetic_Hits"].fillna(0)+out["Num_Clinical_Hits"].fillna(0))>=2; return out.sort_values(["Bridge_Signal","Cohort","Endpoint"],ascending=[False,True,True])

def panel_cols(df):
    exc=set(ENDPTS+["NOMBRE DE LA MUESTRA","CODIGO JC","Tissue","SNP28","rs869402","observaciones EVA","TIENEN RNA","FALTA MUESTRA","cohort","tumour_normal","analysis_group","analysis_role","raw1","raw2","sx1","sx2","snp_code","match_status","match_method","match_candidates","master_join_success","analysis_include_primary","isoform_total_sum","isoform_total_scale","isoform_total_flag","endpoint_missing_count"]+RNA_QC_COLUMNS)
    return sorted([c for c in df.columns if c not in exc and pd.to_numeric(df[c],errors="coerce").notna().sum()>=10])

def compute_immune_response_score(df):
    genes=[c for c in IMMUNE_SIGNATURE_GENES if c in df.columns]
    score=pd.Series(np.nan,index=df.index,dtype=float)
    if len(genes) < 2:
        return score, genes
    zscores=pd.DataFrame({g: z(df[g]) for g in genes})
    score=zscores.mean(axis=1, skipna=True)
    score[zscores.notna().sum(axis=1) < 2] = np.nan
    return score, genes

def exploratory(expr,gl,hc,snp,hap,genes):
    if not genes: return pd.DataFrame({"Note":["No exploratory panel-gene columns were detected."]})
    fs=[] if snp.empty else snp[(snp["Endpoint"].isin(["GSDMB","pyroptotic_isoform_fraction"]))&(snp["P_Value"]<0.05)]["Variant_ID"].dropna().astype(str).unique().tolist()
    fh=[] if hap.empty else hap[(hap["Endpoint"].isin(["GSDMB","pyroptotic_isoform_fraction"]))&(hap["P_Value"]<0.05)]["Haplotype_ID"].dropna().astype(str).unique().tolist()
    if not fs and not fh: return pd.DataFrame({"Note":["No nominally significant GSDMB-core genetic signals were available for exploratory panel-gene testing."]})
    rows=[]
    for cname in ["Breast","Endometrial"]:
        cdf=expr[expr["cohort"]==cname].copy(); codes=set(cdf["snp_code"].dropna())
        for vid in fs:
            g=gl[(gl["Variant_ID"]==vid)&(gl["snp_code"].isin(codes))][[c for c in ["analysis_sample_id","snp_code","Dosage"] if c in gl.columns]].drop_duplicates([c for c in ["analysis_sample_id"] if c in gl.columns] or None).copy(); m=cdf.merge(g,on="snp_code",how="inner")
            for gene in genes:
                t=m[[gene,"Dosage"]].copy(); t[gene]=pd.to_numeric(t[gene],errors="coerce"); t=t.dropna().copy();
                if len(t)>=10 and t["Dosage"].nunique()>=2:
                    sl,_,r,p,_=stats.linregress(t["Dosage"],t[gene]); rows.append({"Analysis_Type":"Exploratory_SNP_PanelGene","Cohort":cname,"Exposure":vid,"Target_Gene":gene,"N":len(t),"Slope":sl,"R":r,"P_Value":p})
        for hid in fh:
            g=hc[(hc["Haplotype_ID"]==hid)&(hc["snp_code"].isin(codes))][[c for c in ["analysis_sample_id","snp_code","Dosage"] if c in hc.columns]].drop_duplicates([c for c in ["analysis_sample_id"] if c in hc.columns] or None).copy(); m=cdf.merge(g,on="snp_code",how="inner")
            for gene in genes:
                t=m[[gene,"Dosage"]].copy(); t[gene]=pd.to_numeric(t[gene],errors="coerce"); t=t.dropna().copy();
                if len(t)>=10 and t["Dosage"].nunique()>=2:
                    sl,_,r,p,_=stats.linregress(t["Dosage"],t[gene]); rows.append({"Analysis_Type":"Exploratory_Haplotype_PanelGene","Cohort":cname,"Exposure":hid,"Target_Gene":gene,"N":len(t),"Slope":sl,"R":r,"P_Value":p})
    return add_fdr(pd.DataFrame(rows),"P_Value",["Analysis_Type","Cohort","Exposure"]) if rows else pd.DataFrame({"Note":["No exploratory panel-gene tests met the minimum data thresholds."]})

def bam_inventory(root,expr):
    roots=[("Breast_Tumour",root/"breast"/"tumour"/"rna"),("Breast_Normal",root/"breast"/"normal"/"rna"),("Endometrial_Tumour",root/"endometrium"/"tumour"/"rna"),("Endometrial_Normal",root/"endometrium"/"normal"/"rna")]
    ec=set(expr["snp_code"].dropna()); pc=set(expr.loc[expr["analysis_include_primary"],"snp_code"].dropna()); rows=[]
    for lab,p in roots:
        if not p.exists(): continue
        for b in sorted(p.glob("*.bam")):
            code=sx(b.stem); rows.append({"bam_group":lab,"bam_path":str(b),"bam_name":b.name,"size_mb":round(b.stat().st_size/(1024*1024),2),"has_bai":b.with_suffix(".bam.bai").exists(),"snp_code":code,"matched_expression_sample":code in ec if code else False,"matched_primary_sample":code in pc if code else False})
    inv=pd.DataFrame(rows)
    if inv.empty: return inv,pd.DataFrame({"Note":["No RNA BAM files were found under the expected project directories."]})
    summ=inv.groupby("bam_group",dropna=False).agg(n_bams=("bam_name","count"),matched_expression=("matched_expression_sample","sum"),matched_primary=("matched_primary_sample","sum"),with_index=("has_bai","sum")).reset_index(); summ["samtools_available"]=shutil.which("samtools") is not None; summ["validation_note"]=np.where(summ["samtools_available"],"Inventory complete; deeper junction/coverage QC can be added with samtools-backed follow-up.","Inventory complete; deeper junction/coverage QC deferred because samtools is not available in this runtime.")
    return inv,summ

def plot_dist(df,out,controls):
    p=df.copy()
    group_map={
        "Breast_Tumour":"Breast tumour",
        "Endometrial_Tumour":"Endometrial tumour",
        "Pooled_Control":"Pooled control",
    }
    colour_map={
        "Breast_Tumour": COHORT_COLORS.get("Breast_Tumour", "#CC79A7"),
        "Endometrial_Tumour": COHORT_COLORS.get("Endometrium_Tumour", "#009E73"),
        "Pooled_Control": TUMOUR_CONTROL_COLORS.get("Pooled control", "#0072B2"),
    }
    if controls:
        p=p[p["analysis_group"].isin(["Breast_Tumour","Endometrial_Tumour","Breast_Normal","Endometrial_Normal"])].copy()
        p["plot_group"]=p["analysis_group"].replace({
            "Breast_Tumour":"Breast_Tumour",
            "Endometrial_Tumour":"Endometrial_Tumour",
            "Breast_Normal":"Pooled_Control",
            "Endometrial_Normal":"Pooled_Control",
        })
        order=[g for g in ["Breast_Tumour","Endometrial_Tumour","Pooled_Control"] if g in p["plot_group"].unique()]
        comparisons=[
            ("Breast_Tumour","Pooled_Control","Breast tumour vs pooled control"),
            ("Endometrial_Tumour","Pooled_Control","Endometrial tumour vs pooled control"),
            ("Breast_Tumour","Endometrial_Tumour","Breast tumour vs endometrial tumour"),
        ]
    else:
        p=p[p["analysis_role"]=="primary_tumour"].copy()
        p["plot_group"]=p["analysis_group"]
        order=[g for g in ["Breast_Tumour","Endometrial_Tumour"] if g in p["plot_group"].unique()]
        comparisons=[("Breast_Tumour","Endometrial_Tumour","Breast tumour vs endometrial tumour")]
    if p.empty or not order:
        return
    endpoints=[
        "G1","G2","G3",
        "G3b","G4","GSDMB",
        "pyroptotic_isoform_fraction","non_pyroptotic_isoform_fraction","G4_vs_total",
    ]
    endpoint_titles={
        "G1":"G1",
        "G2":"G2",
        "G3":"G3",
        "G3b":"G3b",
        "G4":"G4",
        "GSDMB":"GSDMB",
        "pyroptotic_isoform_fraction":"Pyroptotic fraction",
        "non_pyroptotic_isoform_fraction":"Non-pyroptotic fraction",
        "G4_vs_total":"G4 / total",
    }
    fig,ax=plt.subplots(3,3,figsize=(18.5,13.2),constrained_layout=True)
    legend_handles=[Line2D([0],[0],marker="s",linestyle="",markersize=9,markerfacecolor=colour_map[g],markeredgecolor="none",label=group_map.get(g,g)) for g in order]
    fill_map={"Breast_Tumour":"#F3D3DC","Endometrial_Tumour":"#D5F0E8","Pooled_Control":"#D7E7F7"}
    group_counts = summarise_group_counts(p, "plot_group")
    for idx,(a,e) in enumerate(zip(ax.flat,endpoints)):
        t=p[[e,"plot_group"]].dropna().copy()
        if t.empty:
            a.set_visible(False)
            continue
        t[e]=pd.to_numeric(t[e],errors="coerce")
        t=t.dropna(subset=[e]).copy()
        if t.empty:
            a.set_visible(False)
            continue
        scale_to_percent=t[e].abs().max()<=1.5
        t["plot_value"]=t[e]*100.0 if scale_to_percent else t[e]
        y_label="Expression (%)" if scale_to_percent else "Expression"
        if sns is not None:
            sns.violinplot(data=t,x="plot_group",y="plot_value",order=order,ax=a,inner=None,cut=0,linewidth=0.9,palette=[colour_map[g] for g in order])
            sns.boxplot(data=t,x="plot_group",y="plot_value",order=order,ax=a,width=0.22,showfliers=False,boxprops={"facecolor":"white","zorder":3},medianprops={"color":"#222222","linewidth":1.2},whiskerprops={"linewidth":1.0},capprops={"linewidth":1.0})
            sns.stripplot(data=t,x="plot_group",y="plot_value",order=order,ax=a,color="#2f2f2f",size=2.5,alpha=0.5,jitter=0.18,zorder=4)
        else:
            labels=[]; groups=[]; group_cols=[]
            for lab in order:
                vals=pd.to_numeric(t.loc[t["plot_group"].astype(str)==lab,"plot_value"],errors="coerce").dropna().values
                if len(vals):
                    labels.append(group_map.get(lab, lab)); groups.append(vals); group_cols.append(colour_map[lab])
            if groups:
                bp=a.boxplot(groups, labels=labels, patch_artist=True, showfliers=False)
                for patch,lab in zip(bp["boxes"], order):
                    patch.set_facecolor(fill_map[lab])
                for i,vals in enumerate(groups, start=1):
                    a.scatter(np.full(len(vals), i), vals, color=group_cols[i-1], s=8, alpha=0.55)
        a.set_title(endpoint_titles.get(e,e),fontweight="bold",fontsize=11)
        a.set_xlabel("")
        a.set_ylabel(y_label if idx % 3 == 0 else "")
        a.set_xticklabels([f"{group_map.get(g, pretty_group(g))}\n(n={group_counts.get(g,0)})" for g in order], rotation=0, fontsize=8.2)
        a.tick_params(axis="y", labelsize=8.5)
        a.grid(axis="y", linestyle=":", alpha=0.25)
        p_lines=pairwise_mwu_pvals(t, "plot_value", "plot_group", comparisons)
        med_lines=[]
        for grp in order:
            vals=pd.to_numeric(t.loc[t["plot_group"]==grp,"plot_value"], errors="coerce").dropna()
            if len(vals):
                med_lines.append(f"{group_map.get(grp,grp)} median={vals.median():.1f}")
        if p_lines:
            note_lines = p_lines + med_lines[:3]
            a.text(0.03, 0.97, chr(10).join(note_lines), transform=a.transAxes, ha="left", va="top", fontsize=7.0, color="#333333", bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#d0d0d0", alpha=0.9))
    title="Objective 2: GSDMB isoform expression across breast tumour, endometrial tumour, and pooled controls" if controls else "Objective 2: GSDMB isoform expression across breast and endometrial tumours"
    fig.suptitle(title,fontsize=14,fontweight="bold")
    fig.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, 1.01), ncol=min(3,len(legend_handles)), frameon=False, fontsize=9)
    if controls:
        fig.text(0.995, 0.008, "Pooled control combines breast normal and endometrial normal samples. Asterisks are not shown here; pairwise P values are reported in-panel.", ha="right", va="bottom", fontsize=8, color="#555555")
    fig.savefig(out/"20_Isoform_Distributions.png",dpi=300,bbox_inches="tight")
    plt.close(fig)

def heat(df,row,col,val,out,title,row_label_col=None,top_n=20,note=None):
    if df.empty: return
    d=df.copy(); d["_score"]=-np.log10(pd.to_numeric(d[val],errors="coerce").clip(lower=1e-300)); sign=pd.to_numeric(d.get("Slope",d.get("Effect",0)),errors="coerce").fillna(0); d["_signed"]=d["_score"]*np.sign(sign)
    label_series = d[row_label_col] if row_label_col and row_label_col in d.columns else d[row]
    if "Cohort" in d.columns:
        d["_row_label"] = d["Cohort"].astype(str) + " | " + label_series.astype(str)
    else:
        d["_row_label"] = label_series.astype(str)
    top=d.groupby("_row_label")["_score"].max().sort_values(ascending=False).head(top_n).index; d=d[d["_row_label"].isin(top)].copy(); p=d.pivot_table(index="_row_label",columns=col,values="_signed",aggfunc="max")
    if p.empty: return
    ann=d.pivot_table(index="_row_label",columns=col,values=val,aggfunc="min")
    ann_fdr=d.pivot_table(index="_row_label",columns=col,values="FDR",aggfunc="min") if "FDR" in d.columns else pd.DataFrame(index=p.index, columns=p.columns)
    ann = ann.reindex(index=p.index, columns=p.columns)
    ann_fdr = ann_fdr.reindex(index=p.index, columns=p.columns)
    annot = ann.copy().astype(object)
    for rid in annot.index:
        for cid in annot.columns:
            annot.loc[rid, cid] = heat_annotation_symbol(ann.loc[rid, cid], ann_fdr.loc[rid, cid])
    fig_w=max(12, 0.72*len(p.columns)+5)
    fig_h=max(4.5, 0.42*len(p.index)+2)
    fig,a=plt.subplots(figsize=(fig_w,fig_h))
    if sns is not None:
        sns.heatmap(p,cmap="coolwarm",center=0,linewidths=0.5,linecolor="white",ax=a,annot=annot.values,fmt="",annot_kws={"fontsize":8,"fontweight":"bold"})
    else:
        im=a.imshow(p.values, aspect="auto", cmap="coolwarm")
        a.set_xticks(range(len(p.columns))); a.set_xticklabels([short_label(c, 20) for c in p.columns], rotation=35, ha="right", fontsize=8)
        a.set_yticks(range(len(p.index))); a.set_yticklabels([short_label(i, 28) for i in p.index], fontsize=8)
        fig.colorbar(im, ax=a, fraction=0.03, pad=0.02)
    a.set_xticklabels([short_label(c, 24) for c in p.columns], rotation=35, ha="right", fontsize=9)
    a.set_yticklabels([short_label(i, 40) for i in p.index], fontsize=8.5)
    a.set_ylabel("")
    a.set_title(title,fontsize=13,fontweight="bold")
    if note:
        fig.text(0.995, 0.01, note, ha="right", va="bottom", fontsize=8, color="#555555")
    fig.savefig(out,dpi=300,bbox_inches="tight")
    plt.close(fig)

def plot_top_snp_hits(snp,out,top_n=18):
    if snp.empty:
        return
    d=snp.copy()
    d=d[pd.to_numeric(d["P_Value"], errors="coerce").notna()].copy()
    if d.empty:
        return
    d["Primary_rsID"] = d.get("Primary_rsID", pd.Series(index=d.index, dtype=object))
    d["Variant_Display"] = d.get("Variant_Display", d["Variant_ID"])
    d["Plot_Label"] = d["Cohort"].astype(str) + " | " + d["Variant_Display"].astype(str) + " | " + d["Endpoint"].astype(str)
    d["MinusLog10P"] = -np.log10(pd.to_numeric(d["P_Value"], errors="coerce").clip(lower=1e-300))
    d["Effect_Display"] = pd.to_numeric(d.get("Slope"), errors="coerce")
    keep=d.sort_values(["FDR_Sig","P_Value","Carrier_Count"], ascending=[False,True,False]).head(top_n).copy()
    keep=keep.sort_values("MinusLog10P", ascending=True)
    fig,ax=plt.subplots(figsize=(13.8, max(6.5, 0.42*len(keep)+1.8)))
    sc=ax.scatter(
        keep["MinusLog10P"],
        np.arange(len(keep)),
        c=keep["Effect_Display"],
        s=45 + keep["Carrier_Count"].fillna(0)*12,
        cmap="coolwarm",
        edgecolors="#1f1f1f",
        linewidths=0.4,
        alpha=0.9,
    )
    ax.axvline(-np.log10(0.05), color="#7a7a7a", linestyle="--", linewidth=1.0)
    ax.set_yticks(np.arange(len(keep)))
    ax.set_yticklabels([short_label(v, 58) for v in keep["Plot_Label"]], fontsize=8.5)
    ax.set_xlabel("-log10(P value)")
    ax.set_title("Top SNP vs isoform associations", fontsize=13, fontweight="bold")
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    cbar=fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Slope")
    for i,(_,rowv) in enumerate(keep.iterrows()):
        star = "**" if bool(rowv.get("FDR_Sig", False)) else "*" if bool(rowv.get("Nominal_Sig", False)) else ""
        text = f"{rowv.get('Variant_Display', rowv.get('Variant_ID',''))}  n={int(rowv.get('N',0))}, carriers={int(rowv.get('Carrier_Count',0))}{star}"
        ax.text(rowv["MinusLog10P"]+0.03, i, short_label(text, 58), va="center", fontsize=7.5, color="#303030")
    fig.text(0.99, 0.01, "* nominal P<0.05, ** FDR<0.10", ha="right", va="bottom", fontsize=8, color="#555555")
    fig.savefig(out/"20_SNP_Isoform_TopHits.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

def _load_optional_tsv(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, **kwargs)


def _apply_rna_qc_gate(merged: pd.DataFrame, root: Path, bedp: Path | None, manifest_path: Path | None):
    from objective2_rna_qc_utils import (
        RNA_MANIFEST_DECISION_COLUMNS,
        build_analysis_manifest,
        summarise_analysis_manifest,
        summarise_exclusion_reasons,
    )
    if manifest_path and manifest_path.exists():
        manifest = pd.read_csv(manifest_path, sep="	")
        gate_cols = [c for c in DNA_BASELINE_COLUMNS + RNA_QC_COLUMNS + RNA_MANIFEST_DECISION_COLUMNS if c in manifest.columns]
        join_cols = ["snp_code"] + gate_cols
        gate = manifest[join_cols].copy()
        if "RNA_QC_Analysis_Ready" not in gate.columns:
            gate["RNA_QC_Analysis_Ready"] = False
        if "RNA_QC_Exploratory_Ready" not in gate.columns:
            gate["RNA_QC_Exploratory_Ready"] = gate["RNA_QC_Analysis_Ready"]
        gate = gate.sort_values(["RNA_QC_Analysis_Ready", "RNA_QC_Exploratory_Ready"], ascending=[False, False]).drop_duplicates("snp_code")
        merged = merged.drop(columns=[c for c in gate_cols if c in merged.columns], errors="ignore").merge(gate, on="snp_code", how="left")
        for col in RNA_QC_COLUMNS + RNA_MANIFEST_DECISION_COLUMNS:
            if col not in merged.columns:
                merged[col] = np.nan
        merged["RNA_QC_Analysis_Ready"] = merged["RNA_QC_Analysis_Ready"].fillna(False).astype(bool)
        merged["RNA_QC_Exploratory_Ready"] = merged["RNA_QC_Exploratory_Ready"].fillna(merged["RNA_QC_Analysis_Ready"]).astype(bool)
        qc_dir = manifest_path.parent
        sampleqc = _load_optional_tsv(qc_dir / "RNA_Sample_QC.tsv", sep="	")
        targetcov = _load_optional_tsv(qc_dir / "RNA_Target_Coverage.tsv.gz", sep="	", compression="gzip")
        bami = _load_optional_tsv(qc_dir / "RNA_BAM_Inventory.tsv", sep="	")
        bams = _load_optional_tsv(qc_dir / "RNA_BAM_Validation_Summary.tsv", sep="	")
        paneldesign = _load_optional_tsv(qc_dir / "RNA_Panel_Design.tsv", sep="	")
        qc_summary = _load_optional_tsv(qc_dir / "RNA_QC_Summary.tsv", sep="	")
        exclusion_summary = _load_optional_tsv(qc_dir / "RNA_QC_Exclusion_Reasons.tsv", sep="	")
        return merged, manifest, qc_summary, exclusion_summary, bami, targetcov, bams, sampleqc, paneldesign

    bami, targetcov, bams, sampleqc, paneldesign = bam_validation(root, merged, bedp, sx)
    merged = attach_rna_qc(merged, sampleqc)
    manifest = build_analysis_manifest(merged)
    qc_summary = summarise_analysis_manifest(manifest)
    exclusion_summary = summarise_exclusion_reasons(manifest)
    return merged, manifest, qc_summary, exclusion_summary, bami, targetcov, bams, sampleqc, paneldesign


def parse_args():
    ap=argparse.ArgumentParser(description="Objective 2 workbook-first GSDMB isoform analysis")
    ap.add_argument("--expr-xlsx",default=str(D["expr_xlsx"]))
    ap.add_argument("--master",default=str(D["master"]))
    ap.add_argument("--variant-workbook",default=str(D["variant_workbook"]))
    ap.add_argument("--haplotype-input",default=str(D["haplotype_input"]))
    ap.add_argument("--out-dir",default=str(D["out_dir"]))
    ap.add_argument("--rna-bed",default=str(D.get("rna_bed","")))
    ap.add_argument("--rna-qc-manifest",default=str(D.get("rna_qc_manifest","")))
    controls_group=ap.add_mutually_exclusive_group()
    controls_group.add_argument("--include-controls-context",dest="controls_context",action="store_true",help="Show the main figure with pooled controls included.")
    controls_group.add_argument("--no-controls-context",dest="controls_context",action="store_false",help="Restrict the main figure to tumour samples only.")
    ap.set_defaults(controls_context=True)
    ap.add_argument("--include-supplementary-clinical",action="store_true",help="Include supplementary clinical variables in the isoform association tables.")
    ap.add_argument("--strict-rna-qc",action="store_true",help="Restrict primary analyses to strict analysis-ready RNA samples only.")
    return ap.parse_args()


def main():
    a=parse_args(); exprp=Path(a.expr_xlsx); masterp=Path(a.master); varp=Path(a.variant_workbook); happ=Path(a.haplotype_input); out=Path(a.out_dir); bedp=Path(a.rna_bed) if a.rna_bed else None; manifestp=Path(a.rna_qc_manifest) if a.rna_qc_manifest else None; out.mkdir(parents=True,exist_ok=True)
    for p,l in [(exprp,"Expression workbook"),(masterp,"Harmonised master"),(varp,"Variant workbook"),(happ,"Haplotype input")]:
        if not p.exists(): raise FileNotFoundError(f"{l} not found: {p}")
    expr=load_expr(exprp); master=load_master(masterp); merged,audit=match_expr(expr,master)
    merged, rna_manifest, rna_qc_summary, rna_exclusion_summary, bami, targetcov, bams, sampleqc, paneldesign = _apply_rna_qc_gate(merged, Path(__file__).resolve().parent, bedp, manifestp)
    merged["Immune_Response_Score"], immune_genes = compute_immune_response_score(merged)
    strict_ready = merged["RNA_QC_Analysis_Ready"].fillna(False)
    exploratory_ready = merged["RNA_QC_Exploratory_Ready"].fillna(strict_ready)
    primary_gate = strict_ready if a.strict_rna_qc else exploratory_ready
    primary_gate_label = "strict analysis-ready" if a.strict_rna_qc else "exploratory-usable"
    primary=merged[merged["analysis_include_primary"] & primary_gate].copy()
    breast_clin_vars, endo_clin_vars = clinical_var_sets(a.include_supplementary_clinical)
    breast_therapy_vars, endo_therapy_vars = therapy_var_sets()
    back=load_backbone(varp); gl,gw,hd,hf,hc=load_phased(happ,back); mergedg=merged.merge(gw,on="snp_code",how="left"); summ=summarise(mergedg); genes=panel_cols(expr)
    if primary.empty:
        snp=pd.DataFrame(); hap=pd.DataFrame(); clin=pd.DataFrame()
        therapy=pd.DataFrame()
        br=pd.DataFrame({"Note":[f"No primary tumour RNA samples passed the {primary_gate_label} RNA QC gate."]})
        expl=pd.DataFrame({"Note":[f"No exploratory RNA gene analyses were run because no primary tumour samples passed the {primary_gate_label} RNA QC gate."]})
        gene_reports=pd.DataFrame({"Note":[f"No panel-gene reports were generated because no primary tumour samples passed the {primary_gate_label} RNA QC gate."]})
    else:
        snp=assoc_genetic(primary,gl,"Variant_ID","SNP_vs_Isoform"); hap=assoc_genetic(primary,hc,"Haplotype_ID","Haplotype_vs_Isoform"); clin=assoc_clin(primary,breast_clin_vars,endo_clin_vars); br=bridge(snp,hap,clin); expl=exploratory(primary,gl,hc,snp,hap,genes)
        therapy=assoc_clin(primary,breast_therapy_vars,endo_therapy_vars)
    if not snp.empty:
        if "Primary_rsID" not in snp.columns:
            snp["Primary_rsID"] = np.nan
        if "Variant_Display" not in snp.columns:
            snp["Variant_Display"] = snp["Variant_ID"]
    if not hap.empty:
        hap["Haplotype_Display"] = hap["Haplotype_ID"].astype(str) + " (" + (pd.to_numeric(hap["Global_Freq"], errors="coerce")*100).round(1).astype(str) + "%)"
    mancols=[c for c in ["snp_code","NOMBRE DE LA MUESTRA","CODIGO JC","analysis_group","analysis_role","match_status","match_method","sample_id","case_id","sheet","cohort","tumour_normal","Tissue","TIENEN RNA","FALTA MUESTRA","observaciones EVA","isoform_total_scale","isoform_total_sum","isoform_total_flag","endpoint_missing_count","Immune_Response_Score","canon__age","canon__bmi"]+DNA_BASELINE_COLUMNS+RNA_QC_COLUMNS+["RNA_QC_Analysis_Ready","RNA_QC_Exploratory_Ready","RNA_QC_Final_Status","RNA_QC_Exclusion_Reason","RNA_QC_Eligibility_Note"] if c in merged.columns]; man=merged[mancols].sort_values(["analysis_group","snp_code","CODIGO JC"])
    gene_report_dir=out/"panel_gene_reports"
    if not primary.empty:
        gene_reports=export_panel_gene_reports(merged,primary,genes,gl,hc,gene_report_dir,breast_clin_vars,endo_clin_vars,MIN_C,FDR,mancols)
    xlsx=out/"GSDMB_Objective2_Isoform_Results.xlsx"
    with pd.ExcelWriter(xlsx,engine="openpyxl") as w:
        mergedg.sort_values(["analysis_group","snp_code","CODIGO JC"]).to_excel(w,sheet_name="expression_cleaned",index=False); summ.to_excel(w,sheet_name="cohort_summary",index=False); audit.to_excel(w,sheet_name="unmatched_sample_audit",index=False); back.to_excel(w,sheet_name="snp_backbone",index=False); hf.to_excel(w,sheet_name="haplotype_frequencies",index=False); man.to_excel(w,sheet_name="sample_manifest",index=False)
        (snp if not snp.empty else pd.DataFrame({"Note":["No SNP-vs-isoform tests met the minimum thresholds."]})).to_excel(w,sheet_name="snp_isoform_assoc",index=False)
        (hap if not hap.empty else pd.DataFrame({"Note":["No haplotype-vs-isoform tests met the minimum thresholds."]})).to_excel(w,sheet_name="haplotype_isoform_assoc",index=False)
        (clin if not clin.empty else pd.DataFrame({"Note":["No isoform-vs-clinical tests met the minimum thresholds."]})).to_excel(w,sheet_name="clinical_isoform_assoc",index=False)
        (therapy if not therapy.empty else pd.DataFrame({"Note":["No isoform-vs-therapy tests met the minimum thresholds."]})).to_excel(w,sheet_name="therapy_response_assoc",index=False)
        br.to_excel(w,sheet_name="bridge_analysis",index=False); expl.to_excel(w,sheet_name="exploratory_panel_genes",index=False); gene_reports.to_excel(w,sheet_name="panel_gene_reports",index=False); paneldesign.to_excel(w,sheet_name="rna_panel_design",index=False)
        (sampleqc if not sampleqc.empty else pd.DataFrame({"Note":["No sample-level RNA QC rows were available."]})).to_excel(w,sheet_name="rna_sample_qc",index=False)
        (targetcov if not targetcov.empty else pd.DataFrame({"Note":["No target-level RNA coverage rows were available."]})).to_excel(w,sheet_name="rna_target_coverage",index=False)
        bami.to_excel(w,sheet_name="bam_inventory",index=False); bams.to_excel(w,sheet_name="bam_validation_summary",index=False)
        (rna_qc_summary if not rna_qc_summary.empty else pd.DataFrame({"Note":["No RNA QC summary rows were available."]})).to_excel(w,sheet_name="rna_qc_summary",index=False)
        (rna_exclusion_summary if not rna_exclusion_summary.empty else pd.DataFrame({"Note":["No RNA QC exclusion rows were available."]})).to_excel(w,sheet_name="rna_qc_exclusions",index=False)
    plot_dist(merged,out,a.controls_context); plot_rna_qc(sampleqc,out)
    if not snp.empty:
        heat(snp.sort_values("P_Value"),"Variant_ID","Endpoint","P_Value",out/"20_SNP_Isoform_Heatmap.png","Top SNP vs isoform associations",row_label_col="Variant_Display",top_n=24,note="Rows are analysis-set-specific (Breast, Endometrial, or All primary tumours). * nominal P<0.05, ** FDR<0.10. rsIDs are shown where available.")
        plot_top_snp_hits(snp,out)
    if not hap.empty: heat(hap.sort_values("P_Value"),"Haplotype_ID","Endpoint","P_Value",out/"20_Haplotype_Isoform_Heatmap.png","Top haplotype vs isoform associations",row_label_col="Haplotype_Display",top_n=18,note="Rows are analysis-set-specific (Breast, Endometrial, or All primary tumours). * nominal P<0.05, ** FDR<0.10.")
    if not clin.empty: heat(clin.sort_values("P_Value"),"Clin_Label","Endpoint","P_Value",out/"20_Clinical_Isoform_Heatmap.png","Top isoform vs clinical associations",top_n=18,note="Rows are analysis-set-specific (Breast or Endometrial). * nominal P<0.05, ** FDR<0.10.")
    if not therapy.empty: heat(therapy.sort_values("P_Value"),"Clin_Label","Endpoint","P_Value",out/"20_Therapy_Isoform_Heatmap.png","Top isoform vs therapy associations",top_n=18,note="Rows are therapy-specific (breast or endometrial progression treatment). * nominal P<0.05, ** FDR<0.10.")
    print("=== Objective 2 summary ==="); print(f"Workbook rows                    : {len(expr)}"); print(f"Matched RNA rows                : {int(merged['master_join_success'].sum())}"); print(f"RNA QC strict-ready rows        : {int(strict_ready.sum())}"); print(f"RNA QC exploratory-ready rows   : {int(exploratory_ready.sum())}"); print(f"RNA QC exploratory-only rows    : {int((exploratory_ready & ~strict_ready).sum())}"); print(f"RNA QC excluded rows            : {int((~exploratory_ready).sum())}"); print(f"Primary tumour samples analysed : {primary['snp_code'].nunique()} biological RNA samples ({primary_gate_label})"); print(f"Unmatched/ambiguous RNA rows    : {len(audit)}"); print(f"Common-SNP backbone             : {back['Variant_ID'].nunique()} variants"); print(f"Haplotypes retained             : {hf['Haplotype_ID'].nunique() if not hf.empty else 0}"); print(f"SNP vs isoform tests            : {len(snp)}"); print(f"Haplotype vs isoform tests      : {len(hap)}"); print(f"Isoform vs clinical tests       : {len(clin)}"); print(f"Panel-gene report files         : {len(gene_reports)}"); print(f"BAM inventory rows              : {len(bami)}"); print(f"RNA QC manifest                 : {manifestp if manifestp and manifestp.exists() else 'built inline'}"); print(f"Results workbook                : {xlsx}"); print(f"Panel-gene reports              : {gene_report_dir}"); print(f"Figures                         : {out}")

if __name__=="__main__": main()




