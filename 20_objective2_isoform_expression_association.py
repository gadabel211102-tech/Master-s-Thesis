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
import numpy as np
import pandas as pd
try:
    import seaborn as sns
except ImportError:
    sns = None
from scipy import stats
from association_runtime import script20_defaults
from pipeline_utils import combine_gnomad_nfe, extract_rsid
from objective2_rna_qc_utils import DNA_BASELINE_COLUMNS, RNA_QC_COLUMNS, attach_rna_qc, bam_validation, plot_rna_qc
from objective2_gene_report_exports import export_panel_gene_reports
warnings.filterwarnings("ignore")
D=script20_defaults(); MIN_NFE=float(D["min_nfe_af"]); MIN_HAP=float(D["min_hap_freq"]); MIN_C=int(D["min_carriers"]); FDR=float(D["fdr_threshold"])
ISO=["G1","G2","G3","G3b","G4"]; ENDPTS=ISO+["GSDMB","pyroptotic_isoform_fraction","non_pyroptotic_isoform_fraction","G4_vs_total","G2_vs_total"]
TMAP={"breast tumor":("Breast","Tumour","Breast_Tumour","primary_tumour"),"breast tumour":("Breast","Tumour","Breast_Tumour","primary_tumour"),"breast normal":("Breast","Healthy","Breast_Normal","context_control"),"endometrial cancer":("Endometrial","Tumour","Endometrial_Tumour","primary_tumour"),"endometrial normal":("Endometrial","Healthy","Endometrial_Normal","context_control"),"endometrial_normal":("Endometrial","Healthy","Endometrial_Normal","context_control"),"ovarian":("Ovarian","Tumour","Ovarian","excluded"),"ovarian organoids":("Ovarian","Tumour","Ovarian_Organoids","excluded"),"cell line":("CellLine","Unknown","Cell_Line","excluded")}
CB={"BREAST_GRADE_NUMERIC":("continuous","Tumour grade",1),"BREAST_KI67_NUMERIC":("continuous","KI67",0),"canon__her2_copies":("continuous","HER2 copies",0),"BREAST_P53_NUMERIC":("continuous","p53",0),"BREAST_ER_BIN":("binary","ER positive",0),"BREAST_PR_BIN":("binary","PR positive",0),"BREAST_RECURRENCE_DERIVED":("binary","Recurrence / progression",1),"BREAST_ANY_METASTASIS_BIN":("binary","Any metastasis",1),"BREAST_EXITUS_DERIVED":("binary","Death",1),"BREAST_HER2_SUBTYPE":("nominal","HER2 subtype",0),"BREAST_DX_TYPE":("nominal","Histological diagnosis",0)}
CE={"ENDO_FIGO_NUMERIC":("continuous","FIGO stage",1),"ENDO_GRADE_NUMERIC":("continuous","Tumour grade",1),"ENDO_RISK_ORDINAL":("continuous","Risk of recurrence",1),"ENDO_KI67_NUMERIC":("continuous","KI67",0),"ENDO_PDL1_NUMERIC":("continuous","PD-L1",0),"ENDO_CD8_NUMERIC":("continuous","CD8+ TILs",0),"ENDO_EXITUS_DISEASE_BIN":("binary","Disease-specific death",1),"ENDO_EXITUS_BIN":("binary","All-cause death",1),"ENDO_PD_BIN":("binary","Disease progression",1),"ENDO_PTEN_BIN":("binary","PTEN loss/reduced",0),"ENDO_MLH1_BIN":("binary","MLH1 loss/reduced",0),"ENDO_N_STAGE_BIN":("binary","Lymph node involvement",1),"ENDO_LVSI_BIN":("binary","LVSI",0),"ENDO_MYOINV_BIN":("binary","Myometrial invasion =50%",1),"ENDO_MSI_BIN":("binary","MSI-H",0),"ENDO_NEEC_BIN":("binary","Non-endometrioid histology",0),"ENDO_ER_BIN":("binary","ER positive",0),"ENDO_PR_BIN":("binary","PR positive",0),"ENDO_TP53_ABN_BIN":("binary","TP53 IHC abnormal",0),"ENDO_GENE_AMP_BIN":("binary","Gene amplification",0),"canon__molecular_class":("nominal","Molecular classification",0)}

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
    for pat,fmt in [(r"^(SNP_(?:AT|EN|MN|MT-T)_\d+)",lambda m:m.group(1).upper()),(r"^SNP_DNA_(EN|MN|AT)_(\d+)",lambda m:f"SNP_{m.group(1).upper()}_{m.group(2)}"),(r"^SNP_DNA_MT[-_]T_(\d+)",lambda m:f"SNP_MT-T_{m.group(1)}"),(r"^DNA_SNP_(EN|MN|AT)_(\d+)",lambda m:f"SNP_{m.group(1).upper()}_{m.group(2)}"),(r"^DNA_(AT|EN|MN)_(\d+)",lambda m:f"SNP_{m.group(1).upper()}_{m.group(2)}"),(r"^DNA_(?:SNP_)?MT[-_]T_(\d+)",lambda m:f"SNP_MT-T_{m.group(1)}"),(r"^RNA_SNP_(AT|EN|MN|MT-T)_(\d+)",lambda m:f"SNP_{m.group(1).upper()}_{m.group(2)}"),(r"^SNP_RNA_(AT|EN|MN|MT-T)_(\d+)",lambda m:f"SNP_{m.group(1).upper()}_{m.group(2)}"),(r"^RNA_(AT|EN|MN)_(\d+)",lambda m:f"SNP_{m.group(1).upper()}_{m.group(2)}"),(r"^RNA_MT[-_]T_(\d+)",lambda m:f"SNP_MT-T_{m.group(1)}")]:
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
    a=pd.read_excel(path)
    if "gnomADe_NFE_AF" in a.columns or "gnomADg_NFE_AF" in a.columns: a=combine_gnomad_nfe(a); fc="gnomAD_NFE_AF_combined"
    else:
        c=[x for x in a.columns if re.search(r"gnomad.*nfe.*af|nfe.*af",str(x),re.I)]
        if not c: raise ValueError("No gnomAD NFE AF column found in variant workbook")
        fc=c[0]
    a[fc]=pd.to_numeric(a[fc],errors="coerce"); a["POS"]=pd.to_numeric(a["POS"],errors="coerce"); a["REF"]=a["REF"].astype(str).str.strip().str.upper(); a["ALT"]=a["ALT"].astype(str).str.strip().str.upper(); a["SNP_Label"]=a["POS"].astype("Int64").astype(str)+"_"+a["REF"]+">"+a["ALT"]
    a["Variant_ID"]=a["Existing_variation"].map(extract_rsid).fillna(a["SNP_Label"]) if "Existing_variation" in a.columns else a["SNP_Label"]
    b=(a[a[fc]>MIN_NFE].dropna(subset=["POS"]).sort_values(["POS","ALT"]).drop_duplicates("SNP_Label").copy()); b["POS"]=b["POS"].astype(int)
    keep=["POS","REF","ALT","SNP_Label","Variant_ID",fc]+(["Existing_variation"] if "Existing_variation" in b.columns else [])
    return b[keep].rename(columns={fc:"gnomAD_NFE_AF"})
def load_phased(path,back):
    g=pd.read_csv(path,sep="\t",dtype=str); g["POS"]=pd.to_numeric(g["POS"],errors="coerce"); g["REF"]=g["REF"].astype(str).str.strip().str.upper(); g["ALT"]=g["ALT"].astype(str).str.strip().str.upper(); g["SNP_Label"]=g["POS"].astype("Int64").astype(str)+"_"+g["REF"]+">"+g["ALT"]
    g=g[g["SNP_Label"].isin(set(back["SNP_Label"]))].merge(back[["SNP_Label","Variant_ID","gnomAD_NFE_AF"]],on="SNP_Label",how="left").sort_values("POS").reset_index(drop=True)
    if g.empty: raise ValueError("No overlap between backbone SNPs and phased table")
    sc=[c for c in g.columns if c not in {"CHROM","POS","ID","REF","ALT","SNP_Label","Variant_ID","gnomAD_NFE_AF"}]; rows=[]; h=[]
    for c in sc:
        code=sx(c)
        if not code: continue
        h1=[]; h2=[]
        for _,r in g.iterrows():
            a,b=gt_haps(r[c]); h1.append(a); h2.append(b); rows.append({"Sample_phased":c,"snp_code":code,"Variant_ID":r["Variant_ID"],"SNP_Label":r["SNP_Label"],"POS":r["POS"],"Dosage":gt_dose(r[c]),"GT":r[c],"gnomAD_NFE_AF":r["gnomAD_NFE_AF"]})
        h.append({"Sample_phased":c,"snp_code":code,"hap1":"".join(h1),"hap2":"".join(h2),"callable_haplotype":int("N" not in h1 and "N" not in h2)})
    gl=pd.DataFrame(rows).drop_duplicates(["snp_code","Variant_ID"]); gw=gl.pivot_table(index="snp_code",columns="Variant_ID",values="Dosage",aggfunc="first").reset_index(); hd=pd.DataFrame(h)
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
            hc.append({"snp_code":r["snp_code"],"Sample_phased":r["Sample_phased"],"Haplotype_ID":fr["Haplotype_ID"],"Haplotype":fr["Haplotype"],"Global_Freq":fr["Global_Freq"],"Dosage":dose,"Carrier":np.nan if not ok else int(dose>=1),"Callable":ok})
    return gl,gw,hd,hf,pd.DataFrame(hc)
def summarise(df):
    rows=[]
    for g,sub in df.groupby("analysis_group",dropna=False):
        row={"analysis_group":g,"analysis_role":first(sub["analysis_role"]),"n_samples":sub["snp_code"].nunique(dropna=True),"n_rows":len(sub),"n_matched_master":int(sub["master_join_success"].sum()),"n_flagged_isoform_total":int(sub["isoform_total_flag"].sum())}
        for e in ENDPTS:
            v=pd.to_numeric(sub[e],errors="coerce"); row[f"{e}__nonmissing"]=int(v.notna().sum()); row[f"{e}__median"]=v.median(skipna=True); row[f"{e}__mean"]=v.mean(skipna=True)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("analysis_group")

def assoc_genetic(expr,longdf,idcol,pcol):
    rows=[]; cohorts={"Breast":expr[expr["cohort"]=="Breast"].copy(),"Endometrial":expr[expr["cohort"]=="Endometrial"].copy(),"All_Primary_Tumours":expr.copy()}
    for cname,cdf in cohorts.items():
        codes=set(cdf["snp_code"].dropna()); sdf=longdf[longdf["snp_code"].isin(codes)].copy()
        for key,g in sdf.groupby(idcol,dropna=False):
            meta=[c for c in [idcol,"snp_code","Dosage","Carrier","SNP_Label","gnomAD_NFE_AF","Haplotype","Global_Freq"] if c in g.columns]
            m=cdf.merge(g[meta],on="snp_code",how="inner")
            if m.empty: continue
            for e in ENDPTS:
                t=m[[c for c in [e,"Dosage","Carrier","SNP_Label","gnomAD_NFE_AF","Haplotype","Global_Freq"] if c in m.columns]].dropna(subset=[e,"Dosage"]).copy()
                if len(t)<10 or t["Dosage"].nunique()<2 or int((t["Dosage"]>0).sum())<MIN_C: continue
                slope,inter,r,p,se=stats.linregress(t["Dosage"],t[e]); row={"Analysis_Type":pcol,"Cohort":cname,idcol:key,"Endpoint":e,"N":len(t),"Carrier_Count":int((t["Dosage"]>0).sum()),"N_Dosage_0":int((t["Dosage"]==0).sum()),"N_Dosage_1":int((t["Dosage"]==1).sum()),"N_Dosage_2":int((t["Dosage"]==2).sum()),"Slope":slope,"Intercept":inter,"R":r,"SE":se,"P_Value":p}
                if "Carrier" in t.columns:
                    a=pd.to_numeric(t.loc[t["Carrier"]==1,e],errors="coerce").dropna(); b=pd.to_numeric(t.loc[t["Carrier"]==0,e],errors="coerce").dropna(); row["P_Carrier"]=stats.mannwhitneyu(a,b,alternative="two-sided")[1] if len(a)>=MIN_C and len(b)>=MIN_C else np.nan; row["Carrier_Median_Diff"]=a.median()-b.median() if len(a) and len(b) else np.nan
                for c in ["SNP_Label","gnomAD_NFE_AF","Haplotype","Global_Freq"]:
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

def assoc_clin(expr):
    rows=[]
    for cname,cfg in [("Breast",CB),("Endometrial",CE)]:
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

def exploratory(expr,gl,hc,snp,hap,genes):
    if not genes: return pd.DataFrame({"Note":["No exploratory panel-gene columns were detected."]})
    fs=[] if snp.empty else snp[(snp["Endpoint"].isin(["GSDMB","pyroptotic_isoform_fraction"]))&(snp["P_Value"]<0.05)]["Variant_ID"].dropna().astype(str).unique().tolist()
    fh=[] if hap.empty else hap[(hap["Endpoint"].isin(["GSDMB","pyroptotic_isoform_fraction"]))&(hap["P_Value"]<0.05)]["Haplotype_ID"].dropna().astype(str).unique().tolist()
    if not fs and not fh: return pd.DataFrame({"Note":["No nominally significant GSDMB-core genetic signals were available for exploratory panel-gene testing."]})
    rows=[]
    for cname in ["Breast","Endometrial"]:
        cdf=expr[expr["cohort"]==cname].copy(); codes=set(cdf["snp_code"].dropna())
        for vid in fs:
            g=gl[(gl["Variant_ID"]==vid)&(gl["snp_code"].isin(codes))][["snp_code","Dosage"]].copy(); m=cdf.merge(g,on="snp_code",how="inner")
            for gene in genes:
                t=m[[gene,"Dosage"]].copy(); t[gene]=pd.to_numeric(t[gene],errors="coerce"); t=t.dropna().copy();
                if len(t)>=10 and t["Dosage"].nunique()>=2:
                    sl,_,r,p,_=stats.linregress(t["Dosage"],t[gene]); rows.append({"Analysis_Type":"Exploratory_SNP_PanelGene","Cohort":cname,"Exposure":vid,"Target_Gene":gene,"N":len(t),"Slope":sl,"R":r,"P_Value":p})
        for hid in fh:
            g=hc[(hc["Haplotype_ID"]==hid)&(hc["snp_code"].isin(codes))][["snp_code","Dosage"]].copy(); m=cdf.merge(g,on="snp_code",how="inner")
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
    if controls:
        p=p[p["analysis_group"].isin(["Breast_Tumour","Endometrial_Tumour","Breast_Normal","Endometrial_Normal"])].copy()
        p["plot_group"]=p["analysis_group"].replace({
            "Breast_Normal":"Pooled_Control",
            "Endometrial_Normal":"Pooled_Control",
        })
        order=[g for g in ["Breast_Tumour","Endometrial_Tumour","Pooled_Control"] if g in p["plot_group"].unique()]
    else:
        p=p[p["analysis_role"]=="primary_tumour"].copy()
        p["plot_group"]=p["analysis_group"]
        order=[g for g in ["Breast_Tumour","Endometrial_Tumour"] if g in p["plot_group"].unique()]
    if p.empty or not order: return
    fig,ax=plt.subplots(2,3,figsize=(16.5,8.8),constrained_layout=True)
    display_map={
        "Breast_Tumour":"Breast\nTumour",
        "Endometrial_Tumour":"Endometrial\nTumour",
        "Pooled_Control":"Pooled\nControl",
    }
    display_order=[display_map.get(g, pretty_group(g)) for g in order]
    pair_labels=[
        ("Breast_Tumour","Pooled_Control","Breast tum. vs pooled ctrl"),
        ("Endometrial_Tumour","Pooled_Control","Endom. tum. vs pooled ctrl"),
        ("Breast_Tumour","Endometrial_Tumour","Breast vs endom. tum."),
    ]
    for a,e in zip(ax.flat,["G1","G2","G3","G3b","G4","GSDMB"]):
        t=p[[e,"plot_group"]].dropna().copy()
        if t.empty: a.set_visible(False); continue
        if sns is not None:
            sns.boxplot(data=t,x="plot_group",y=e,order=order,ax=a,color="#d9e8f5",fliersize=0,width=0.62)
            sns.stripplot(data=t,x="plot_group",y=e,order=order,ax=a,color="#2f5d80",size=2.8,alpha=0.6,jitter=0.22)
        else:
            labels=[]; groups=[]
            for lab in order:
                vals=pd.to_numeric(t.loc[t["plot_group"].astype(str)==lab,e],errors="coerce").dropna().values
                if len(vals):
                    labels.append(lab); groups.append(vals)
            if groups:
                a.boxplot(groups, labels=labels, patch_artist=True, boxprops=dict(facecolor="#d9e8f5"))
                for i,vals in enumerate(groups, start=1):
                    a.scatter(np.full(len(vals), i), vals, color="#2f5d80", s=8, alpha=0.6)
        a.set_title(e,fontweight="bold",fontsize=11)
        a.set_xlabel("")
        a.set_xticklabels(display_order, rotation=0, fontsize=8.5)
        a.tick_params(axis="y", labelsize=8.5)
        a.grid(axis="y", linestyle=":", alpha=0.25)
        p_lines=pairwise_mwu_pvals(t, e, "plot_group", pair_labels if controls else [pair_labels[-1]])
        if p_lines:
            a.text(
                0.03, 0.97, "\n".join(p_lines),
                transform=a.transAxes, ha="left", va="top",
                fontsize=7.3, color="#333333",
                bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#d0d0d0", alpha=0.88),
            )
    title="Objective 2: GSDMB isoform and total-expression distributions"
    if controls:
        title += " (with pooled controls)"
    fig.suptitle(title,fontsize=14,fontweight="bold")
    fig.savefig(out/"20_Isoform_Distributions.png",dpi=300,bbox_inches="tight")
    plt.close(fig)

def heat(df,row,col,val,out,title):
    if df.empty: return
    d=df.copy(); d["_score"]=-np.log10(pd.to_numeric(d[val],errors="coerce").clip(lower=1e-300)); sign=pd.to_numeric(d.get("Slope",d.get("Effect",0)),errors="coerce").fillna(0); d["_signed"]=d["_score"]*np.sign(sign)
    top=d.groupby(row)["_score"].max().sort_values(ascending=False).head(20).index; d=d[d[row].isin(top)].copy(); p=d.pivot_table(index=row,columns=col,values="_signed",aggfunc="max")
    if p.empty: return
    fig_w=max(12, 0.72*len(p.columns)+5)
    fig_h=max(4.5, 0.42*len(p.index)+2)
    fig,a=plt.subplots(figsize=(fig_w,fig_h))
    if sns is not None:
        sns.heatmap(p,cmap="coolwarm",center=0,linewidths=0.5,linecolor="white",ax=a)
    else:
        im=a.imshow(p.values, aspect="auto", cmap="coolwarm")
        a.set_xticks(range(len(p.columns))); a.set_xticklabels([short_label(c, 20) for c in p.columns], rotation=35, ha="right", fontsize=8)
        a.set_yticks(range(len(p.index))); a.set_yticklabels([short_label(i, 28) for i in p.index], fontsize=8)
        fig.colorbar(im, ax=a, fraction=0.03, pad=0.02)
    a.set_xticklabels([short_label(c, 20) for c in p.columns], rotation=35, ha="right", fontsize=8)
    a.set_yticklabels([short_label(i, 28) for i in p.index], fontsize=8)
    a.set_title(title,fontsize=13,fontweight="bold")
    fig.savefig(out,dpi=300,bbox_inches="tight")
    plt.close(fig)

def parse_args():
    ap=argparse.ArgumentParser(description="Objective 2 workbook-first GSDMB isoform analysis"); ap.add_argument("--expr-xlsx",default=str(D["expr_xlsx"])); ap.add_argument("--master",default=str(D["master"])); ap.add_argument("--variant-workbook",default=str(D["variant_workbook"])); ap.add_argument("--haplotype-input",default=str(D["haplotype_input"])); ap.add_argument("--out-dir",default=str(D["out_dir"])); ap.add_argument("--rna-bed",default=str(D.get("rna_bed",""))); ap.add_argument("--include-controls-context",action="store_true"); return ap.parse_args()

def main():
    a=parse_args(); exprp=Path(a.expr_xlsx); masterp=Path(a.master); varp=Path(a.variant_workbook); happ=Path(a.haplotype_input); out=Path(a.out_dir); bedp=Path(a.rna_bed) if a.rna_bed else None; out.mkdir(parents=True,exist_ok=True)
    for p,l in [(exprp,"Expression workbook"),(masterp,"Harmonised master"),(varp,"Variant workbook"),(happ,"Haplotype input")]:
        if not p.exists(): raise FileNotFoundError(f"{l} not found: {p}")
    expr=load_expr(exprp); master=load_master(masterp); merged,audit=match_expr(expr,master); bami,targetcov,bams,sampleqc,paneldesign=bam_validation(Path(__file__).resolve().parent,merged,bedp,sx); merged=attach_rna_qc(merged,sampleqc); primary=merged[merged["analysis_include_primary"]].copy()
    if primary.empty: raise ValueError("No primary tumour RNA samples could be matched to the harmonised master.")
    back=load_backbone(varp); gl,gw,hd,hf,hc=load_phased(happ,back); mergedg=merged.merge(gw,on="snp_code",how="left"); summ=summarise(mergedg); snp=assoc_genetic(primary,gl,"Variant_ID","SNP_vs_Isoform"); hap=assoc_genetic(primary,hc,"Haplotype_ID","Haplotype_vs_Isoform"); clin=assoc_clin(primary); br=bridge(snp,hap,clin); genes=panel_cols(expr); expl=exploratory(primary,gl,hc,snp,hap,genes)
    mancols=[c for c in ["snp_code","NOMBRE DE LA MUESTRA","CODIGO JC","analysis_group","analysis_role","match_status","match_method","sample_id","case_id","sheet","cohort","tumour_normal","Tissue","TIENEN RNA","FALTA MUESTRA","observaciones EVA","isoform_total_scale","isoform_total_sum","isoform_total_flag","endpoint_missing_count","canon__age","canon__bmi"]+RNA_QC_COLUMNS if c in merged.columns]; man=merged[mancols].sort_values(["analysis_group","snp_code","CODIGO JC"])
    gene_report_dir=out/"panel_gene_reports"; gene_reports=export_panel_gene_reports(merged,primary,genes,gl,hc,gene_report_dir,CB,CE,MIN_C,FDR,mancols)
    xlsx=out/"GSDMB_Objective2_Isoform_Results.xlsx"
    with pd.ExcelWriter(xlsx,engine="openpyxl") as w:
        mergedg.sort_values(["analysis_group","snp_code","CODIGO JC"]).to_excel(w,sheet_name="expression_cleaned",index=False); summ.to_excel(w,sheet_name="cohort_summary",index=False); audit.to_excel(w,sheet_name="unmatched_sample_audit",index=False); back.to_excel(w,sheet_name="snp_backbone",index=False); hf.to_excel(w,sheet_name="haplotype_frequencies",index=False); man.to_excel(w,sheet_name="sample_manifest",index=False)
        (snp if not snp.empty else pd.DataFrame({"Note":["No SNP-vs-isoform tests met the minimum thresholds."]})).to_excel(w,sheet_name="snp_isoform_assoc",index=False)
        (hap if not hap.empty else pd.DataFrame({"Note":["No haplotype-vs-isoform tests met the minimum thresholds."]})).to_excel(w,sheet_name="haplotype_isoform_assoc",index=False)
        (clin if not clin.empty else pd.DataFrame({"Note":["No isoform-vs-clinical tests met the minimum thresholds."]})).to_excel(w,sheet_name="clinical_isoform_assoc",index=False)
        br.to_excel(w,sheet_name="bridge_analysis",index=False); expl.to_excel(w,sheet_name="exploratory_panel_genes",index=False); gene_reports.to_excel(w,sheet_name="panel_gene_reports",index=False); paneldesign.to_excel(w,sheet_name="rna_panel_design",index=False); (sampleqc if not sampleqc.empty else pd.DataFrame({"Note":["No sample-level RNA QC rows were available."]})).to_excel(w,sheet_name="rna_sample_qc",index=False); (targetcov if not targetcov.empty else pd.DataFrame({"Note":["No target-level RNA coverage rows were available."]})).to_excel(w,sheet_name="rna_target_coverage",index=False); bami.to_excel(w,sheet_name="bam_inventory",index=False); bams.to_excel(w,sheet_name="bam_validation_summary",index=False)
    plot_dist(merged,out,bool(a.include_controls_context)); plot_rna_qc(sampleqc,out)
    if not snp.empty: heat(snp.sort_values("P_Value"),"Variant_ID","Endpoint","P_Value",out/"20_SNP_Isoform_Heatmap.png","Top SNP vs isoform associations")
    if not hap.empty: heat(hap.sort_values("P_Value"),"Haplotype_ID","Endpoint","P_Value",out/"20_Haplotype_Isoform_Heatmap.png","Top haplotype vs isoform associations")
    if not clin.empty: heat(clin.sort_values("P_Value"),"Clin_Label","Endpoint","P_Value",out/"20_Clinical_Isoform_Heatmap.png","Top isoform vs clinical associations")
    print("=== Objective 2 summary ==="); print(f"Workbook rows                    : {len(expr)}"); print(f"Matched RNA rows                : {int(merged['master_join_success'].sum())}"); print(f"Primary tumour samples analysed : {primary['snp_code'].nunique()}"); print(f"Unmatched/ambiguous RNA rows    : {len(audit)}"); print(f"Common-SNP backbone             : {back['Variant_ID'].nunique()} variants"); print(f"Haplotypes retained             : {hf['Haplotype_ID'].nunique() if not hf.empty else 0}"); print(f"SNP vs isoform tests            : {len(snp)}"); print(f"Haplotype vs isoform tests      : {len(hap)}"); print(f"Isoform vs clinical tests       : {len(clin)}"); print(f"Panel-gene report files         : {len(gene_reports)}"); print(f"BAM inventory rows              : {len(bami)}"); print(f"Results workbook                : {xlsx}"); print(f"Panel-gene reports              : {gene_report_dir}"); print(f"Figures                         : {out}")


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
        manifest = pd.read_csv(manifest_path, sep="\t")
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
        sampleqc = _load_optional_tsv(qc_dir / "RNA_Sample_QC.tsv", sep="\t")
        targetcov = _load_optional_tsv(qc_dir / "RNA_Target_Coverage.tsv.gz", sep="\t", compression="gzip")
        bami = _load_optional_tsv(qc_dir / "RNA_BAM_Inventory.tsv", sep="\t")
        bams = _load_optional_tsv(qc_dir / "RNA_BAM_Validation_Summary.tsv", sep="\t")
        paneldesign = _load_optional_tsv(qc_dir / "RNA_Panel_Design.tsv", sep="\t")
        qc_summary = _load_optional_tsv(qc_dir / "RNA_QC_Summary.tsv", sep="\t")
        exclusion_summary = _load_optional_tsv(qc_dir / "RNA_QC_Exclusion_Reasons.tsv", sep="\t")
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
    ap.add_argument("--include-controls-context",action="store_true")
    ap.add_argument("--strict-rna-qc",action="store_true",help="Restrict primary analyses to strict analysis-ready RNA samples only.")
    return ap.parse_args()


def main():
    a=parse_args(); exprp=Path(a.expr_xlsx); masterp=Path(a.master); varp=Path(a.variant_workbook); happ=Path(a.haplotype_input); out=Path(a.out_dir); bedp=Path(a.rna_bed) if a.rna_bed else None; manifestp=Path(a.rna_qc_manifest) if a.rna_qc_manifest else None; out.mkdir(parents=True,exist_ok=True)
    for p,l in [(exprp,"Expression workbook"),(masterp,"Harmonised master"),(varp,"Variant workbook"),(happ,"Haplotype input")]:
        if not p.exists(): raise FileNotFoundError(f"{l} not found: {p}")
    expr=load_expr(exprp); master=load_master(masterp); merged,audit=match_expr(expr,master)
    merged, rna_manifest, rna_qc_summary, rna_exclusion_summary, bami, targetcov, bams, sampleqc, paneldesign = _apply_rna_qc_gate(merged, Path(__file__).resolve().parent, bedp, manifestp)
    strict_ready = merged["RNA_QC_Analysis_Ready"].fillna(False)
    exploratory_ready = merged["RNA_QC_Exploratory_Ready"].fillna(strict_ready)
    primary_gate = strict_ready if a.strict_rna_qc else exploratory_ready
    primary_gate_label = "strict analysis-ready" if a.strict_rna_qc else "exploratory-usable"
    primary=merged[merged["analysis_include_primary"] & primary_gate].copy()
    back=load_backbone(varp); gl,gw,hd,hf,hc=load_phased(happ,back); mergedg=merged.merge(gw,on="snp_code",how="left"); summ=summarise(mergedg); genes=panel_cols(expr)
    if primary.empty:
        snp=pd.DataFrame(); hap=pd.DataFrame(); clin=pd.DataFrame()
        br=pd.DataFrame({"Note":[f"No primary tumour RNA samples passed the {primary_gate_label} RNA QC gate."]})
        expl=pd.DataFrame({"Note":[f"No exploratory RNA gene analyses were run because no primary tumour samples passed the {primary_gate_label} RNA QC gate."]})
        gene_reports=pd.DataFrame({"Note":[f"No panel-gene reports were generated because no primary tumour samples passed the {primary_gate_label} RNA QC gate."]})
    else:
        snp=assoc_genetic(primary,gl,"Variant_ID","SNP_vs_Isoform"); hap=assoc_genetic(primary,hc,"Haplotype_ID","Haplotype_vs_Isoform"); clin=assoc_clin(primary); br=bridge(snp,hap,clin); expl=exploratory(primary,gl,hc,snp,hap,genes)
    mancols=[c for c in ["snp_code","NOMBRE DE LA MUESTRA","CODIGO JC","analysis_group","analysis_role","match_status","match_method","sample_id","case_id","sheet","cohort","tumour_normal","Tissue","TIENEN RNA","FALTA MUESTRA","observaciones EVA","isoform_total_scale","isoform_total_sum","isoform_total_flag","endpoint_missing_count","canon__age","canon__bmi"]+DNA_BASELINE_COLUMNS+RNA_QC_COLUMNS+["RNA_QC_Analysis_Ready","RNA_QC_Exploratory_Ready","RNA_QC_Final_Status","RNA_QC_Exclusion_Reason","RNA_QC_Eligibility_Note"] if c in merged.columns]; man=merged[mancols].sort_values(["analysis_group","snp_code","CODIGO JC"])
    gene_report_dir=out/"panel_gene_reports"
    if not primary.empty:
        gene_reports=export_panel_gene_reports(merged,primary,genes,gl,hc,gene_report_dir,CB,CE,MIN_C,FDR,mancols)
    xlsx=out/"GSDMB_Objective2_Isoform_Results.xlsx"
    with pd.ExcelWriter(xlsx,engine="openpyxl") as w:
        mergedg.sort_values(["analysis_group","snp_code","CODIGO JC"]).to_excel(w,sheet_name="expression_cleaned",index=False); summ.to_excel(w,sheet_name="cohort_summary",index=False); audit.to_excel(w,sheet_name="unmatched_sample_audit",index=False); back.to_excel(w,sheet_name="snp_backbone",index=False); hf.to_excel(w,sheet_name="haplotype_frequencies",index=False); man.to_excel(w,sheet_name="sample_manifest",index=False)
        (snp if not snp.empty else pd.DataFrame({"Note":["No SNP-vs-isoform tests met the minimum thresholds."]})).to_excel(w,sheet_name="snp_isoform_assoc",index=False)
        (hap if not hap.empty else pd.DataFrame({"Note":["No haplotype-vs-isoform tests met the minimum thresholds."]})).to_excel(w,sheet_name="haplotype_isoform_assoc",index=False)
        (clin if not clin.empty else pd.DataFrame({"Note":["No isoform-vs-clinical tests met the minimum thresholds."]})).to_excel(w,sheet_name="clinical_isoform_assoc",index=False)
        br.to_excel(w,sheet_name="bridge_analysis",index=False); expl.to_excel(w,sheet_name="exploratory_panel_genes",index=False); gene_reports.to_excel(w,sheet_name="panel_gene_reports",index=False); paneldesign.to_excel(w,sheet_name="rna_panel_design",index=False)
        (sampleqc if not sampleqc.empty else pd.DataFrame({"Note":["No sample-level RNA QC rows were available."]})).to_excel(w,sheet_name="rna_sample_qc",index=False)
        (targetcov if not targetcov.empty else pd.DataFrame({"Note":["No target-level RNA coverage rows were available."]})).to_excel(w,sheet_name="rna_target_coverage",index=False)
        bami.to_excel(w,sheet_name="bam_inventory",index=False); bams.to_excel(w,sheet_name="bam_validation_summary",index=False)
        (rna_qc_summary if not rna_qc_summary.empty else pd.DataFrame({"Note":["No RNA QC summary rows were available."]})).to_excel(w,sheet_name="rna_qc_summary",index=False)
        (rna_exclusion_summary if not rna_exclusion_summary.empty else pd.DataFrame({"Note":["No RNA QC exclusion rows were available."]})).to_excel(w,sheet_name="rna_qc_exclusions",index=False)
    plot_dist(merged,out,bool(a.include_controls_context)); plot_rna_qc(sampleqc,out)
    if not snp.empty: heat(snp.sort_values("P_Value"),"Variant_ID","Endpoint","P_Value",out/"20_SNP_Isoform_Heatmap.png","Top SNP vs isoform associations")
    if not hap.empty: heat(hap.sort_values("P_Value"),"Haplotype_ID","Endpoint","P_Value",out/"20_Haplotype_Isoform_Heatmap.png","Top haplotype vs isoform associations")
    if not clin.empty: heat(clin.sort_values("P_Value"),"Clin_Label","Endpoint","P_Value",out/"20_Clinical_Isoform_Heatmap.png","Top isoform vs clinical associations")
    print("=== Objective 2 summary ==="); print(f"Workbook rows                    : {len(expr)}"); print(f"Matched RNA rows                : {int(merged['master_join_success'].sum())}"); print(f"RNA QC strict-ready rows        : {int(strict_ready.sum())}"); print(f"RNA QC exploratory-ready rows   : {int(exploratory_ready.sum())}"); print(f"RNA QC exploratory-only rows    : {int((exploratory_ready & ~strict_ready).sum())}"); print(f"RNA QC excluded rows            : {int((~exploratory_ready).sum())}"); print(f"Primary tumour samples analysed : {primary['snp_code'].nunique()} ({primary_gate_label})"); print(f"Unmatched/ambiguous RNA rows    : {len(audit)}"); print(f"Common-SNP backbone             : {back['Variant_ID'].nunique()} variants"); print(f"Haplotypes retained             : {hf['Haplotype_ID'].nunique() if not hf.empty else 0}"); print(f"SNP vs isoform tests            : {len(snp)}"); print(f"Haplotype vs isoform tests      : {len(hap)}"); print(f"Isoform vs clinical tests       : {len(clin)}"); print(f"Panel-gene report files         : {len(gene_reports)}"); print(f"BAM inventory rows              : {len(bami)}"); print(f"RNA QC manifest                 : {manifestp if manifestp and manifestp.exists() else 'built inline'}"); print(f"Results workbook                : {xlsx}"); print(f"Panel-gene reports              : {gene_report_dir}"); print(f"Figures                         : {out}")

if __name__=="__main__": main()














