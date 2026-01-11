"""
Exercise 10.2 — Identify top SNP–Gene correlations

Scop:
- incarcati matricele SNP si Expression
- calculati corelatii intre fiecare SNP si fiecare gena
- filtrati |r| > 0.5
- exportati snp_gene_pairs_<handle>.csv
"""

from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
from scipy import stats

HANDLE = "rosestoica"

SNP_CSV = Path(f"data/work/{HANDLE}/lab10/snp_matrix_{HANDLE}.csv")
EXP_CSV = Path(f"data/work/{HANDLE}/lab10/expression_matrix_{HANDLE}.csv")

OUT_DIR = Path(f"labs/10_integrative/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_CSV = OUT_DIR / f"snp_gene_pairs_{HANDLE}.csv"
OUT_HEATMAP = OUT_DIR / f"correlation_heatmap_{HANDLE}.png"

CORRELATION_THRESHOLD = 0.3  # folosim 0.3 pentru a avea mai multe rezultate pe date simulate

def load_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Incarca matricele SNP si Expression.
    """
    df_snp = pd.read_csv(SNP_CSV, index_col=0)
    df_exp = pd.read_csv(EXP_CSV, index_col=0)
    
    common_samples = df_snp.columns.intersection(df_exp.columns)
    df_snp = df_snp[common_samples]
    df_exp = df_exp[common_samples]
    
    print(f" SNP matrix: {df_snp.shape}")
    print(f" Expression matrix: {df_exp.shape}")
    print(f" Common samples: {len(common_samples)}")
    
    return df_snp, df_exp


def compute_cross_correlations(
    df_snp: pd.DataFrame,
    df_exp: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculeaza corelatia Pearson intre fiecare SNP si fiecare gena.
    Returneaza un DataFrame cu toate perechile.
    """
    results = []
    
    snp_ids = df_snp.index.tolist()
    gene_ids = df_exp.index.tolist()
    
    total_pairs = len(snp_ids) * len(gene_ids)
    print(f" Calculare corelatii pentru {total_pairs} perechi SNP-gene...")
    
    for snp_id in snp_ids:
        snp_values = df_snp.loc[snp_id].values.astype(float)
        
        for gene_id in gene_ids:
            gene_values = df_exp.loc[gene_id].values.astype(float)
            
            r, p_value = stats.pearsonr(snp_values, gene_values)
            
            results.append({
                "snp": snp_id,
                "gene": gene_id,
                "correlation": r,
                "abs_correlation": abs(r),
                "p_value": p_value
            })
    
    df_results = pd.DataFrame(results)
    return df_results


def filter_significant_pairs(
    df_corr: pd.DataFrame,
    threshold: float = 0.5
) -> pd.DataFrame:
    """
    Filtreaza perechile cu |r| >= threshold.
    """
    df_filtered = df_corr[df_corr["abs_correlation"] >= threshold].copy()
    df_filtered = df_filtered.sort_values("abs_correlation", ascending=False)
    df_filtered = df_filtered.reset_index(drop=True)
    
    return df_filtered


def plot_top_correlations(
    df_filtered: pd.DataFrame,
    df_snp: pd.DataFrame,
    df_exp: pd.DataFrame,
    out_path: Path,
    top_n: int = 20
) -> None:
    """
    Genereaza un heatmap al top corelatiilor.
    """
    import matplotlib.pyplot as plt
    
    top_pairs = df_filtered.head(top_n)
    
    if len(top_pairs) == 0:
        print("[WARN] Nu exista perechi semnificative pentru heatmap.")
        return

    top_snps = top_pairs["snp"].unique()[:10]
    top_genes = top_pairs["gene"].unique()[:10]

    corr_matrix = np.zeros((len(top_snps), len(top_genes)))
    
    for i, snp in enumerate(top_snps):
        for j, gene in enumerate(top_genes):
            snp_vals = df_snp.loc[snp].values.astype(float)
            gene_vals = df_exp.loc[gene].values.astype(float)
            r, _ = stats.pearsonr(snp_vals, gene_vals)
            corr_matrix[i, j] = r
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    im = ax.imshow(corr_matrix, cmap="RdBu_r", aspect="auto", vmin=-1, vmax=1)
    
    ax.set_xticks(range(len(top_genes)))
    ax.set_yticks(range(len(top_snps)))
    ax.set_xticklabels(top_genes, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(top_snps, fontsize=9)
    
    ax.set_xlabel("Genes", fontsize=12)
    ax.set_ylabel("SNPs", fontsize=12)
    ax.set_title("Top SNP-Gene Correlations Heatmap", fontsize=14, fontweight="bold")
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Pearson Correlation", fontsize=11)
    
    for i in range(len(top_snps)):
        for j in range(len(top_genes)):
            text = ax.text(j, i, f"{corr_matrix[i, j]:.2f}",
                          ha="center", va="center", color="black", fontsize=8)
    
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    
    print(f"[OK] Heatmap salvat: {out_path}")

if __name__ == "__main__":
    df_snp, df_exp = load_data()
    
    df_corr = compute_cross_correlations(df_snp, df_exp)
    print(f" Total perechi calculate: {len(df_corr)}")
    
    print(f"\n Statistici corelatii:")
    print(f"    Mean |r|: {df_corr['abs_correlation'].mean():.3f}")
    print(f"    Max |r|:  {df_corr['abs_correlation'].max():.3f}")
    print(f"    Min |r|:  {df_corr['abs_correlation'].min():.3f}")
    
    df_filtered = filter_significant_pairs(df_corr, threshold=CORRELATION_THRESHOLD)
    print(f"\n Perechi cu |r| >= {CORRELATION_THRESHOLD}: {len(df_filtered)}")

    df_filtered.to_csv(OUT_CSV, index=False)
    print(f"[OK] Salvat: {OUT_CSV}")

    plot_top_correlations(df_filtered, df_snp, df_exp, OUT_HEATMAP)
    
    print("\n" + "="*60)
    print("SUMAR CROSS-OMICS ANALYSIS")
    print("="*60)
    print(f"SNPs analizate:           {len(df_snp)}")
    print(f"Gene analizate:           {len(df_exp)}")
    print(f"Total perechi testate:    {len(df_corr)}")
    print(f"Perechi semnificative:    {len(df_filtered)} (|r| >= {CORRELATION_THRESHOLD})")
    print(f"SNPs cu corelatii:        {df_filtered['snp'].nunique()}")
    print(f"Gene cu corelatii:        {df_filtered['gene'].nunique()}")
    print("="*60)
