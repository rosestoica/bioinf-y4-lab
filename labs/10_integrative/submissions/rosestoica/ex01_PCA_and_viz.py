"""
Exercise 10 — PCA Single-Omics vs Joint

Scop:
- incarcati SNP si Expression
- normalizati fiecare strat (z-score)
- rulati PCA pe:
    1) strat SNP
    2) strat Expression
    3) strat Joint (concat)
- generati 3 figuri PNG
- comparati vizual distributia probelor
"""

from pathlib import Path
from typing import Tuple, Dict

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

HANDLE = "rosestoica"

SNP_CSV = Path(f"data/work/{HANDLE}/lab10/snp_matrix_{HANDLE}.csv")
EXP_CSV = Path(f"data/work/{HANDLE}/lab10/expression_matrix_{HANDLE}.csv")
META_CSV = Path(f"data/work/{HANDLE}/lab10/metadata_{HANDLE}.csv")

OUT_DIR = Path(f"labs/10_integrative/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

def load_and_align_data() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Incarca datele SNP si Expression si le aliniaza pe samples comune.
    """
    df_snp = pd.read_csv(SNP_CSV, index_col=0)
    df_exp = pd.read_csv(EXP_CSV, index_col=0)
    df_meta = pd.read_csv(META_CSV)
    
    common_samples = df_snp.columns.intersection(df_exp.columns)
    
    df_snp = df_snp[common_samples]
    df_exp = df_exp[common_samples]
    
    df_meta = df_meta[df_meta["sample_id"].isin(common_samples)]
    df_meta = df_meta.set_index("sample_id").loc[common_samples].reset_index()
    
    print(f" SNP matrix: {df_snp.shape} (features x samples)")
    print(f" Expression matrix: {df_exp.shape} (features x samples)")
    print(f" Common samples: {len(common_samples)}")
    
    return df_snp, df_exp, df_meta


def normalize_zscore(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalizeaza datele cu z-score (pe features/randuri).
    """
    df_norm = (df.T - df.T.mean()) / df.T.std()
    return df_norm.T.fillna(0)


def run_pca(df: pd.DataFrame, n_components: int = 2) -> Tuple[np.ndarray, PCA]:
    """
    Ruleaza PCA pe date (samples pe randuri, features pe coloane).
    """
    X = df.T.values
    
    pca = PCA(n_components=n_components)
    proj = pca.fit_transform(X)
    
    return proj, pca


def plot_pca(
    proj: np.ndarray,
    pca: PCA,
    labels: np.ndarray,
    title: str,
    out_path: Path,
    colors: Dict[str, str] = None
) -> None:
    """
    Genereaza un scatter plot PCA colorat dupa labels.
    """
    if colors is None:
        colors = {"Group_A": "#E53935", "Group_B": "#1E88E5", "Group_C": "#43A047"}
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    for label in np.unique(labels):
        mask = labels == label
        ax.scatter(
            proj[mask, 0], proj[mask, 1],
            c=colors.get(label, "gray"),
            label=label,
            s=80,
            alpha=0.7,
            edgecolors="white",
            linewidths=0.5
        )
    
    var_explained = pca.explained_variance_ratio_ * 100
    ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}% variance)", fontsize=12)
    ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}% variance)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend(loc="best", fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    
    print(f"[OK] Salvat: {out_path}")

if __name__ == "__main__":

    df_snp, df_exp, df_meta = load_and_align_data()
    labels = df_meta["group"].values
    

    print("\n Normalizare z-score...")
    df_snp_norm = normalize_zscore(df_snp)
    df_exp_norm = normalize_zscore(df_exp)
    
    df_joint = pd.concat([df_snp_norm, df_exp_norm], axis=0)
    print(f" Joint matrix: {df_joint.shape} (features x samples)")
    joint_csv = OUT_DIR / f"multiomics_concat_{HANDLE}.csv"
    df_joint.to_csv(joint_csv)
    print(f"[OK] Salvat: {joint_csv}")
    
    print("\n Rulare PCA...")
    
    proj_snp, pca_snp = run_pca(df_snp_norm)
    plot_pca(
        proj_snp, pca_snp, labels,
        title="PCA on SNP Data (Single-Omics)",
        out_path=OUT_DIR / f"pca_snp_{HANDLE}.png"
    )
    
    proj_exp, pca_exp = run_pca(df_exp_norm)
    plot_pca(
        proj_exp, pca_exp, labels,
        title="PCA on Expression Data (Single-Omics)",
        out_path=OUT_DIR / f"pca_expression_{HANDLE}.png"
    )
    

    proj_joint, pca_joint = run_pca(df_joint)
    plot_pca(
        proj_joint, pca_joint, labels,
        title="PCA on Integrated Multi-Omics (SNP + Expression)",
        out_path=OUT_DIR / f"pca_joint_{HANDLE}.png"
    )

    print("COMPARATIE PCA - VARIANCE")
    print("="*60)
    print(f"{'Dataset':<20} {'PC1':>10} {'PC2':>10} {'Total':>10}")
    print("-"*60)
    print(f"{'SNP':<20} {pca_snp.explained_variance_ratio_[0]*100:>9.1f}% {pca_snp.explained_variance_ratio_[1]*100:>9.1f}% {sum(pca_snp.explained_variance_ratio_)*100:>9.1f}%")
    print(f"{'Expression':<20} {pca_exp.explained_variance_ratio_[0]*100:>9.1f}% {pca_exp.explained_variance_ratio_[1]*100:>9.1f}% {sum(pca_exp.explained_variance_ratio_)*100:>9.1f}%")
    print(f"{'Joint (Integrated)':<20} {pca_joint.explained_variance_ratio_[0]*100:>9.1f}% {pca_joint.explained_variance_ratio_[1]*100:>9.1f}% {sum(pca_joint.explained_variance_ratio_)*100:>9.1f}%")
