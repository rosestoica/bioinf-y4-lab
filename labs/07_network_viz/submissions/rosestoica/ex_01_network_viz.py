"""
Exercise 8 — Visualization of Co-Expression Networks + Hub Genes

TODO:
- Load the expression matrix and module mapping from Lab 6
- Rebuild (or load) the adjacency matrix
- Construct the graph from adjacency
- Color nodes by module
- Compute hub genes (top degree)
- Visualize and export the network figure (.png)
- Export hub genes to CSV (optional)
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


# --------------------------
# Config — complete with your values
# --------------------------
HANDLE = "rosestoica"

# Input files
BASE_DIR = Path(__file__).resolve().parents[3]
EXPR_CSV = BASE_DIR / "data" / "work" / HANDLE / "lab06" / "expression_matrix.csv"
MODULES_CSV = BASE_DIR / "labs" / "06_wgcna" / "submissions" / HANDLE / f"modules_{HANDLE}.csv"

# Optional: if you saved adjacency in Lab 6, load it here
PRECOMPUTED_ADJ_CSV: Optional[Path] = None

# Parameters for adjacency reconstruction (if PRECOMPUTED_ADJ_CSV is None)
CORR_METHOD = "spearman"   # Spearman - mai robust pentru date non-normale
USE_ABS_CORR = True        # Folosim valori absolute ale corelatiei
ADJ_THRESHOLD = 0.7        # Prag corelatie (consistent cu Lab 6)
VARIANCE_THRESHOLD = 1.0   # Prag varianta pentru filtrare gene
TOP_GENES = 500            # Top N gene dupa varianta
WEIGHTED = False           # Adiacenta binara

# Visualization parameters
SEED = 42
TOPK_HUBS = 10
NODE_BASE_SIZE = 60
EDGE_ALPHA = 0.15

# Outputs
OUT_DIR = Path(__file__).resolve().parent
OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"


# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    """Check that a file exists, raise error if not."""
    if not path.exists():
        raise FileNotFoundError(f"Fisierul nu exista: {path}")
    print(f"   ✓ Gasit: {path.name}")


def read_expression_matrix(path: Path) -> pd.DataFrame:
    """
    Citeste matricea de expresie din CSV.
    - index = gene names
    - columns = sample IDs
    """
    df = pd.read_csv(path, index_col=0)
    return df


def read_modules_csv(path: Path) -> Dict[str, int]:
    """
    Citeste CSV-ul cu mapping gene -> modul.
    Columns: Gene, Module
    Returns: dict gene -> module_id
    """
    df = pd.read_csv(path)
    return dict(zip(df["Gene"], df["Module"]))


def correlation_to_adjacency(expr: pd.DataFrame,
                             method: str,
                             use_abs: bool,
                             threshold: float,
                             weighted: bool) -> pd.DataFrame:
    """
    Calculeaza matricea de adiacenta din corelatii.
    - compute correlation matrix pe expr (gene x gene)
    - aplica abs() daca use_abs=True
    - aplica threshold pentru a construi adiacenta
    - elimina diagonala
    """
    # Transpunem pentru ca .corr() calculeaza corelatia intre coloane
    corr = expr.T.corr(method=method)
    
    if use_abs:
        corr = corr.abs()
    
    # Setam diagonala la 0
    np.fill_diagonal(corr.values, 0)
    
    # Construim adiacenta
    if weighted:
        adj = corr.where(corr >= threshold, 0)
    else:
        adj = (corr >= threshold).astype(int)
    
    return adj


def graph_from_adjacency(A: pd.DataFrame) -> nx.Graph:
    """
    Converteste matricea de adiacenta in graf NetworkX.
    Elimina nodurile izolate.
    """
    G = nx.from_pandas_adjacency(A)
    isolates = list(nx.isolates(G))
    if isolates:
        G.remove_nodes_from(isolates)
    return G


def color_map_from_modules(nodes: Iterable[str], gene2module: Dict[str, int]) -> list:
    """
    Asigneaza o culoare fiecarui nod bazat pe modulul sau.
    Folosim paleta 'tab10' din matplotlib.
    """
    cmap = plt.cm.get_cmap('tab10')
    colors = []
    for node in nodes:
        module_id = gene2module.get(node, -1)
        if module_id >= 0:
            colors.append(cmap(module_id % 10))
        else:
            colors.append('lightgray')
    return colors


def compute_hubs(G: nx.Graph, topk: int) -> pd.DataFrame:
    """
    Calculeaza hub genes bazat pe degree si betweenness centrality.
    Returneaza top-k gene.
    """
    # Calculeaza degree
    degrees = dict(G.degree())
    
    # Calculeaza betweenness centrality
    betweenness = nx.betweenness_centrality(G)
    
    # Creeaza DataFrame
    df = pd.DataFrame({
        'Gene': list(degrees.keys()),
        'Degree': list(degrees.values()),
        'Betweenness': [betweenness[node] for node in degrees.keys()]
    })
    
    # Sorteaza dupa degree si returneaza top-k
    df = df.sort_values('Degree', ascending=False).head(topk)
    df = df.reset_index(drop=True)
    
    return df


def log_and_filter(df: pd.DataFrame,
                   variance_threshold: float,
                   top_n: int = None) -> pd.DataFrame:
    """
    Preprocesare:
    - aplica log2(x+1)
    - filtreaza genele cu varianta scazuta
    - optional: selecteaza top N gene dupa varianta
    """
    df_log = np.log2(df + 1)
    gene_variance = df_log.var(axis=1)
    df_filtered = df_log[gene_variance >= variance_threshold]
    
    if top_n is not None and len(df_filtered) > top_n:
        gene_variance_filtered = df_filtered.var(axis=1)
        top_genes = gene_variance_filtered.nlargest(top_n).index
        df_filtered = df_filtered.loc[top_genes]
    
    return df_filtered


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    print("=== Network Visualization & Hub Genes ===")
    print(f"Handle: {HANDLE}")
    print()
    
    # 1. Verify input files exist
    print("1. Verificare fisiere input...")
    ensure_exists(EXPR_CSV)
    ensure_exists(MODULES_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 2. Load expression matrix and module mapping
    print("\n2. Incarcare date...")
    expr = read_expression_matrix(EXPR_CSV)
    print(f"   Matrice expresie: {expr.shape[0]} gene × {expr.shape[1]} probe")
    
    gene2module = read_modules_csv(MODULES_CSV)
    print(f"   Module incarcate: {len(gene2module)} gene in {len(set(gene2module.values()))} module")

    # 3. Preprocesare si reconstruire adiacenta
    print("\n3. Preprocesare si construire adiacenta...")
    if PRECOMPUTED_ADJ_CSV and PRECOMPUTED_ADJ_CSV.exists():
        A = pd.read_csv(PRECOMPUTED_ADJ_CSV, index_col=0)
        # Filtreaza doar genele din module
        module_genes = [g for g in A.index if g in gene2module]
        A = A.loc[module_genes, module_genes]
        print(f"   Adiacenta incarcata din fisier: {A.shape}")
    else:
        # Preprocesare similar cu Lab 6
        expr_filtered = log_and_filter(expr, VARIANCE_THRESHOLD, TOP_GENES)
        print(f"   Gene dupa filtrare: {len(expr_filtered)}")
        
        # Calculeaza adiacenta
        A = correlation_to_adjacency(expr_filtered, CORR_METHOD, USE_ABS_CORR, ADJ_THRESHOLD, WEIGHTED)
        print(f"   Adiacenta calculata: {A.shape}")

    # 4. Build graph
    print("\n4. Construire graf...")
    G = graph_from_adjacency(A)
    print(f"   Graf: {G.number_of_nodes()} noduri, {G.number_of_edges()} muchii")
    
    if G.number_of_nodes() == 0:
        print("   Eroare: Graful nu are noduri. Ajustati pragurile.")
        exit(1)

    # 5. Compute colors by module
    print("\n5. Asignare culori dupa modul...")
    node_colors = color_map_from_modules(G.nodes(), gene2module)
    n_colored = sum(1 for n in G.nodes() if n in gene2module)
    print(f"   Noduri colorate: {n_colored}/{G.number_of_nodes()}")

    # 6. Compute hub genes
    print(f"\n6. Calculare hub genes (top {TOPK_HUBS})...")
    hubs_df = compute_hubs(G, TOPK_HUBS)
    print("   Hub genes:")
    for _, row in hubs_df.iterrows():
        print(f"      {row['Gene']}: degree={row['Degree']}, betweenness={row['Betweenness']:.4f}")
    
    # Dimensiuni noduri bazate pe degree
    degrees = dict(G.degree())
    max_degree = max(degrees.values()) if degrees else 1
    node_sizes = [NODE_BASE_SIZE + (degrees[n] / max_degree) * 200 for n in G.nodes()]
    
    hub_genes = set(hubs_df['Gene'].tolist())

    # 7. Compute layout and draw graph
    print("\n7. Vizualizare retea...")
    plt.figure(figsize=(14, 10))
    
    # Layout
    pos = nx.spring_layout(G, seed=SEED, k=2/np.sqrt(G.number_of_nodes()))
    
    # Draw edges
    nx.draw_networkx_edges(G, pos, alpha=EDGE_ALPHA, edge_color='gray', width=0.5)
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8)
    
    # Draw labels only for hubs
    hub_labels = {n: n for n in G.nodes() if n in hub_genes}
    nx.draw_networkx_labels(G, pos, labels=hub_labels, font_size=8, font_weight='bold')
    
    plt.title(f"Gene Co-Expression Network ({HANDLE})\n{G.number_of_nodes()} nodes, {G.number_of_edges()} edges, {len(set(gene2module.values()))} modules")
    plt.axis('off')
    plt.tight_layout()

    # 8. Save network figure
    print(f"\n8. Salvare figura...")
    plt.savefig(OUT_PNG, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"   Salvat: {OUT_PNG}")

    # 9. Save hub genes to CSV
    print(f"\n9. Salvare hub genes...")
    hubs_df.to_csv(OUT_HUBS, index=False)
    print(f"   Salvat: {OUT_HUBS}")

    print("\n=== Vizualizare completa ===")
