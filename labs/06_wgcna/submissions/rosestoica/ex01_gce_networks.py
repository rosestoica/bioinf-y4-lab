"""
Exercițiu Gene Co-Expression Networks (GCEs) — Construirea rețelei și detectarea modulelor

Obiectiv:
- Să construiți o rețea de co-expresie dintr-o matrice de expresie RNA-Seq
- Să detectați module (comunități) de gene folosind un algoritm de tip Louvain (sau alternativ)

Instrucțiuni (în laborator):
1) Pregătire date
   - Descărcați și pregătiți matricea de expresie (ex: GSE115469) într-un CSV cu:
     * rânduri = gene (index), coloane = probe (sample IDs)
   - Salvați fișierul la: data/work/<handle>/lab06/expression_matrix.csv

2) Preprocesare
   - log2(x + 1)
   - filtrare gene cu varianță scăzută

3) Corelație → Adiacență
   - completați funcția `correlation_matrix`
   - funcția `adjacency_from_correlation` este deja implementată

4) Graf + Module
   - construiți graful cu NetworkX
   - detectați modulele (Louvain sau alternativă)
   - exportați mapping-ul gene → modul în submissions/<handle>/modules_<handle>.csv

Notă:
- Documentați în <github_handle>_notes.md: metrica de corelație, pragul, observații scurte.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable
import gzip
import urllib.request
import io

import numpy as np
import pandas as pd
import networkx as nx


# --------------------------
# Config — completați după nevoie
# --------------------------
HANDLE = "rosestoica"
DATA_DIR = Path(__file__).resolve().parents[4] / "data" / "work" / HANDLE / "lab06"
INPUT_CSV = DATA_DIR / "expression_matrix.csv"
OUTPUT_DIR = Path(__file__).resolve().parent
OUTPUT_CSV = OUTPUT_DIR / f"modules_{HANDLE}.csv"

# GEO Dataset - GSE48350 (Human brain aging study)
# Dataset cu ~173 probe, potrivit pentru analize de co-expresie
GEO_ACCESSION = "GSE48350"
GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48350/matrix/GSE48350_series_matrix.txt.gz"

CORR_METHOD = "spearman"   # Spearman este mai robust pentru date non-normale
VARIANCE_THRESHOLD = 1.0   # prag pentru filtrare gene cu varianță scăzută
ADJ_THRESHOLD = 0.7        # prag pentru |cor| (valori >= 0.7 sunt considerate puternic corelate)
TOP_GENES = 500            # Limităm la top N gene cu varianță mare pentru performanță
USE_ABS_CORR = True        # True => folosiți |cor| la prag
MAKE_UNDIRECTED = True     # rețelele de co-expresie sunt de obicei neorientate


def download_geo_matrix(geo_url: str, output_path: Path) -> None:
    """
    Descarcă matricea de expresie de la GEO (NCBI) și o salvează local.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if output_path.exists():
        print(f"   Fișierul există deja: {output_path}")
        return
    
    print(f"   Descărcare de la: {geo_url}")
    
    # Descarcă fișierul gzip
    with urllib.request.urlopen(geo_url) as response:
        compressed_data = response.read()
    
    # Decomprimă și parsează
    with gzip.open(io.BytesIO(compressed_data), 'rt') as f:
        lines = f.readlines()
    
    # Găsește începutul datelor (după !series_matrix_table_begin)
    start_idx = None
    end_idx = None
    for i, line in enumerate(lines):
        if line.startswith('!series_matrix_table_begin'):
            start_idx = i + 1
        elif line.startswith('!series_matrix_table_end'):
            end_idx = i
            break
    
    if start_idx is None or end_idx is None:
        raise ValueError("Nu s-a găsit tabelul de date în fișierul GEO")
    
    # Extrage datele
    data_lines = lines[start_idx:end_idx]
    
    # Salvează ca CSV
    with open(output_path, 'w') as f:
        for line in data_lines:
            # Convertește tab-separated la comma-separated
            f.write(line.replace('\t', ',').replace('"', ''))
    
    print(f"   Salvat în: {output_path}")


def read_expression_matrix(path: Path) -> pd.DataFrame:
    """Citește matricea de expresie din CSV. Genele sunt pe rânduri, probele pe coloane."""
    df = pd.read_csv(path, index_col=0)
    return df


def log_and_filter(df: pd.DataFrame,
                   variance_threshold: float,
                   top_n: int = None) -> pd.DataFrame:
    """
    Preprocesare:
    - aplică log2(x+1)
    - filtrează genele cu varianță scăzută
    - opțional: selectează top N gene după varianță
    """
    # Aplică log2(x + 1) pentru normalizare
    df_log = np.log2(df + 1)
    
    # Calculează varianța pentru fiecare genă (pe rânduri)
    gene_variance = df_log.var(axis=1)
    
    # Filtrează genele cu varianță sub prag
    df_filtered = df_log[gene_variance >= variance_threshold]
    
    print(f"   Gene rămase după filtrare varianță: {len(df_filtered)} din {len(df_log)}")
    
    # Selectează top N gene după varianță pentru performanță
    if top_n is not None and len(df_filtered) > top_n:
        gene_variance_filtered = df_filtered.var(axis=1)
        top_genes = gene_variance_filtered.nlargest(top_n).index
        df_filtered = df_filtered.loc[top_genes]
        print(f"   Selectate top {top_n} gene după varianță")
    
    return df_filtered

def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    """
    Calculează matricea de corelație între gene (rânduri).
    
    Args:
        df: DataFrame cu gene pe rânduri și probe pe coloane
        method: "pearson" sau "spearman"
        use_abs: dacă True, returnează valori absolute ale corelației
    
    Returns:
        Matrice de corelație (gene x gene)
    """
    # Transpunem pentru că .corr() calculează corelația între coloane
    corr = df.T.corr(method=method)
    
    if use_abs:
        corr = corr.abs()
    
    # Setăm diagonala la 0 pentru a evita self-loops
    np.fill_diagonal(corr.values, 0)
    
    return corr


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float,
                               weighted: bool = False) -> pd.DataFrame:
    """
    Construiți matricea de adiacență din corelații.
    - binară: A_ij = 1 dacă corr_ij >= threshold, altfel 0
    - ponderată: A_ij = corr_ij dacă corr_ij >= threshold, altfel 0
    """
    if weighted:
        # Păstrăm valorile corelației unde depășesc pragul
        adj = corr.where(corr >= threshold, 0)
    else:
        # Matrice binară: 1 dacă corelația >= prag, altfel 0
        adj = (corr >= threshold).astype(int)
    
    return adj


def graph_from_adjacency(A: pd.DataFrame,
                         undirected: bool = True) -> nx.Graph:
    if undirected:
        G = nx.from_pandas_adjacency(A)
    else:
        G = nx.from_pandas_adjacency(A, create_using=nx.DiGraph)
    isolates = list(nx.isolates(G))
    if isolates:
        G.remove_nodes_from(isolates)
    return G


def detect_modules_louvain_or_greedy(G: nx.Graph) -> Dict[str, int]:
    """
    Detectează comunități (module) și întoarce un dict gene -> modul_id.
    Variante:
      - încercați louvain_communities(G, seed=42) dacă e disponibil
      - altfel greedy_modularity_communities(G)
    """
    gene_to_module: Dict[str, int] = {}
    
    try:
        # Încercăm algoritmul Louvain (disponibil în NetworkX >= 2.7)
        communities = nx.community.louvain_communities(G, seed=42)
        print("Folosim algoritmul Louvain pentru detectarea modulelor.")
    except AttributeError:
        # Fallback la greedy modularity dacă Louvain nu este disponibil
        communities = list(nx.community.greedy_modularity_communities(G))
        print("Folosim algoritmul Greedy Modularity pentru detectarea modulelor.")
    
    # Construim mapping-ul gene -> modul_id
    for module_id, community in enumerate(communities):
        for gene in community:
            gene_to_module[gene] = module_id
    
    return gene_to_module


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_modules = (
        pd.DataFrame({"Gene": list(mapping.keys()), "Module": list(mapping.values())})
        .sort_values(["Module", "Gene"])
    )
    df_modules.to_csv(out_csv, index=False)


if __name__ == "__main__":
    print(f"=== Gene Co-Expression Network Analysis ===")
    print(f"Dataset: {GEO_ACCESSION} de la NCBI GEO")
    print(f"Output: {OUTPUT_CSV}")
    print()
    
    # 0. Descarcă datele de la GEO dacă nu există
    print("0. Descărcare date de la NCBI GEO...")
    download_geo_matrix(GEO_URL, INPUT_CSV)
    
    # 1. Citește matricea de expresie
    print("\n1. Citire matrice de expresie...")
    expr_df = read_expression_matrix(INPUT_CSV)
    print(f"   Dimensiuni inițiale: {expr_df.shape[0]} gene × {expr_df.shape[1]} probe")
    
    # 2. Preprocesare: log2 și filtrare
    print("\n2. Preprocesare (log2 + filtrare varianță + top genes)...")
    expr_filtered = log_and_filter(expr_df, VARIANCE_THRESHOLD, top_n=TOP_GENES)
    
    # 3. Calculează matricea de corelație
    print(f"\n3. Calculare matrice de corelație ({CORR_METHOD})...")
    corr_matrix = correlation_matrix(expr_filtered, method=CORR_METHOD, use_abs=USE_ABS_CORR)
    print(f"   Dimensiuni matrice corelație: {corr_matrix.shape}")
    
    # 4. Construiește matricea de adiacență
    print(f"\n4. Construire matrice de adiacență (prag={ADJ_THRESHOLD})...")
    adj_matrix = adjacency_from_correlation(corr_matrix, threshold=ADJ_THRESHOLD, weighted=False)
    n_edges = (adj_matrix.values.sum() - np.trace(adj_matrix.values)) // 2
    print(f"   Număr muchii potențiale: {n_edges}")
    
    # 5. Construiește graful
    print("\n5. Construire graf NetworkX...")
    G = graph_from_adjacency(adj_matrix, undirected=MAKE_UNDIRECTED)
    print(f"   Graf creat cu {G.number_of_nodes()} noduri și {G.number_of_edges()} muchii.")
    
    # 6. Detectează module
    print("\n6. Detectare module (comunități)...")
    if G.number_of_nodes() > 0:
        gene_to_module = detect_modules_louvain_or_greedy(G)
        print(f"   S-au detectat {len(set(gene_to_module.values()))} module.")
        
        # 7. Salvează rezultatele
        print(f"\n7. Salvare rezultate...")
        save_modules_csv(gene_to_module, OUTPUT_CSV)
        print(f"   Am salvat mapping-ul gene→modul în: {OUTPUT_CSV}")
    else:
        print("   Avertisment: Graful nu are noduri. Ajustați pragurile.")
    
    print("\n=== Analiză completă ===")
