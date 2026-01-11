"""
Exercise 9.1 — Drug–Gene Bipartite Network & Drug Similarity Network

Scop:
- sa construiti o retea bipartita drug–gene plecand de la un CSV
- sa proiectati layer-ul de medicamente folosind similaritatea dintre seturile de gene
- sa exportati un fisier cu muchiile de similaritate intre medicamente

TODO:
- incarcati datele drug–gene
- construiti dict-ul drug -> set de gene tinta
- construiti graful bipartit drug–gene (NetworkX)
- calculati similaritatea dintre medicamente (ex. Jaccard)
- construiti graful drug similarity
- exportati tabelul cu muchii: drug1, drug2, weight
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, Tuple, List

import itertools

import networkx as nx
import pandas as pd

HANDLE = "rosestoica"

DRUG_GENE_CSV = Path(f"data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv")

OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_SUMMARY = OUT_DIR / f"drug_summary_{HANDLE}.csv"
OUT_DRUG_SIMILARITY = OUT_DIR / f"drug_similarity_{HANDLE}.csv"
OUT_GRAPH_DRUG_GENE = OUT_DIR / f"network_drug_gene_{HANDLE}.gpickle"


def ensure_exists(path: Path) -> None:
    """
    Verifica ca fisierul exista, altfel ridica FileNotFoundError.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Nu am gasit fisierul: {path}")


def load_drug_gene_table(path: Path) -> pd.DataFrame:
    """
    Citeste CSV-ul si valideaza ca exista coloanele necesare.
    """
    df = pd.read_csv(path)
    required_cols = {"drug", "gene"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV-ul trebuie sa contina coloanele: {required_cols}")
    return df


def build_drug2genes(df: pd.DataFrame) -> Dict[str, Set[str]]:
    """
    Construieste un dict: drug -> set de gene tinta.
    """
    drug2genes = df.groupby("drug")["gene"].apply(set).to_dict()
    return drug2genes


def build_bipartite_graph(drug2genes: Dict[str, Set[str]]) -> nx.Graph:
    """
    Construieste graful bipartit drug-gene.
    """
    G = nx.Graph()
    
    for drug in drug2genes:
        G.add_node(drug, bipartite="drug")
    
    all_genes = set()
    for genes in drug2genes.values():
        all_genes.update(genes)
    
    for gene in all_genes:
        G.add_node(gene, bipartite="gene")
    
    for drug, genes in drug2genes.items():
        for gene in genes:
            G.add_edge(drug, gene)
    
    return G


def summarize_drugs(drug2genes: Dict[str, Set[str]]) -> pd.DataFrame:
    """
    Construieste un DataFrame cu sumarul medicamentelor.
    """
    data = [
        {"drug": drug, "num_targets": len(genes)}
        for drug, genes in drug2genes.items()
    ]
    df = pd.DataFrame(data)
    df = df.sort_values("num_targets", ascending=False).reset_index(drop=True)
    return df


def jaccard_similarity(s1: Set[str], s2: Set[str]) -> float:
    """
    Calculati similaritatea Jaccard intre doua seturi de gene:
    J(A, B) = |A ∩ B| / |A ∪ B|
    """
    if not s1 and not s2:
        return 0.0
    inter = len(s1 & s2)
    union = len(s1 | s2)
    return inter / union if union > 0 else 0.0


def compute_drug_similarity_edges(
    drug2genes: Dict[str, Set[str]],
    min_sim: float = 0.0,
) -> List[Tuple[str, str, float]]:
    """
    Calculeaza similaritatea Jaccard intre toate perechile de medicamente.
    Returneaza muchiile cu similaritate >= min_sim.
    """
    edges = []
    drugs = list(drug2genes.keys())
    
    for drug1, drug2 in itertools.combinations(drugs, 2):
        sim = jaccard_similarity(drug2genes[drug1], drug2genes[drug2])
        if sim >= min_sim:
            edges.append((drug1, drug2, sim))
    
    edges.sort(key=lambda x: x[2], reverse=True)
    return edges


def edges_to_dataframe(edges: List[Tuple[str, str, float]]) -> pd.DataFrame:
    """
    Transforma lista de muchii intr-un DataFrame.
    """
    df = pd.DataFrame(edges, columns=["drug1", "drug2", "similarity"])
    return df

if __name__ == "__main__":
    ensure_exists(DRUG_GENE_CSV)
    print(f"[INFO] Fisier gasit: {DRUG_GENE_CSV}")

    df = load_drug_gene_table(DRUG_GENE_CSV)
    print(f"[INFO] Incarcat {len(df)} interactiuni drug-gene")

    drug2genes = build_drug2genes(df)
    print(f"[INFO] Numar medicamente: {len(drug2genes)}")

    G_bipartite = build_bipartite_graph(drug2genes)
    drug_nodes = [n for n, d in G_bipartite.nodes(data=True) if d.get("bipartite") == "drug"]
    gene_nodes = [n for n, d in G_bipartite.nodes(data=True) if d.get("bipartite") == "gene"]
    print(f"[INFO] Graf bipartit: {len(drug_nodes)} drugs, {len(gene_nodes)} genes, {G_bipartite.number_of_edges()} edges")

    import pickle
    with open(OUT_GRAPH_DRUG_GENE, "wb") as f:
        pickle.dump(G_bipartite, f)
    print(f"[INFO] Graf salvat in: {OUT_GRAPH_DRUG_GENE}")

    df_summary = summarize_drugs(drug2genes)
    df_summary.to_csv(OUT_DRUG_SUMMARY, index=False)
    print(f"[INFO] Sumar medicamente salvat in: {OUT_DRUG_SUMMARY}")
    print(df_summary.head(10).to_string(index=False))

    MIN_SIMILARITY = 0.1  # pastram doar muchiile cu similaritate >= 10%
    edges = compute_drug_similarity_edges(drug2genes, min_sim=MIN_SIMILARITY)
    print(f"\n[INFO] Muchii de similaritate (Jaccard >= {MIN_SIMILARITY}): {len(edges)}")
    
    df_similarity = edges_to_dataframe(edges)
    df_similarity.to_csv(OUT_DRUG_SIMILARITY, index=False)
    print(f"[INFO] Similaritati salvate in: {OUT_DRUG_SIMILARITY}")
    
    print("\n[INFO] Top 10 perechi de medicamente similare:")
    print(df_similarity.head(10).to_string(index=False))
