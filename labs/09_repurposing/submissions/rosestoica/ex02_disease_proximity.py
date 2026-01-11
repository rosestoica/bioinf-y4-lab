"""
Exercise 9.2 — Disease Proximity and Drug Ranking

Scop:
- sa calculati distanta medie dintre fiecare medicament si un set de gene asociate unei boli
- sa ordonati medicamentele in functie de proximitate (network-based prioritization)

TODO-uri principale:
- incarcati graful bipartit drug–gene (din exercitiul 9.1)
- incarcati lista de disease genes
- pentru fiecare medicament, calculati distanta minima / medie pana la genele bolii
- exportati un tabel cu medicamente si scorul lor de proximitate
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, List, Tuple

import networkx as nx
import pandas as pd

HANDLE = "rosestoica"


GRAPH_DRUG_GENE = Path(f"labs/09_repurposing/submissions/{HANDLE}/network_drug_gene_{HANDLE}.gpickle")
DRUG_GENE_CSV = Path(f"data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv")


DISEASE_GENES_TXT = Path(f"data/work/{HANDLE}/lab09/disease_genes_{HANDLE}.txt")

OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_PRIORITY = OUT_DIR / f"drug_priority_{HANDLE}.csv"

def ensure_exists(path: Path) -> None:
    """
    Verifica ca fisierul exista, altfel ridica FileNotFoundError.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Nu am gasit fisierul: {path}")


def load_bipartite_graph_or_build() -> nx.Graph:
    import pickle
    
    if GRAPH_DRUG_GENE.is_file():
        with open(GRAPH_DRUG_GENE, "rb") as f:
            G = pickle.load(f)
        print(f"[INFO] Graf incarcat din: {GRAPH_DRUG_GENE}")
        return G
    
    ensure_exists(DRUG_GENE_CSV)
    df = pd.read_csv(DRUG_GENE_CSV)
    
    G = nx.Graph()
    drug2genes = df.groupby("drug")["gene"].apply(set).to_dict()
    
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
    
    print(f"[INFO] Graf reconstruit din: {DRUG_GENE_CSV}")
    return G


def load_disease_genes(path: Path) -> Set[str]:
    """
    Incarca fisierul text cu gene (una pe linie).
    """
    with open(path, "r") as f:
        genes = {line.strip() for line in f if line.strip()}
    return genes


def get_drug_nodes(B: nx.Graph) -> List[str]:
    """
    Extrage lista nodurilor de tip 'drug'.
    """
    drugs = [n for n, d in B.nodes(data=True) if d.get("bipartite") == "drug"]
    return drugs


def compute_drug_disease_distance(
    B: nx.Graph,
    drug: str,
    disease_genes: Set[str],
    mode: str = "mean",
    max_dist: int = 5,
) -> float:
    """
    Calculeaza distanta dintre un medicament si genele bolii.
    - mode="mean": media distantelor
    - mode="min": distanta minima
    """
    distances = []
    
    for gene in disease_genes:
        if gene not in B.nodes:
            continue
        try:
            dist = nx.shortest_path_length(B, drug, gene)
            distances.append(dist)
        except nx.NetworkXNoPath:
            # Nu exista drum -> penalizam
            distances.append(max_dist + 1)
    
    if not distances:
        return float(max_dist + 1)
    
    if mode == "min":
        return min(distances)
    else:  # mean
        return sum(distances) / len(distances)


def rank_drugs_by_proximity(
    B: nx.Graph,
    disease_genes: Set[str],
    mode: str = "mean",
) -> pd.DataFrame:
    """
    Calculeaza ranking-ul medicamentelor dupa proximitatea fata de genele bolii.
    """
    drugs = get_drug_nodes(B)
    
    results = []
    for drug in drugs:
        dist = compute_drug_disease_distance(B, drug, disease_genes, mode=mode)
        results.append({"drug": drug, "distance": dist})
    
    df = pd.DataFrame(results)
    df = df.sort_values("distance", ascending=True).reset_index(drop=True)

    df["proximity_score"] = 1 / (df["distance"] + 0.1)
    
    return df

if __name__ == "__main__":
    ensure_exists(DISEASE_GENES_TXT)
    print(f"[INFO] Fisier disease genes gasit: {DISEASE_GENES_TXT}")

    G = load_bipartite_graph_or_build()
    drug_nodes = get_drug_nodes(G)
    print(f"[INFO] Graf cu {len(drug_nodes)} medicamente si {G.number_of_nodes() - len(drug_nodes)} gene")

    disease_genes = load_disease_genes(DISEASE_GENES_TXT)
    print(f"[INFO] Disease genes: {len(disease_genes)} gene")
    print(f"    Gene: {disease_genes}")
    
    genes_in_graph = disease_genes & set(G.nodes)
    print(f"[INFO] Disease genes prezente in graf: {len(genes_in_graph)}")

    print("\n[INFO] Calculam proximitatea medicamentelor fata de genele bolii...")
    df_priority = rank_drugs_by_proximity(G, disease_genes, mode="mean")

    df_priority.to_csv(OUT_DRUG_PRIORITY, index=False)
    print(f"[INFO] Ranking salvat in: {OUT_DRUG_PRIORITY}")
