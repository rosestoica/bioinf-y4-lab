"""
Exercise 8 — Supervised ML pipeline pentru expresie genica (Random Forest)

TODO-uri principale:
- Incarcati matricea de expresie (ex. subset TP53 / GTEx) pentru HANDLE-ul vostru
- Separati features (gene) si label (ultima coloana)
- Encodati etichetele
- Impartiti in train/test
- Antrenati un RandomForestClassifier (model de baza)
- Evaluati: classification_report + matrice de confuzie (salvate)
- Calculati importanta trasaturilor si salvati in CSV
- (Optional) Aplicati KMeans pe X si comparati clustere vs etichete reale
"""

from __future__ import annotations
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

HANDLE = "rosestoica"

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]  
DATA_CSV = PROJECT_ROOT / f"data/work/{HANDLE}/lab08/expression_matrix_{HANDLE}.csv"

TEST_SIZE = 0.2
RANDOM_STATE = 42
N_ESTIMATORS = 200
TOPK_FEATURES = 20

OUT_DIR = SCRIPT_DIR
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_CONFUSION = OUT_DIR / f"confusion_rf_{HANDLE}.png"
OUT_REPORT = OUT_DIR / f"classification_report_{HANDLE}.txt"
OUT_FEATIMP = OUT_DIR / f"feature_importance_{HANDLE}.csv"
OUT_CLUSTER_CROSSTAB = OUT_DIR / f"cluster_crosstab_{HANDLE}.csv"

def ensure_exists(path: Path) -> None:
    """
    Verifica daca fisierul de input exista.
    Daca nu, ridica FileNotFoundError cu un mesaj clar.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Nu am gasit fisierul: {path}")


def load_dataset(path: Path) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Citeste CSV-ul cu pandas.
    Presupunem ca ultima coloana este label-ul (ex. 'Label').
    X = toate coloanele numerice (mai putin sample_id si Label)
    y = ultima coloana (Label)
    """
    df = pd.read_csv(path)

    X = df.iloc[:, 1:-1]  # Toate coloanele gene
    y = df.iloc[:, -1]     # Coloana Label
    return X, y


def encode_labels(y: pd.Series) -> Tuple[np.ndarray, LabelEncoder]:
    """
    Foloseste LabelEncoder pentru a converti etichetele string in valori numerice.
    Returneaza y_encoded si encoder-ul.
    """
    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    return y_enc, le


def train_random_forest(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
    n_estimators: int,
    random_state: int,
) -> RandomForestClassifier:
    """
    Initializeaza si antreneaza un RandomForestClassifier.
    Returneaza modelul antrenat.
    """
    rf = RandomForestClassifier(
        n_estimators=n_estimators,
        random_state=random_state,
        n_jobs=-1
    )
    rf.fit(X_train, y_train)
    return rf


def evaluate_model(
    model: RandomForestClassifier,
    X_test: pd.DataFrame,
    y_test: np.ndarray,
    label_encoder: LabelEncoder,
    out_png: Path,
    out_txt: Path,
) -> None:
    """
    Calculeaza predictiile pe X_test, genereaza classification_report
    si salveaza matricea de confuzie ca imagine .png.
    """
    y_pred = model.predict(X_test)

    target_names = label_encoder.classes_
    report = classification_report(y_test, y_pred, target_names=target_names)
    print("\n=== Classification Report ===")
    print(report)
    out_txt.write_text(report)
    print(f"[INFO] Raport salvat in: {out_txt}")

    cm = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(6, 5))
    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="Blues",
        xticklabels=target_names,
        yticklabels=target_names,
    )
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title("Random Forest — Confusion Matrix")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"[INFO] Matrice de confuzie salvata in: {out_png}")


def compute_feature_importance(
    model: RandomForestClassifier,
    feature_names: pd.Index,
    out_csv: Path,
) -> pd.DataFrame:
    """
    Extrage feature_importances_ din model, construieste un DataFrame
    cu coloanele 'Feature' si 'Importance', sorteaza descrescator
    si salveaza in CSV.
    """
    importances = model.feature_importances_
    df_imp = pd.DataFrame(
        {"Feature": feature_names, "Importance": importances}
    ).sort_values("Importance", ascending=False)
    df_imp.to_csv(out_csv, index=False)
    print(f"[INFO] Feature importance salvata in: {out_csv}")
    print(f"\nTop 10 gene cele mai importante:")
    print(df_imp.head(10).to_string(index=False))
    return df_imp


def run_kmeans_and_crosstab(
    X: pd.DataFrame,
    y: np.ndarray,
    label_encoder: LabelEncoder,
    n_clusters: int,
    out_csv: Path,
) -> None:
    """
    Ruleaza KMeans cu n_clusters egal cu numarul de clase,
    construieste un crosstab intre eticheta reala si cluster
    si salveaza crosstab-ul in CSV.
    """
    kmeans = KMeans(n_clusters=n_clusters, random_state=RANDOM_STATE, n_init="auto")
    clusters = kmeans.fit_predict(X.values)

    df = pd.DataFrame(
        {"Label": label_encoder.inverse_transform(y), "Cluster": clusters}
    )
    ctab = pd.crosstab(df["Label"], df["Cluster"])
    ctab.to_csv(out_csv)
    print(f"\n[INFO] Crosstab salvat in: {out_csv}")
    print("\nCrosstab Label vs Cluster:")
    print(ctab)


if __name__ == "__main__":
    print("="*60)
    print("Exercise 8 — Supervised ML Pipeline (Random Forest)")
    print("="*60)
    print(f"Handle: {HANDLE}")
    print(f"Data: {DATA_CSV}")
    print(f"Output dir: {OUT_DIR}")
    print()

    ensure_exists(DATA_CSV)
    print(f"[OK] Fisierul de date exista: {DATA_CSV}")

    X, y = load_dataset(DATA_CSV)
    print(f"[OK] Date incarcate: {X.shape[0]} probe, {X.shape[1]} gene")
    print(f"     Clase: {y.unique().tolist()}")

    y_enc, le = encode_labels(y)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_enc,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        stratify=y_enc,
    )
    print(f"[OK] Train/Test split: {len(X_train)} train, {len(X_test)} test")

    print(f"\n[INFO] Antrenez Random Forest cu {N_ESTIMATORS} estimatori...")
    rf = train_random_forest(X_train, y_train, N_ESTIMATORS, RANDOM_STATE)
    print("[OK] Model antrenat!")
    
    evaluate_model(rf, X_test, y_test, le, OUT_CONFUSION, OUT_REPORT)

    feat_imp_df = compute_feature_importance(rf, X.columns, OUT_FEATIMP)

    n_classes = len(le.classes_)
    run_kmeans_and_crosstab(X, y_enc, le, n_clusters=n_classes, out_csv=OUT_CLUSTER_CROSSTAB)

    print("\n" + "="*60)
    print("[DONE] Pipeline complet! Verificati outputurile in:")
    print(f"  - {OUT_CONFUSION}")
    print(f"  - {OUT_REPORT}")
    print(f"  - {OUT_FEATIMP}")
    print(f"  - {OUT_CLUSTER_CROSSTAB}")
    print("="*60)
