"""
Exercise 8b — Logistic Regression vs Random Forest pe expresie genica

Scop:
- sa antrenam si sa comparam doua modele:
  - Logistic Regression (multiclass, liniar)
  - Random Forest (non-liniar, bazat pe arbori)
- sa vedem daca performanta si erorile sunt similare sau diferite

TODO:
- Incarcati expresia pentru HANDLE
- Impartiti in X (gene) si y (Label)
- Encodati etichetele
- Impartiti in train/test
- Scalati features pentru logistic regression
- Antrenati RF si Logistic Regression
- Comparati classification_report pentru ambele modele
"""

from __future__ import annotations
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, StandardScaler

HANDLE = "rosestoica"

DATA_CSV = Path(f"data/work/{HANDLE}/lab08/expression_matrix_{HANDLE}.csv")

TEST_SIZE = 0.2
RANDOM_STATE = 42
N_ESTIMATORS = 200
MAX_ITER_LOGREG = 1000

OUT_DIR = Path(f"labs/08_ML_flower/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_REPORT_TXT = OUT_DIR / f"rf_vs_logreg_report_{HANDLE}.txt"

def ensure_exists(path: Path) -> None:
    """
    Verifica ca fisierul exista, altfel ridica exceptie.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Nu am gasit fisierul: {path}")


def load_dataset(path: Path) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Citeste CSV si returneaza X (features) si y (labels).
    X = toate coloanele mai putin sample_id si Label
    y = coloana Label
    """
    df = pd.read_csv(path)
    X = df.drop(columns=["sample_id", "Label"])
    y = df["Label"]
    return X, y


def encode_labels(y: pd.Series) -> Tuple[np.ndarray, LabelEncoder]:
    """
    Encodeaza etichetele folosind LabelEncoder.
    """
    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    return y_enc, le


def train_models(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
) -> Tuple[RandomForestClassifier, LogisticRegression, StandardScaler]:
    """
    Antreneaza RandomForest si LogisticRegression.
    LogisticRegression necesita scaling pentru convergenta mai buna.
    """
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)

    rf = RandomForestClassifier(
        n_estimators=N_ESTIMATORS,
        random_state=RANDOM_STATE,
        n_jobs=-1
    )
    rf.fit(X_train, y_train)

    logreg = LogisticRegression(
        multi_class="multinomial",
        max_iter=MAX_ITER_LOGREG,
        n_jobs=-1,
        random_state=RANDOM_STATE
    )
    logreg.fit(X_train_scaled, y_train)

    return rf, logreg, scaler


def compare_models(
    rf: RandomForestClassifier,
    logreg: LogisticRegression,
    scaler: StandardScaler,
    X_test: pd.DataFrame,
    y_test: np.ndarray,
    label_encoder: LabelEncoder,
    out_txt: Path,
) -> None:
    """
    Compara Random Forest si Logistic Regression.
    Genereaza classification_report pentru ambele si salveaza rezultatele.
    """
    X_test_scaled = scaler.transform(X_test)

    y_pred_rf = rf.predict(X_test)
    y_pred_logreg = logreg.predict(X_test_scaled)

    target_names = label_encoder.classes_

    report_rf = classification_report(y_test, y_pred_rf, target_names=target_names)
    report_logreg = classification_report(y_test, y_pred_logreg, target_names=target_names)

    print("=== Random Forest ===")
    print(report_rf)
    print("\n=== Logistic Regression ===")
    print(report_logreg)

    combined = (
        "=== Random Forest ===\n"
        + report_rf
        + "\n\n=== Logistic Regression ===\n"
        + report_logreg
    )
    out_txt.write_text(combined)
    print(f"\n[INFO] Raportul a fost salvat in: {out_txt}")


if __name__ == "__main__":
    ensure_exists(DATA_CSV)

    X, y = load_dataset(DATA_CSV)
    print(f"[INFO] Dataset incarcat: {X.shape[0]} samples, {X.shape[1]} features")
    print(f"[INFO] Clase: {y.value_counts().to_dict()}")

    y_enc, le = encode_labels(y)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_enc,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        stratify=y_enc,
    )
    print(f"[INFO] Train: {len(X_train)}, Test: {len(X_test)}")

    print("[INFO] Antrenare Random Forest si Logistic Regression...")
    rf, logreg, scaler = train_models(X_train, y_train)

    compare_models(rf, logreg, scaler, X_test, y_test, le, OUT_REPORT_TXT)
