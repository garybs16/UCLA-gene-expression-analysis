"""Evaluate whether transcriptomic embeddings retain biological label signal."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
)
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def pc_columns(columns: pd.Index) -> list[str]:
    return sorted(
        [col for col in columns if re.fullmatch(r"pc\d+", col)],
        key=lambda col: int(col[2:]),
    )


def save_confusion_matrix(y_true: pd.Series, y_pred: np.ndarray, labels: list[str], output_path: Path) -> None:
    matrix = confusion_matrix(y_true, y_pred, labels=labels, normalize="true")

    fig, ax = plt.subplots(figsize=(9, 7))
    image = ax.imshow(matrix, cmap="Blues", vmin=0, vmax=1)
    ax.set_title("Mechanism-of-Action Prediction Confusion Matrix")
    ax.set_xlabel("Predicted label")
    ax.set_ylabel("True label")
    ax.set_xticks(np.arange(len(labels)), labels=labels, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(labels)), labels=labels)
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("Class-normalized share")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def evaluate(embeddings: pd.DataFrame, label_col: str, output_dir: Path) -> dict[str, float]:
    pcs = pc_columns(embeddings.columns)
    if not pcs:
        raise ValueError("no PC columns found in embeddings")
    if label_col not in embeddings.columns:
        raise ValueError(f"label column not found: {label_col}")

    data = embeddings.dropna(subset=[label_col]).copy()
    x = data[pcs]
    y = data[label_col].astype(str)

    x_train, x_test, y_train, y_test = train_test_split(
        x,
        y,
        test_size=0.25,
        random_state=42,
        stratify=y,
    )

    model = Pipeline(
        steps=[
            ("scale", StandardScaler()),
            (
                "classifier",
                LogisticRegression(
                    class_weight="balanced",
                    max_iter=2_000,
                    random_state=42,
                ),
            ),
        ]
    )
    baseline = DummyClassifier(strategy="most_frequent")

    model.fit(x_train, y_train)
    baseline.fit(x_train, y_train)

    predictions = model.predict(x_test)
    baseline_predictions = baseline.predict(x_test)

    metrics = {
        "accuracy": float(accuracy_score(y_test, predictions)),
        "balanced_accuracy": float(balanced_accuracy_score(y_test, predictions)),
        "macro_f1": float(f1_score(y_test, predictions, average="macro")),
        "weighted_f1": float(f1_score(y_test, predictions, average="weighted")),
        "baseline_accuracy": float(accuracy_score(y_test, baseline_predictions)),
        "baseline_balanced_accuracy": float(balanced_accuracy_score(y_test, baseline_predictions)),
        "classes": float(y.nunique()),
        "train_rows": float(len(x_train)),
        "test_rows": float(len(x_test)),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [{"metric": metric, "value": value} for metric, value in metrics.items()]
    ).to_csv(output_dir / "moa_prediction_metrics.csv", index=False)

    report = classification_report(y_test, predictions, output_dict=True, zero_division=0)
    pd.DataFrame(report).transpose().to_csv(output_dir / "moa_classification_report.csv")

    labels = sorted(y.unique())
    save_confusion_matrix(
        y_true=y_test,
        y_pred=predictions,
        labels=labels,
        output_path=output_dir / "moa_confusion_matrix.png",
    )

    return metrics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate mechanism-of-action signal in embeddings.")
    parser.add_argument("--embeddings", default="embeddings.csv")
    parser.add_argument("--label-col", default="moa-fine")
    parser.add_argument("--output-dir", default="outputs/evaluation")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    embeddings = pd.read_csv(args.embeddings)
    metrics = evaluate(
        embeddings=embeddings,
        label_col=args.label_col,
        output_dir=Path(args.output_dir),
    )
    for metric, value in metrics.items():
        print(f"{metric}: {value:.4f}")


if __name__ == "__main__":
    main()
