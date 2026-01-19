import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import auc  # type: ignore

from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold  # type: ignore
from sklearn.model_selection import (  # Import from sklearn, not xgboost
    RandomizedSearchCV,
    train_test_split,
)
from xgboost import XGBClassifier

from python.normalization import (
    presence_absence_transform,
    rclr_transform,
    total_sum_scaling,
)


def prevalence_filtering(
    count_table_df: pd.DataFrame, prevalence_threshold: float
) -> pd.DataFrame:
    """Remove features from feature table with prevalence below threshold."""
    feature_prevalence = (
        (count_table_df > 0).sum(axis=1) / count_table_df.shape[1] * 100
    )

    # Remove features with prevalence below threshold
    return count_table_df.loc[feature_prevalence >= prevalence_threshold, :]


def plot_roc_curve(
    model: XGBClassifier, X_test: pd.DataFrame, y_test: pd.Series, filename: str
) -> None:
    """Plot ROC curve for the given model and test data."""

    y_pred_proba = model.predict_proba(X_test)[:, 1]
    fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
    roc_auc = auc(fpr, tpr)

    plt.figure(figsize=(7, 6))
    plt.plot(fpr, tpr, color="blue", lw=2, label=f"ROC curve (AUC = {roc_auc:.3f})")
    plt.plot([0, 1], [0, 1], color="gray", linestyle="--", lw=1, label="Random chance")
    plt.xlabel("False Positive Rate (1 - Specificity)")
    plt.ylabel("True Positive Rate (Sensitivity)")
    plt.title("Receiver Operating Characteristic (ROC) Curve")
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(filename, format="pdf")
    plt.close()
    print(f"ROC curve saved as '{filename}'")


def train_xgboost_model(
    count_table_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    primary_label: str,
    normalization_method: str,
    output_dir: str,
) -> None:
    """Train an XGBoost classifier on the provided data.

    Args:
        count_table_df (pd.DataFrame): Feature count table with samples as rows and features as columns.
        metadata_df (pd.DataFrame): Metadata containing the target labels.
        primary_label (str): Column name in metadata_df to use as the target variable.

    Returns:
        XGBClassifier: Trained XGBoost model.
        dict: Cross-validation scores including ROC AUC, PR AUC, accuracy, and F1 score.
    """

    # Do prevalence filtering
    prevalence_threshold = 0.2

    count_table_df = prevalence_filtering(count_table_df, prevalence_threshold)

    if normalization_method == "rclr":
        data = rclr_transform(count_table_df, qiime_orientation=True)
    elif normalization_method == "presence_absence":
        data = presence_absence_transform(count_table_df, qiime_orientation=True)
    elif normalization_method == "tss":
        data = total_sum_scaling(count_table_df, qiime_orientation=True)
    else:
        raise ValueError(
            f"Unsupported normalization method: {normalization_method}. "
            "Choose 'rclr' or 'presence_absence' or 'tss'."
        )

    # Create output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    #########################################################
    # Extract target variable from metadata
    #########################################################
    if primary_label not in metadata_df.columns:
        raise ValueError(f"Primary label '{primary_label}' not found in metadata.")
    target = metadata_df[primary_label]

    # Ensure samples align between data and metadata
    common_samples = data.index.intersection(metadata_df.index)
    if len(common_samples) == 0:
        raise ValueError("No common samples found between count table and metadata.")

    data = data.loc[common_samples]
    target = target.loc[common_samples]

    # Split data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(
        data, target, test_size=0.2, random_state=42, stratify=target
    )

    # Class imbalance (compute from TRAIN only)
    pos = int((y_train == 1).sum())
    neg = int((y_train == 0).sum())
    scale_pos_weight = max(1.0, neg / max(1, pos))

    model_params = dict(
        tree_method="hist",
        missing=np.nan,
        eval_metric="auc",
        random_state=42,
        scale_pos_weight=scale_pos_weight,
    )

    # Initialize base model
    base_model = XGBClassifier(**model_params)

    param_dist = {
        "n_estimators": [300, 500, 800, 1200],
        "learning_rate": np.linspace(0.01, 0.2, 10),
        "max_depth": [2, 3, 4, 5, 6],
        "min_child_weight": [1, 2, 3, 5, 8],
        "subsample": np.linspace(0.6, 1.0, 9),
        "colsample_bytree": np.linspace(0.4, 1.0, 13),
        "gamma": [0, 0.5, 1, 2],
        "reg_alpha": [0, 1e-3, 1e-2, 1e-1, 1, 10],
        "reg_lambda": [0.1, 0.5, 1, 2, 5, 10],
    }

    # Cross-validation setup
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    # Randomized search for hyperparameter tuning
    # Use a dictionary of scoring metrics to get all metrics in cv_results_
    scoring_metrics = {
        "roc_auc": "roc_auc",
        "average_precision": "average_precision",
        "accuracy": "accuracy",
        "f1": "f1",
        "precision": "precision",
        "recall": "recall",
    }

    search = RandomizedSearchCV(
        estimator=base_model,
        param_distributions=param_dist,
        n_iter=60,
        scoring=scoring_metrics,
        n_jobs=-1,
        cv=cv,
        refit="roc_auc",  # Refit using ROC AUC as the primary metric
        verbose=1,
        random_state=42,
    )

    # Fit the model
    print("Starting model training with hyperparameter tuning...")
    search.fit(X_train, y_train)
    print("Model training completed.")

    # Print best results
    print("Best CV AUROC: %.3f" % search.best_score_)
    print("Best params:", search.best_params_)

    # Get the best model
    best_model = search.best_estimator_

    # Cross-validation results
    cv_results = search.cv_results_
    cv_roc_auc = cv_results["mean_test_roc_auc"]
    cv_pr_auc = cv_results["mean_test_average_precision"]
    cv_accuracy = cv_results["mean_test_accuracy"]
    cv_f1 = cv_results["mean_test_f1"]
    cv_precision = cv_results["mean_test_precision"]
    cv_recall = cv_results["mean_test_recall"]

    cv_results_df = pd.DataFrame(
        {
            "ROC AUC": cv_roc_auc,
            "PR AUC": cv_pr_auc,
            "Accuracy": cv_accuracy,
            "F1-score": cv_f1,
            "Precision": cv_precision,
            "Recall": cv_recall,
        },
        index=cv_results["params"],
    )
    cv_results_df.to_csv(f"{output_dir}/cv_results.csv")

    # Path to AUC ROC plot file
    roc_plot_file = f"{output_dir}/roc_curve.pdf"

    # Plot ROC curve
    plot_roc_curve(best_model, X_test, y_test, roc_plot_file)

    # Evaluate on test set
    y_pred = best_model.predict(X_test)
    y_pred_proba = best_model.predict_proba(X_test)[:, 1]

    # Use the model's class ordering consistently across metrics
    classes = best_model.classes_

    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred, labels=classes)
    confusion_matrix_df = pd.DataFrame(
        cm,
        index=classes,
        columns=classes,
    )
    confusion_matrix_df.to_csv(f"{output_dir}/confusion_matrix.csv")

    roc_auc = roc_auc_score(y_test, y_pred_proba)
    pr_auc = average_precision_score(y_test, y_pred_proba)
    precision = precision_score(
        y_test, y_pred, labels=classes, average=None, zero_division=0
    )
    recall = recall_score(y_test, y_pred, labels=classes, average=None, zero_division=0)
    f1 = f1_score(y_test, y_pred, labels=classes, average=None, zero_division=0)
    accuracy = accuracy_score(y_test, y_pred)

    overall_report_df = pd.DataFrame(
        {"Value": [roc_auc, precision, recall, f1, pr_auc, accuracy]},
        index=["ROC AUC", "Precision", "Recall", "F1-score", "PR AUC", "Accuracy"],
    )

    overall_report_df.to_csv(f"{output_dir}/overall_report.csv")
