import numpy as np
import pandas as pd
from sklearn.metrics import (  # type: ignore
    accuracy_score,
    average_precision_score,
    f1_score,
    roc_auc_score,
)
from sklearn.model_selection import (  # type: ignore
    RandomizedSearchCV,  # Import from sklearn, not xgboost
    StratifiedKFold,
    train_test_split,
)
from xgboost import XGBClassifier

from python.normalization import presence_absence_transform, rclr_transform


def prevalence_filtering(
    count_table_df: pd.DataFrame, prevalence_threshold: float
) -> pd.DataFrame:
    """Remove features from feature table with prevalence below threshold."""
    feature_prevalence = (
        (count_table_df > 0).sum(axis=1) / count_table_df.shape[1] * 100
    )

    # Remove features with prevalence below threshold
    return count_table_df.loc[feature_prevalence >= prevalence_threshold, :]


def train_xgboost_model(
    count_table_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    primary_label: str,
    normalization_method: str,
) -> tuple[XGBClassifier, dict]:
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
    else:
        raise ValueError(
            f"Unsupported normalization method: {normalization_method}. "
            "Choose 'rclr' or 'presence_absence'."
        )

    # Extract target variable from metadata
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
    search = RandomizedSearchCV(
        estimator=base_model,
        param_distributions=param_dist,
        n_iter=60,
        scoring="roc_auc",
        n_jobs=-1,
        cv=cv,
        refit=True,
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

    # Evaluate on test set
    y_pred = best_model.predict(X_test)
    y_pred_proba = best_model.predict_proba(X_test)[:, 1]

    # Calculate metrics
    roc_auc = roc_auc_score(y_test, y_pred_proba)
    pr_auc = average_precision_score(y_test, y_pred_proba)
    accuracy = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)

    scores = {
        "roc_auc": roc_auc,
        "pr_auc": pr_auc,
        "accuracy": accuracy,
        "f1_score": f1,
    }

    return best_model, scores
