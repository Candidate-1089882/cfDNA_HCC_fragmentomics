import numpy as np
import pandas as pd
import argparse
from sklearn.model_selection import StratifiedKFold, RandomizedSearchCV
from sklearn.feature_selection import SelectKBest, mutual_info_classif, RFE
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.metrics import f1_score, roc_auc_score, make_scorer
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
import numpy as np
from sklearn.metrics import f1_score

# ----- Setup: parse arguments from bash script -------
parser = argparse.ArgumentParser()
parser.add_argument('--Featurematrix', required=True, help='Path to feature matrix')
parser.add_argument('--Target', required=True, help='Path to target matrix')
parser.add_argument('--ValidationFeatures', required=True, help='Path to validation feature matrix')
parser.add_argument('--ValidationTarget', required=True, help='Path to validation target matrix')
parser.add_argument('--output_dir', required=True, help='Path to preferred output directory')
args = parser.parse_args()

Featurematrix = args.Featurematrix
Target = args.Target
ValidationFeatures = args.ValidationFeatures
ValidationTarget = args.ValidationTarget
output_dir = args.output_dir

# Load data
X_df = pd.read_csv(Featurematrix, na_values=['NA'])
X_df = X_df.loc[:, ~X_df.columns.str.startswith(('chrX','chrY','chrM'))]
X_df = X_df.dropna(axis=1, how='all') # Drop columns entirely NA
X_df = X_df.iloc[:, 1:]  # skip first col with sample names
X = X_df.values

y_df = pd.read_csv(Target)
y = y_df.iloc[:, 1].astype(int).values

common_columns = X_df.columns
X_val_df = pd.read_csv(ValidationFeatures, na_values=['NA'])
X_val_df = X_val_df.loc[:, ~X_val_df.columns.str.startswith(('chrX','chrY','chrM'))]
X_val_df = X_val_df.dropna(axis=1, how='all')
X_val_df = X_val_df.iloc[:, 1:]  # skip first col
X_val_df = X_val_df.reindex(columns=common_columns)  # enforce same feature order
X_val = X_val_df.values

y_val_df = pd.read_csv(ValidationTarget)
y_val = y_val_df.iloc[:, 1].astype(int).values

# Feature selectors
filter_selector = SelectKBest(mutual_info_classif, k=30)
rfe_rf = RFE(estimator=RandomForestClassifier(n_estimators=100), n_features_to_select=30)
rfe_xgb = RFE(estimator=XGBClassifier(eval_metric='logloss'), n_features_to_select=30)

# Define selected model-selector-PCA combos
model_selector_pca_combos = [
    ('XGBoost', XGBClassifier(eval_metric='logloss'), rfe_xgb, False),  # Wrapper + no PCA
    ('RandomForest', RandomForestClassifier(n_estimators=100, random_state=42), rfe_rf, False),  # Wrapper + no PCA
    ('LogisticRegression', LogisticRegression(max_iter=1000), filter_selector,  True),  # Filter + PCA
    ('SVM_Linear', SVC(kernel='linear', probability=True, random_state=42), filter_selector,  True),  # Filter + PCA
]

# Parameter distributions for RandomizedSearchCV
param_distributions = {
    'XGBoost': {
        'model__n_estimators': [100, 200, 300],
        'model__learning_rate': [0.01, 0.05, 0.1, 0.2],
        'model__max_depth': [3, 5, 7, 10],
        'model__subsample': [0.6, 0.8, 1.0],
        'model__colsample_bytree': [0.6, 0.8, 1.0],
        'model__gamma': [0, 1, 5],
    },
    'RandomForest': {
        'model__n_estimators': [100, 200, 300],
        'model__max_depth': [None, 5, 10, 20, 30],
        'model__min_samples_split': [2, 5, 10],
        'model__min_samples_leaf': [1, 2, 4],
        'model__bootstrap': [True, False]
    },
    'LogisticRegression': {
        'model__C': np.logspace(-3, 2, 10),
        'model__penalty': ['l2'], 
        'model__solver': ['liblinear']
    },
    'SVM_Linear': {
        'model__C': np.logspace(-3, 2, 10)
    }
}

def find_best_f1_threshold(y_true, y_proba):
    best_threshold = 0.5
    best_f1 = 0
    thresholds = np.linspace(0.1, 0.9, 81)  # steps of 0.01 between 0.1 and 0.9

    for thresh in thresholds:
        y_pred_thresh = (y_proba >= thresh).astype(int)
        score = f1_score(y_true, y_pred_thresh)
        if score > best_f1:
            best_f1 = score
            best_threshold = thresh
    return best_threshold, best_f1

def weighted_score(y_true, y_pred_proba, **kwargs):
    y_pred = (y_pred_proba >= 0.5).astype(int)
    f1 = f1_score(y_true, y_pred)
    auc = roc_auc_score(y_true, y_pred_proba)
    return 0.3 * f1 + 0.7 * auc 

def train_and_validate(X_train, y_train, X_val, y_val, selector, model, param_dist, use_pca):
    print(f"\nStarting for model {model.__class__.__name__} with selector {selector} and PCA={use_pca}")

    # Split training data into internal train/test split
    X_train_split, X_test_split, y_train_split, y_test_split = train_test_split(
        X_train, y_train, test_size=0.2, stratify=y_train, random_state=42)
    
    # Pipeline setup as before
    steps = [
        ('imputer', SimpleImputer(strategy='mean')),
        ('scaler', StandardScaler())
    ]
    if selector is not None:
        steps.append(('feature_selection', selector))
    if use_pca:
        steps.append(('pca', PCA(n_components=0.95)))
    steps.append(('model', model))
    pipe = Pipeline(steps)

    # Fit on training split
    pipe.fit(X_train_split, y_train_split)

    # Predict on internal test split
    y_test_pred = pipe.predict(X_test_split)
    y_test_proba = pipe.predict_proba(X_test_split)[:, 1] if hasattr(pipe, "predict_proba") else None

    print("\nInternal test split performance:")
    print(classification_report(y_test_split, y_test_pred, digits=4))

    # Now do hyperparameter tuning on full training set (as before)
    cv_inner = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
    scorer = make_scorer(weighted_score, needs_proba=True)

    search = RandomizedSearchCV(
        pipe,
        param_distributions=param_dist,
        n_iter=30,
        cv=cv_inner,
        scoring=scorer,
        n_jobs=-1,
        random_state=42,
        verbose=1
    )

    print(f"\nStarting hyperparameter tuning for model {model.__class__.__name__} with selector {selector} and PCA={use_pca}")
    search.fit(X_train, y_train)
    print("Best params:", search.best_params_)
    best_model = search.best_estimator_

    #Check for overfitting
    # Rebuild full pipeline with best parameters from search
    best_pipe = Pipeline(steps)
    best_pipe = clone(search.best_estimator_)  # from sklearn.base import clone
    best_pipe.fit(X_train_split, y_train_split)
    y_test_pred_best = best_pipe.predict(X_test_split)
    y_test_proba_best = best_pipe.predict_proba(X_test_split)[:, 1]

    print("\nInternal test split performance WITH BEST PARAMETERS:")
    print(classification_report(y_test_split, y_test_pred_best, digits=4))

    #Check selected features
    feature_df = None
    if use_pca:
        # Extract PCA loadings
        pca = best_model.named_steps['pca']
        feature_selector = best_model.named_steps['feature_selection']
        selected_mask = feature_selector.get_support()
        selected_features = common_columns[selected_mask]

        loadings = pd.DataFrame(
            pca.components_.T,
            index=selected_features,
            columns=[f'PC{i+1}' for i in range(pca.n_components_)]
        )
        
        # Add magnitude of loading for PC1 as a proxy importance
        loadings['PC1_abs_loading'] = loadings['PC1'].abs()
        feature_df = loadings.sort_values(by='PC1_abs_loading', ascending=False).drop(columns='PC1_abs_loading')

        print("\nPCA loadings extracted as feature contributions.")
    else:
        # Extract model feature importances or coefficients
        feature_selector = best_model.named_steps['feature_selection']
        selected_mask = feature_selector.get_support()
        selected_features = common_columns[selected_mask]
        model_step = best_model.named_steps['model']

        if hasattr(model_step, 'feature_importances_'):
            importances = model_step.feature_importances_
            feature_df = pd.DataFrame({
                "feature": selected_features,
                "importance": importances
            }).sort_values(by="importance", ascending=False)

        elif hasattr(model_step, 'coef_'):
            coefs = model_step.coef_.flatten()
            feature_df = pd.DataFrame({
                "feature": selected_features,
                "importance": np.abs(coefs),
                "coefficient": coefs
            }).sort_values(by="importance", ascending=False)

        print("\nModel-based feature importances extracted.")

    # Save feature importances
    importance_path = f"{output_dir}/{model.__class__.__name__}_{type(selector).__name__}_PCA{use_pca}_feature_importances_C.csv"
    feature_df.to_csv(importance_path, index=False)
    print(f"\nFeature importances saved to: {importance_path}") 

    # Evaluate on validation set
    if hasattr(best_model, "predict_proba"):
        y_val_proba = best_model.predict_proba(X_val)[:, 1]
        best_thresh, best_f1 = find_best_f1_threshold(y_val, y_val_proba)
        print(f"\nBest F1 threshold on validation set: {best_thresh:.3f} with F1={best_f1:.4f}")

        y_val_pred = (y_val_proba >= best_thresh).astype(int)
    else:
        y_val_pred = best_model.predict(X_val)
        best_thresh = 0.5
        best_f1 = f1_score(y_val, y_val_pred)
    auc = roc_auc_score(y_val, y_val_proba) if hasattr(best_model, "predict_proba") else None

    print(f"\nValidation set performance with threshold {best_thresh:.3f}:")
    print(classification_report(y_val, y_val_pred, digits=4))

    # Save predictions
    preds_df = pd.DataFrame({
        'y_true': y_val,
        'y_pred': y_val_pred,
        'y_proba': y_val_proba if hasattr(best_model, "predict_proba") else [None] * len(y_val)
    })
    preds_df.to_csv(f"{output_dir}/{model.__class__.__name__}_{type(selector).__name__}_PCA{use_pca}_predictions_C.csv", index=False)

    return search.best_params_, best_f1, auc, best_thresh


def main():
    records = []
    for model_name, model, selector, use_pca in model_selector_pca_combos:
        # pick selector
        selector = selector
        params = param_distributions[model_name]

        best_params, val_f1, val_auc, best_tresh = train_and_validate(X, y, X_val, y_val, clone(selector), clone(model), params, use_pca)
        records.append({
            'Model': model_name,
            'Selector': selector,
            'PCA': use_pca,
            'BestParams': best_params,
            'Validation_F1': val_f1,
            'Validation_AUC': val_auc,
            'F1 Treshold' : best_tresh
        })

    results_df = pd.DataFrame(records)
    results_df.to_csv(f"{output_dir}/validation_results_C.csv", index=False)
    print(f"Validation results saved to {output_dir}/validation_results_C.csv")

if __name__ == "__main__":
    main()
