import numpy as np
import pandas as pd
import argparse
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.feature_selection import SelectKBest, mutual_info_classif, RFE
from sklearn.linear_model import LogisticRegression, LassoCV
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.inspection import permutation_importance
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.metrics import make_scorer, f1_score, roc_auc_score
from sklearn.impute import SimpleImputer

# ----- Setup: parse arguments from bash script -------
parser = argparse.ArgumentParser()
parser.add_argument('--Featurematrix', required=True, help='Path to feature matrix')
parser.add_argument('--Target', required=True, help='Path to target matrix')
parser.add_argument('--output_dir', required=True, help='Path to preferred output directory')
args = parser.parse_args()

Featurematrix = args.Featurematrix
Target = args.Target
output_dir = args.output_dir

X_df = pd.read_csv(Featurematrix, na_values=['NA'])
X_df = X_df.loc[:, ~X_df.columns.str.startswith(('chrX','chrY','chrM'))]
X_df = X_df.dropna(axis=1, how='all') # Drop columns that are entirely NA
X_df = X_df.iloc[:, 1:]  # skip first column with sample names
X = X_df.values

y_df = pd.read_csv(Target)
y = y_df.iloc[:, 1].astype(int).values #2nd column is label column

# ----- Setup : define needed classes and functions ------
class LassoSelector(BaseEstimator, TransformerMixin):
    def __init__(self, max_iter=10000, cv=3, random_state=42):
        self.max_iter = max_iter
        self.cv = cv
        self.random_state = random_state

    def fit(self, X, y):
        self.model_ = LassoCV(cv=self.cv, random_state=self.random_state, max_iter=self.max_iter)
        self.model_.fit(X, y)
        self.selected_mask_ = self.model_.coef_ != 0
        return self

    def transform(self, X):
        return X[:, self.selected_mask_]

class PermutationImportanceSelector(BaseEstimator, TransformerMixin):
    def __init__(self, estimator, k=30, scoring='f1', random_state=42):
        self.estimator = estimator
        self.k = k
        self.scoring = scoring
        self.random_state = random_state

    def fit(self, X, y):
        self.estimator_ = clone(self.estimator)
        self.estimator_.fit(X, y)
        result = permutation_importance(
            self.estimator_, X, y,
            scoring=self.scoring,
            n_repeats=10,
            random_state=self.random_state
        )
        self.top_indices_ = np.argsort(result.importances_mean)[::-1][:self.k]
        return self

    def transform(self, X):
        return X[:, self.top_indices_]


def run_nested_cv(X, y, feature_selector, model, param_grid, use_pca=False, n_components=10):
    outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42) #to preserve class distribution
    inner_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42) #to preserve class distribution

    f1_binary = make_scorer(f1_score, average='binary')

    scores_f1 = []
    scores_auc = []

    for train_idx, test_idx in outer_cv.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        steps = [
            ('imputer', SimpleImputer(strategy='mean')),  # Impute missing values per feature using mean
            ('scaler', StandardScaler())
        ]
        if feature_selector:
            steps.append(('feature_selection', feature_selector))
        if use_pca:
            steps.append(('pca', PCA(n_components=n_components)))

        # Add model as a placeholder for now
        steps.append(('model', model))
        pipe = Pipeline(steps)

        # Define parameter grid using model prefix in the pipeline
        search = GridSearchCV(pipe, param_grid=param_grid, cv=inner_cv, scoring=f1_binary, n_jobs=-1)
        search.fit(X_train, y_train)

        y_pred = search.predict(X_test)
        score_f1 = f1_score(y_test, y_pred, average='binary')
        scores_f1.append(score_f1)

        y_proba = search.predict_proba(X_test)[:, 1]
        score_auc = roc_auc_score(y_test, y_proba)
        scores_auc.append(score_auc)

    return np.mean(scores_f1), np.std(scores_f1), np.mean(scores_auc), np.std(scores_auc)

# ----- Main: test all combinations of feature selection methods and ML models -------

def main(X, y):
    results = []

    # Filter Method: Mutual Information
    filter_selector = SelectKBest(mutual_info_classif, k=30)

    # Wrapper: RFE for linear methods, Permuation Importance for SVM with rbf kernel, 
    rfe_lr = RFE(estimator=LogisticRegression(penalty='l2', solver='liblinear', max_iter=500), n_features_to_select=30)
    rfe_svm_lin = RFE(estimator=SVC(kernel='linear'), n_features_to_select=30)
    perm_svm_rbf = PermutationImportanceSelector(SVC(kernel='rbf', probability=True), k=30)
    rfe_rf = RFE(estimator=RandomForestClassifier(n_estimators=100), n_features_to_select=30)
    rfe_xgb = RFE(estimator=XGBClassifier(eval_metric='logloss'), n_features_to_select=30)


    # Embedded: Lasso
    embedded_selector = LassoSelector()

    # Models and their parameter grids
    models_and_params = {
        'LogisticRegression': (
            LogisticRegression(max_iter=1000),
            {'model__C': [0.01, 0.1, 1, 10]}
        ),
        'SVM_Linear': (
            SVC(kernel='linear', probability=True),
            {'model__C': [0.1, 1, 10]}
        ),
        'SVM_NonLinear': (
            SVC(kernel='rbf', probability=True),
            {'model__C': [0.1, 1, 10], 'model__gamma': ['scale', 'auto']}
        ),
        'RandomForest': (
            RandomForestClassifier(n_estimators=100),
            {'model__max_depth': [None, 5, 10]}
        ),
        'XGBoost': (
            XGBClassifier(eval_metric='logloss'),
            {'model__max_depth': [3, 5, 7], 'model__learning_rate': [0.01, 0.1]}
        )
    }

    selectors = {
        'Filter': filter_selector,
        'Wrapper_LR': rfe_lr,
        'Wrapper_RF': rfe_rf,
        'Wrapper_SVM_Lin': rfe_svm_lin,
        'Wrapper_SVM_RFB': perm_svm_rbf,
        'Wrapper_XGB': rfe_xgb,
        'Embedded_Lasso': embedded_selector
    }

    for sel_name, selector in selectors.items():
        for model_name, (model, param_grid) in models_and_params.items():
            # Skip incompatible model-selector pairs
            if sel_name == 'Wrapper_LR' and model_name != 'LogisticRegression':
                continue
            if sel_name == 'Wrapper_RF' and model_name != 'RandomForest':
                continue
            if sel_name == 'Wrapper_SVM_Lin' and model_name != 'SVM_Linear':
                continue
            if sel_name == 'Wrapper_SVM_RFB' and model_name != 'SVM_NonLinear':
                continue
            if sel_name == 'Wrapper_XGB' and model_name != 'XGBoost':
                continue
            
            print(f"Running feature selector: {sel_name}, model: {model_name}, PCA: No PCA", flush=True)
            f1_mean, f1_std, auc_mean, auc_std = run_nested_cv(X, y, selector, model, param_grid, use_pca=False)
            results.append((sel_name, model_name, 'No PCA', f1_mean, f1_std, auc_mean, auc_std))
            print(f"Completed feature selector: {sel_name}, model: {model_name}, PCA: No PCA", flush=True)

            print(f"Running feature selector: {sel_name}, model: {model_name}, PCA: With PCA", flush=True)
            f1_mean, f1_std, auc_mean, auc_std = run_nested_cv(X, y, selector, model, param_grid, use_pca=True, n_components=10)
            results.append((sel_name, model_name, 'With PCA', f1_mean, f1_std, auc_mean, auc_std))
            print(f"Completed feature selector: {sel_name}, model: {model_name}, PCA: With PCA", flush=True)

    results_df = pd.DataFrame(results, columns=['Selector', 'Model', 'PCA', 'Mean F1', 'Std F1', "Mean AUC", "Std AUC"])
    return results_df

# ----- Run -----
results_df = main(X, y)
results_df.to_csv(f"{output_dir}/model_selection_results.csv", index=False)
print(f"Results saved to {output_dir}/model_selection_results.csv")