import matplotlib.pyplot as plt
from lists import get_data
from lists import build_exclude_names
import csv
import numpy as np
import seaborn as sns
from sklearn.metrics import accuracy_score, precision_recall_fscore_support, confusion_matrix, roc_curve, auc
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split, GridSearchCV
from collections import defaultdict

metadata = {}
with open("SQ_loci_meta.csv", newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        locus = row['cluster']
        pathway = row['Pathway']
        species = row['Species']
        metadata[locus] = {'pathway': pathway, 'species': species}

exclude_names = build_exclude_names("all_ibs.tsv", "all_hc.tsv")
ibs_data = get_data("all_ibs.tsv", exclude_names)
hc_data = get_data("all_hc.tsv", exclude_names)

all_loci = sorted(set(ibs_data.keys()) | set(hc_data.keys()))
selected_loci = all_loci[:30]

def prepare_vectorized_data(ibs_data, hc_data, loci, threshold=6.5):
    X, y = [], []
    n_ibs_samples = len(next(iter(ibs_data.values())))
    n_hc_samples = len(next(iter(hc_data.values())))

    # IBS
    for i in range(n_ibs_samples):
        sample_vec = []
        for locus in loci:
            vals = ibs_data.get(locus, [0] * n_ibs_samples)
            val = vals[i] if i < len(vals) else 0
            sample_vec.append(val if val > threshold else 0)
        X.append(sample_vec)
        y.append(1)

    # HC
    for i in range(n_hc_samples):
        sample_vec = []
        for locus in loci:
            vals = hc_data.get(locus, [0] * n_hc_samples)
            val = vals[i] if i < len(vals) else 0
            sample_vec.append(val if val > threshold else 0)
        X.append(sample_vec)
        y.append(0)

    return np.array(X), np.array(y)


X, y = prepare_vectorized_data(ibs_data, hc_data, selected_loci, threshold=6.5)

def find_best_test_size(X, y, test_sizes, random_state=42):
    results = []
    rf = RandomForestClassifier(n_estimators=100, random_state=random_state)
    for ts in test_sizes:
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=ts, stratify=y, random_state=random_state)
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state)
        scores = cross_val_score(rf, X_train, y_train, cv=cv, scoring='accuracy')
        mean_score = scores.mean()
        results.append({'test_size': ts, 'cv_mean_accuracy': mean_score})
        print(f"Test size: {ts:.2f}, CV mean accuracy: {mean_score:.4f}")

    best = max(results, key=lambda x: x['cv_mean_accuracy'])
    print(f"\nBest test_size: {best['test_size']:.2f} with CV accuracy: {best['cv_mean_accuracy']:.4f}")
    return best['test_size'], results

test_sizes = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
best_size, all_results = find_best_test_size(X, y, test_sizes)
print(f"Selected best test size: {best_size}")

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=best_size, stratify=y, random_state=42)

param_grid = {
    'n_estimators': [100, 200, 500],
    'max_depth': [None, 10, 30, 50, 100],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4],
    'bootstrap': [True, False],
    'criterion': ['gini', 'entropy']
}

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

rf_model = RandomForestClassifier(random_state=42)
grid_search = GridSearchCV(estimator=rf_model, param_grid=param_grid, cv=cv, n_jobs=-1, scoring='accuracy')

grid_search.fit(X_train, y_train)

best_rf = grid_search.best_estimator_
best_rf.fit(X_train, y_train)

y_pred = best_rf.predict(X_test)
y_proba = best_rf.predict_proba(X_test)[:, 1]

accuracy = accuracy_score(y_test, y_pred)
precision, recall, f1, _ = precision_recall_fscore_support(y_test, y_pred, average='binary')

print(f"Test Accuracy: {accuracy:.3f}")
print(f"Precision: {precision:.3f}")
print(f"Recall: {recall:.3f}")
print(f"F1-score: {f1:.3f}")

cm = confusion_matrix(y_test, y_pred)
cm_percent = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis] * 100

plt.figure(figsize=(6,5))
sns.heatmap(cm_percent, annot=True, fmt='.1f', cmap='Blues', cbar_kws={'format': '%.0f%%'})
plt.title('Confusion Matrix (Percentage)')
plt.xlabel('Predicted Label\n0 = HC, 1 = IBS')
plt.ylabel('True Label\n0 = HC, 1 = IBS')
plt.show()

fpr, tpr, thresholds = roc_curve(y_test, y_proba)
roc_auc = auc(fpr, tpr)
plt.figure(figsize=(6,5))
plt.plot(fpr, tpr, label=f'ROC curve (AUC = {roc_auc:.3f})')
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel('False Positive Rate (HC as negative class)')
plt.ylabel('True Positive Rate (IBS as positive class)')
plt.title('ROC Curve')
plt.legend(loc='lower right')
plt.show()
