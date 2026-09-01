import numpy as np
import pandas as pd
from sklearn.metrics import (precision_recall_fscore_support, roc_auc_score, precision_recall_curve, roc_curve, auc)
from sklearn.preprocessing import label_binarize

def evaluate(holdout_name, model_name, probs_file, meta_file, classes_file):
    probs = np.load(probs_file)
    meta = pd.read_csv(meta_file)
    with open(classes_file) as f:
        classes = f.read().splitlines()

    y_true = meta["true_label"].values
    y_pred_idx = probs.argmax(axis=1)
    class_to_idx = {c: i for i, c in enumerate(classes)}
    y_pred = np.array([classes[i] for i in y_pred_idx])

    # only evaluate rows the model actually covered (nonzero prob row for graph;
    # VAE always covers everyone it has a fingerprint for)
    covered = probs.sum(axis=1) > 0
    y_true_c, y_pred_c, probs_c = y_true[covered], y_pred[covered], probs[covered]

    accuracy = (y_true_c == y_pred_c).mean()

    # precision/recall/F1 - macro (unweighted across classes - penalizes
    # ignoring rare classes) and weighted (accounts for class size)
    y_true_idx = np.array([class_to_idx.get(t, -1) for t in y_true_c])
    valid = y_true_idx >= 0  # true label must be a valid known class (it will be, but guard anyway)
    p_macro, r_macro, f1_macro, _ = precision_recall_fscore_support(
        y_true_c[valid], y_pred_c[valid], average="macro", zero_division=0)
    p_weighted, r_weighted, f1_weighted, _ = precision_recall_fscore_support(
        y_true_c[valid], y_pred_c[valid], average="weighted", zero_division=0)

    # ROC-AUC: one-vs-rest, macro-averaged. Only over classes that actually
    # appear in this holdout's true labels (an OvR AUC is undefined for a
    # class with zero positive examples to evaluate against)
    present_classes = sorted(set(y_true_c[valid]))
    present_idx = [class_to_idx[c] for c in present_classes]
    y_true_bin = label_binarize(y_true_c[valid], classes=present_classes)
    probs_present = probs_c[valid][:, present_idx]
    try:
        auc_macro = roc_auc_score(y_true_bin, probs_present, average="macro", multi_class="ovr")
    except ValueError as e:
        auc_macro = np.nan
        print(f"  AUC macro failed: {e}")

    print(f"\n{holdout_name} / {model_name}")
    print(f"  Coverage: {covered.sum()} / {len(y_true)}")
    print(f"  Accuracy: {accuracy:.3f}")
    print(f"  Precision (macro/weighted): {p_macro:.3f} / {p_weighted:.3f}")
    print(f"  Recall    (macro/weighted): {r_macro:.3f} / {r_weighted:.3f}")
    print(f"  F1        (macro/weighted): {f1_macro:.3f} / {f1_weighted:.3f}")
    print(f"  ROC-AUC (macro, one-vs-rest, {len(present_classes)} classes present in holdout): {auc_macro:.3f}")

    # micro-average PR and ROC curves (flatten all classes into one binary problem)
    y_true_bin_all = label_binarize(y_true_c[valid], classes=classes)
    precision, recall, _ = precision_recall_curve(y_true_bin_all.ravel(), probs_c[valid].ravel())
    fpr, tpr, _ = roc_curve(y_true_bin_all.ravel(), probs_c[valid].ravel())
    roc_auc_micro = auc(fpr, tpr)

    return {"accuracy": accuracy, "precision_macro": p_macro, "recall_macro": r_macro, "f1_macro": f1_macro,
            "precision_weighted": p_weighted, "recall_weighted": r_weighted, "f1_weighted": f1_weighted,
            "auc_macro_ovr": auc_macro, "auc_micro": roc_auc_micro,
            "pr_curve": (precision, recall), "roc_curve": (fpr, tpr)}


results = {}
results[("Class", "VAE")] = evaluate("Class", "VAE", "vae_probs_class_v2.npy", "vae_holdout_meta_class_v2.csv", "vae_classes_class_v2.txt")
results[("Class", "Graph")] = evaluate("Class", "Graph", "graph_probs_class_v2.npy", "graph_holdout_meta_class_v2.csv", "graph_classes_class_v2.txt")
results[("Pathway", "VAE")] = evaluate("Pathway", "VAE", "vae_probs_pathway_v2.npy", "vae_holdout_meta_pathway_v2.csv", "vae_classes_pathway_v2.txt")
results[("Pathway", "Graph")] = evaluate("Pathway", "Graph", "graph_probs_pathway_v2.npy", "graph_holdout_meta_pathway_v2.csv", "graph_classes_pathway_v2.txt")

summary_rows = []
for (task, model), r in results.items():
    summary_rows.append({"task": task, "model": model, **{k: v for k, v in r.items() if k not in ("pr_curve", "roc_curve")}})
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv("metrics_summary.csv", index=False)
print(summary_df.to_string(index=False))

import pickle
with open("curves.pkl", "wb") as f:
    pickle.dump({k: (v["pr_curve"], v["roc_curve"]) for k, v in results.items()}, f)
