import pickle
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 11

with open("curves.pkl", "rb") as f:
    curves = pickle.load(f)
summary = pd.read_csv("metrics_summary.csv")


fig, axes = plt.subplots(2, 2, figsize=(13, 11))
colors = {"VAE": "#4C72B0", "Graph": "#C44E52"}

for col, task in enumerate(["Class", "Pathway"]):
    ax_pr = axes[0, col]
    ax_roc = axes[1, col]
    for model in ["VAE", "Graph"]:
        (precision, recall), (fpr, tpr) = curves[(task, model)]
        auc_val = summary[(summary.task == task) & (summary.model == model)]["auc_micro"].values[0]
        ax_pr.plot(recall, precision, label=f"{model}", color=colors[model], linewidth=2)
        ax_roc.plot(fpr, tpr, label=f"{model} (AUC={auc_val:.3f})", color=colors[model], linewidth=2)
    ax_roc.plot([0, 1], [0, 1], 'k--', linewidth=0.8, alpha=0.5)
    ax_pr.set_xlabel("Recall"); ax_pr.set_ylabel("Precision")
    ax_pr.set_title(f"{task} - Precision-Recall (micro-avg, one-vs-rest)", fontsize=11, fontweight='bold')
    ax_pr.legend(fontsize=9)
    ax_roc.set_xlabel("False Positive Rate"); ax_roc.set_ylabel("True Positive Rate")
    ax_roc.set_title(f"{task} - ROC (micro-avg, one-vs-rest)", fontsize=11, fontweight='bold')
    ax_roc.legend(fontsize=9)

plt.tight_layout()
plt.savefig("fig6_pr_roc_curves.png", dpi=200, bbox_inches='tight')
plt.close()
print("saved fig6")


fig, axes = plt.subplots(1, 2, figsize=(14, 6))
metrics_to_plot = ["accuracy", "f1_macro", "f1_weighted", "auc_macro_ovr"]
metric_labels = ["Accuracy", "F1 (macro)\n<- rare-class sensitive", "F1 (weighted)", "ROC-AUC\n(macro OvR)"]

for ax, task in zip(axes, ["Class", "Pathway"]):
    x = np.arange(len(metrics_to_plot))
    width = 0.35
    vae_vals = [summary[(summary.task==task)&(summary.model=="VAE")][m].values[0] for m in metrics_to_plot]
    graph_vals = [summary[(summary.task==task)&(summary.model=="Graph")][m].values[0] for m in metrics_to_plot]
    ax.bar(x - width/2, vae_vals, width, label="VAE", color=colors["VAE"])
    ax.bar(x + width/2, graph_vals, width, label="Graph", color=colors["Graph"])
    for i, (v, g) in enumerate(zip(vae_vals, graph_vals)):
        ax.text(i - width/2, v + 0.02, f"{v:.2f}", ha='center', fontsize=8)
        ax.text(i + width/2, g + 0.02, f"{g:.2f}", ha='center', fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(metric_labels, fontsize=9)
    ax.set_ylim(0, 1.08)
    ax.set_title(f"{task} prediction - metric comparison\n(stratified, class-imbalance-aware holdout)", fontsize=11, fontweight='bold')
    ax.legend()

plt.tight_layout()
plt.savefig("fig7_metric_comparison.png", dpi=200, bbox_inches='tight')
plt.close()
print("saved fig7")
