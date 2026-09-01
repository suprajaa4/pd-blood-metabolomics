import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
with open("curves.pkl", "rb") as f:
    curves = pickle.load(f)

summary = pd.read_csv("metrics_summary.csv")

colors = {
    "VAE": "#4C72B0",
    "Graph": "#C44E52"
}

tasks = ["Class", "Pathway"]
models = ["VAE", "Graph"]


fig, axes = plt.subplots(2, 2, figsize=(13, 11))

for i, task in enumerate(tasks):

    pr_ax = axes[0, i]
    roc_ax = axes[1, i]

    for model in models:

        (precision, recall), (fpr, tpr) = curves[(task, model)]

        row = summary[
            (summary["task"] == task) &
            (summary["model"] == model)
        ].iloc[0]

        auc = row["auc_micro"]

        pr_ax.plot(
            recall,
            precision,
            label=model,
            color=colors[model]
        )

        roc_ax.plot(
            fpr,
            tpr,
            label=f"{model} (AUC={auc:.3f})",
            color=colors[model]
        )

    # PR plot
    pr_ax.set_title(f"{task} - Precision-Recall")
    pr_ax.set_xlabel("Recall")
    pr_ax.set_ylabel("Precision")
    pr_ax.legend()

    # ROC plot
    roc_ax.plot([0, 1], [0, 1], "k--")
    roc_ax.set_title(f"{task} - ROC")
    roc_ax.set_xlabel("False Positive Rate")
    roc_ax.set_ylabel("True Positive Rate")
    roc_ax.legend()


plt.tight_layout()
plt.savefig("fig6_pr_roc_curves.png", dpi=200)
plt.close()

print("Saved fig6")



metrics = [
    "accuracy",
    "f1_macro",
    "f1_weighted",
    "auc_macro_ovr"
]

labels = [
    "Accuracy",
    "F1 Macro",
    "F1 Weighted",
    "ROC-AUC"
]

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, task in zip(axes, tasks):

    task_data = summary[summary["task"] == task]

    vae = task_data[task_data["model"] == "VAE"].iloc[0]
    graph = task_data[task_data["model"] == "Graph"].iloc[0]

    vae_values = vae[metrics].values
    graph_values = graph[metrics].values

    x = np.arange(len(metrics))
    width = 0.35

    ax.bar(
        x - width / 2,
        vae_values,
        width,
        label="VAE",
        color=colors["VAE"]
    )

    ax.bar(
        x + width / 2,
        graph_values,
        width,
        label="Graph",
        color=colors["Graph"]
    )

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Score")
    ax.set_title(f"{task} Prediction")
    ax.legend()


plt.tight_layout()
plt.savefig("fig7_metric_comparison.png", dpi=200)
plt.close()

print("Saved fig7")
