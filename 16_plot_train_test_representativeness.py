"""
Does the holdout (test) set's class balance actually mirror the training
set's? For each class, compares its proportion of the training set against
its proportion of the holdout set - a stratified split should land close to
the diagonal; a poorly-representative one (like the original random holdout)
would scatter far off it, especially for rare classes.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, entropy

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 11

entities = pd.read_csv("data/processed/entities_full.csv")
primary_pathway = pd.read_csv("data/processed/entity_primary_pathway.csv")
entities = entities.merge(primary_pathway[["entity_id", "primary_pathway"]], on="entity_id", how="left")

class_train = pd.read_csv("data/processed/class_training_set_v2.csv")
pathway_train = pd.read_csv("data/processed/pathway_training_set_v2.csv")
class_holdout_ids = set(pd.read_csv("data/processed/holdout_ids_v2_stratified.csv")["entity_id"])
pathway_holdout_ids = set(pd.read_csv("data/processed/pathway_holdout_ids_v2_stratified.csv")["entity_id"])


def representativeness(label_col, train_df, holdout_ids, task_name, ax):
    holdout_df = entities[entities["entity_id"].isin(holdout_ids)]

    train_counts = train_df[label_col].value_counts()
    holdout_counts = holdout_df[label_col].value_counts()

    all_classes = sorted(set(train_counts.index) | set(holdout_counts.index))
    train_prop = (train_counts.reindex(all_classes, fill_value=0) / train_counts.sum())
    holdout_prop = (holdout_counts.reindex(all_classes, fill_value=0) / holdout_counts.sum())

    # log-log scatter (add a small floor so zero-count classes are still visible, at the axis edge)
    floor = 0.5 / max(train_counts.sum(), holdout_counts.sum())
    x = np.log10(train_prop.values + floor)
    y = np.log10(holdout_prop.values + floor)
    sizes = np.log1p(train_counts.reindex(all_classes, fill_value=0).values) * 15 + 10

    ax.scatter(x, y, s=sizes, alpha=0.55, c="#4C72B0", edgecolor="white", linewidth=0.4)
    lims = [min(x.min(), y.min()) - 0.2, max(x.max(), y.max()) + 0.2]
    ax.plot(lims, lims, 'k--', linewidth=1, alpha=0.6, label="Perfect representativeness (y=x)")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel("log10(proportion of TRAINING set)")
    ax.set_ylabel("log10(proportion of TEST/holdout set)")

    # correlation + JSD as quantitative summaries
    both_present = (train_prop > 0) & (holdout_prop > 0)
    r, _ = pearsonr(train_prop[both_present], holdout_prop[both_present])
    m = 0.5 * (train_prop.values + holdout_prop.values)
    m_safe = np.where(m > 0, m, 1)
    jsd = 0.5 * entropy(train_prop.values, m_safe) + 0.5 * entropy(holdout_prop.values, m_safe)

    n_zero_holdout = (holdout_counts.reindex(all_classes, fill_value=0) == 0).sum()
    ax.set_title(f"{task_name}\nPearson r={r:.3f} | JSD={jsd:.4f} | "
                 f"{n_zero_holdout}/{len(all_classes)} classes absent from test set",
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='upper left')
    return r, jsd, n_zero_holdout


fig, axes = plt.subplots(1, 2, figsize=(13, 6.5))
r_c, jsd_c, miss_c = representativeness("Class", class_train, class_holdout_ids, "Class", axes[0])
r_p, jsd_p, miss_p = representativeness("primary_pathway", pathway_train, pathway_holdout_ids, "Pathway", axes[1])

plt.suptitle("How representative is the test set of the training set's class balance?\n"
             "(each point = one class; point size = how common that class is in training)", fontsize=12, y=1.06)
plt.tight_layout()
plt.savefig("results/model_comparison/fig8_train_test_representativeness.png", dpi=200, bbox_inches='tight')
plt.close()

print(f"Class:   Pearson r={r_c:.3f}, JSD={jsd_c:.4f}, {miss_c} classes with zero test examples")
print(f"Pathway: Pearson r={r_p:.3f}, JSD={jsd_p:.4f}, {miss_p} classes with zero test examples")
print("saved fig8_train_test_representativeness.png")
