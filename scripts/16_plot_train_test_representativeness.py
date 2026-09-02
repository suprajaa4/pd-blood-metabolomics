import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

entities = pd.read_csv("data/processed/entities_full.csv")
primary_pathway = pd.read_csv("data/processed/entity_primary_pathway.csv")
entities = entities.merge(
    primary_pathway[["entity_id", "primary_pathway"]],
    on="entity_id",
    how="left"
)

class_train = pd.read_csv("data/processed/class_training_set_v2.csv")
pathway_train = pd.read_csv("data/processed/pathway_training_set_v2.csv")
class_holdout_ids = set(pd.read_csv("data/processed/holdout_ids_v2_stratified.csv")["entity_id"])
pathway_holdout_ids = set(pd.read_csv("data/processed/pathway_holdout_ids_v2_stratified.csv")["entity_id"])
class_test = entities[entities["entity_id"].isin(class_holdout_ids)]
pathway_test = entities[entities["entity_id"].isin(pathway_holdout_ids)]
class_train_labels = set(class_train["Class"].dropna())
class_test_labels = set(class_test["Class"].dropna())
pathway_train_labels = set(pathway_train["primary_pathway"].dropna())
pathway_test_labels = set(pathway_test["primary_pathway"].dropna())

fig, axes = plt.subplots(1, 2, figsize=(12, 6))

venn2(
    [class_train_labels, class_test_labels],
    set_labels=("Train labels", "Test labels"),
    ax=axes[0]
)
axes[0].set_title("Class labels")

venn2(
    [pathway_train_labels, pathway_test_labels],
    set_labels=("Train labels", "Test labels"),
    ax=axes[1]
)
axes[1].set_title("Pathway labels")
plt.tight_layout()
plt.savefig("results/model_comparison/fig8_train_test_label_venn.png", dpi=200, bbox_inches="tight")
plt.close()
print("class overlap:", len(class_train_labels & class_test_labels))
print("pathway overlap:", len(pathway_train_labels & pathway_test_labels))
