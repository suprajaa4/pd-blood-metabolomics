import numpy as np
import pandas as pd
DATA_DIR = "data/processed/"
OUTPUT_DIR = "results/graph/"
entities = pd.read_csv(DATA_DIR + "entities_full.csv")
primary_pathway = pd.read_csv(DATA_DIR + "entity_primary_pathway.csv")
entities = entities.merge(primary_pathway[["entity_id", "primary_pathway"]], on="entity_id", how="left")
pd_scores = pd.read_csv(DATA_DIR + "entity_pd_scores.csv")
entities = entities.merge(
    pd_scores[["entity_id", "PD_score", "consistently_altered", "majority_direction", "is_known_drug",
               "annotation_score", "has_formula", "has_SMILES", "has_HMDB", "has_KEGG", "has_PubChem"]],
    on="entity_id", how="left")
print(f"Merged PD_score - {entities['PD_score'].notna().sum()} / {len(entities)} entities scored, "
      f"{entities['consistently_altered'].sum():.0f} consistently altered")

class_holdout_ids = set(pd.read_csv(DATA_DIR + "holdout_ids_v2_stratified.csv")["entity_id"])
pathway_holdout_ids = set(pd.read_csv(DATA_DIR + "pathway_holdout_ids_v2_stratified.csv")["entity_id"])

LAYER_FILES = {
    "study": "edges_study.csv", "adduct": "edges_adduct.csv", "platform": "edges_platform.csv",
    "structural": "edges_structural.csv", "pathway": "edges_pathway.csv", "class": "edges_class.csv",
    "mz": "edges_mz.csv",
}
_layer_cache = {}

def build_edge_table(relations):
    key = tuple(sorted(relations))
    if key in _layer_cache:
        return _layer_cache[key]
    dfs = [pd.read_csv(DATA_DIR + LAYER_FILES[r]) for r in relations]
    combined = pd.concat(dfs, ignore_index=True)
    combined_sym = pd.concat([combined, combined.rename(columns={"from": "to", "to": "from"})], ignore_index=True)
    edge_table = combined_sym.groupby(["from", "to"], as_index=False)["weight"].sum()
    _layer_cache[key] = edge_table
    return edge_table


def predict(target_ids, known_label_lookup, edge_table):
    targets_df = pd.DataFrame({"to": target_ids})
    nbr_edges = edge_table.merge(targets_df, on="to")
    nbr_edges["neighbor_label"] = nbr_edges["from"].map(known_label_lookup)
    nbr_edges = nbr_edges.dropna(subset=["neighbor_label"])
    votes = nbr_edges.groupby(["to", "neighbor_label"])["weight"].sum().reset_index(name="vote")
    totals = votes.groupby("to")["vote"].sum().reset_index(name="total_vote")
    best = votes.loc[votes.groupby("to")["vote"].idxmax()].merge(totals, on="to")
    best["confidence"] = best["vote"] / best["total_vote"]
    best = best.rename(columns={"to": "entity_id", "neighbor_label": "predicted_label"})
    return pd.DataFrame({"entity_id": target_ids}).merge(
        best[["entity_id", "predicted_label", "confidence"]], on="entity_id", how="left")


def run_task(task_name, label_col, holdout_ids, exclude_layer, unlabeled_mask):
    relations = [r for r in LAYER_FILES if r != exclude_layer]
    edge_table = build_edge_table(relations)
    print(f"\n{task_name}  (layers used: {relations})\n")

    full_lookup = dict(zip(entities["entity_id"], entities[label_col]))
    masked_lookup = full_lookup.copy()
    for eid in holdout_ids:
        masked_lookup.pop(eid, None)

    holdout_preds = predict(list(holdout_ids), masked_lookup, edge_table)
    holdout_preds["true_label"] = holdout_preds["entity_id"].map(full_lookup)
    holdout_preds = holdout_preds.dropna(subset=["predicted_label"])
    correct = holdout_preds["predicted_label"] == holdout_preds["true_label"]
    print(f"Holdout validation (n={len(holdout_preds)}/{len(holdout_ids)} covered)")
    print(f"Top-1 accuracy: {correct.mean():.3f}")

    unlabeled_ids = entities.loc[unlabeled_mask, "entity_id"].tolist()
    real_preds = predict(unlabeled_ids, full_lookup, edge_table)
    real_preds = real_preds.merge(entities[["entity_id", "canonical_metabolite"]], on="entity_id")
    real_preds = real_preds.rename(columns={"predicted_label": f"predicted_{label_col}"})
    print(f"Predictions for {len(unlabeled_ids)} metabolites with no known {label_col}: "
          f"coverage {real_preds[f'predicted_{label_col}'].notna().sum()} / {len(unlabeled_ids)}")
    return holdout_preds, real_preds


class_holdout_results, class_real_preds = run_task(
    "CLASS PREDICTION", "Class", class_holdout_ids, "class", entities["Class"].isna())
class_real_preds.to_csv(OUTPUT_DIR + "predicted_class_results.csv", index=False)

pathway_holdout_results, pathway_real_preds = run_task(
    "PATHWAY PREDICTION", "primary_pathway", pathway_holdout_ids, "pathway", entities["primary_pathway"].isna())
pathway_real_preds.to_csv(OUTPUT_DIR + "predicted_pathway_results.csv", index=False)



#  which poorly-annotated metabolites are central to PD-associated subnetworks?
print(f"\npoorly-annotated metabolites central to the PD subnetwork\n")

pd_subnet_ids = set(entities.loc[entities["consistently_altered"] == True, "entity_id"])
full_edges = build_edge_table(list(LAYER_FILES.keys()))
sub_edges = full_edges[full_edges["from"].isin(pd_subnet_ids) & full_edges["to"].isin(pd_subnet_ids)]

G = nx.Graph()
G.add_nodes_from(pd_subnet_ids)
for _, row in sub_edges.iterrows():
    G.add_edge(row["from"], row["to"], weight=row["weight"])
# betweenness needs DISTANCE not similarity - invert the weight
for u, v, d in G.edges(data=True):
    d["distance"] = 1 / (d["weight"] + 1e-6)

print(f"PD subnetwork: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
betweenness = nx.betweenness_centrality(G, weight="distance", normalized=True)

hub_df = pd.DataFrame({"entity_id": list(betweenness.keys()), "betweenness": list(betweenness.values())})
hub_df = hub_df.merge(entities[["entity_id", "canonical_metabolite", "PD_score", "annotation_score",
                                 "has_formula", "has_SMILES", "has_HMDB", "has_KEGG", "has_PubChem", "Class"]],
                       on="entity_id", how="left")
hub_df["degree"] = hub_df["entity_id"].map(dict(G.degree()))
hub_df = hub_df[hub_df["degree"] > 0].sort_values("betweenness", ascending=False)

print("\nMost central metabolites in the PD subnetwork (any annotation level):")
print(hub_df[["canonical_metabolite", "betweenness", "degree", "PD_score", "annotation_score"]].head(10).to_string(index=False))

hidden_hubs = hub_df[hub_df["annotation_score"] <= 2].sort_values("betweenness", ascending=False)
print("\n'Hidden hub' candidates - HIGH centrality, LOW annotation (annotate these first):")
print(hidden_hubs[["canonical_metabolite", "betweenness", "degree", "PD_score", "annotation_score",
                    "has_formula", "has_SMILES", "has_HMDB", "has_KEGG", "has_PubChem", "Class"]]
      .head(10).to_string(index=False))
hidden_hubs.to_csv(OUTPUT_DIR + "hidden_hubs.csv", index=False)


# which metabolites connect multiple dysregulated pathways?

print(f"\nmetabolites connecting multiple dysregulated pathways\n")

membership = pd.read_csv(DATA_DIR + "metabolite_pathway_membership.csv")
altered_membership = membership[membership["entity_id"].isin(pd_subnet_ids)]
pathway_count = altered_membership.groupby("entity_id")["Pathway"].nunique().reset_index(name="n_pathways")
pathway_count = pathway_count.merge(entities[["entity_id", "canonical_metabolite", "PD_score"]], on="entity_id")
pathway_count = pathway_count.sort_values("n_pathways", ascending=False)

print("(a) Direct membership: consistently-altered metabolites in the MOST distinct pathways:")
print(pathway_count.head(10).to_string(index=False))

pathway_edges = build_edge_table(["pathway"])
pathway_sub_edges = pathway_edges[pathway_edges["from"].isin(pd_subnet_ids) & pathway_edges["to"].isin(pd_subnet_ids)]
G_pathway = nx.Graph()
G_pathway.add_nodes_from(pd_subnet_ids)
for _, row in pathway_sub_edges.iterrows():
    G_pathway.add_edge(row["from"], row["to"], weight=row["weight"])
for u, v, d in G_pathway.edges(data=True):
    d["distance"] = 1 / (d["weight"] + 1e-6)

bridge_betweenness = nx.betweenness_centrality(G_pathway, weight="distance", normalized=True)
bridge_df = pd.DataFrame({"entity_id": list(bridge_betweenness.keys()), "pathway_betweenness": list(bridge_betweenness.values())})
bridge_df = bridge_df.merge(entities[["entity_id", "canonical_metabolite", "PD_score"]], on="entity_id")
bridge_df["pathway_degree"] = bridge_df["entity_id"].map(dict(G_pathway.degree()))
bridge_df = bridge_df[bridge_df["pathway_degree"] > 0].sort_values("pathway_betweenness", ascending=False)

print("\n(b) Network bridging : highest betweenness WITHIN the pathway-only layer:")
print(bridge_df[["canonical_metabolite", "pathway_betweenness", "pathway_degree", "PD_score"]].head(10).to_string(index=False))

pathway_count.to_csv(OUTPUT_DIR + "pathway_membership_count.csv", index=False)
bridge_df.to_csv(OUTPUT_DIR + "pathway_bridging.csv", index=False)



#  which chemical classes contain the highest-confidence PD_score metabolites?

print(f"\n classes with the highest-confidence PD association\n")

class_ranking = entities[entities["consistently_altered"] == True].dropna(subset=["Class"]).groupby("Class").agg(
    n_metabolites=("entity_id", "size"),
    mean_PD_score=("PD_score", "mean"),
    max_PD_score=("PD_score", "max"),
).reset_index()
class_ranking["top_metabolite"] = class_ranking["Class"].apply(
    lambda c: entities[(entities["Class"] == c) & (entities["consistently_altered"] == True)]
    .sort_values("PD_score", ascending=False)["canonical_metabolite"].iloc[0])
class_ranking = class_ranking[class_ranking["n_metabolites"] >= 2].sort_values("mean_PD_score", ascending=False)

print(class_ranking.head(10).to_string(index=False))
class_ranking.to_csv(OUTPUT_DIR + "class_ranking.csv", index=False)

# functionally similar metabolites
print(f"\nfunctionally similar metabolites (pathway + structural)\n")

functional_edges = build_edge_table(["pathway", "structural"])
functional_sub = functional_edges[functional_edges["from"].isin(pd_subnet_ids) | functional_edges["to"].isin(pd_subnet_ids)]

name_lookup = dict(zip(entities["entity_id"], entities["canonical_metabolite"]))
functional_sub = functional_sub.copy()
functional_sub["from_name"] = functional_sub["from"].map(name_lookup)
functional_sub["to_name"] = functional_sub["to"].map(name_lookup)
functional_sub = functional_sub.sort_values("weight", ascending=False)

print("top 15 strongest functional-similarity connections (pathway-sharing + structural), among consistently-altered metabolites:")
print(functional_sub[["from_name", "to_name", "weight"]].head(15).to_string(index=False))
functional_sub.to_csv(OUTPUT_DIR + "functional_similarity_pairs.csv", index=False)

