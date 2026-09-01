import pandas as pd
import numpy as np

np.random.seed(123)

membership = pd.read_csv("data/processed/metabolite_pathway_membership.csv")
print(f"Metabolite-pathway membership rows: {len(membership)}, distinct metabolites: {membership['entity_id'].nunique()}")

# A metabolite can sit in multiple pathways - for a single-label task (matching
# how Class prediction works) take its FIRST listed pathway as the ground
# truth, same simplification principle used throughout. Multi-label pathway
# assignment would be a legitimate future extension, not done here.
primary_pathway = membership.groupby("entity_id").first().reset_index()[["entity_id", "canonical_metabolite", "Pathway"]]
primary_pathway = primary_pathway.rename(columns={"Pathway": "primary_pathway"})
primary_pathway.to_csv("data/processed/entity_primary_pathway.csv", index=False)
print(f"entity_primary_pathway.csv: {len(primary_pathway)} metabolites")

# dedicated holdout for PATHWAY validation - sampled from the pathway-known
# population specifically (this is NOT the same population as the Class
# holdout - only 645 metabolites have a known pathway at all)
n_holdout = int(0.2 * len(primary_pathway))
pathway_holdout_ids = primary_pathway.sample(n=n_holdout, random_state=123)["entity_id"].tolist()
pd.DataFrame({"entity_id": pathway_holdout_ids}).to_csv("data/processed/pathway_holdout_ids.csv", index=False)
print(f"pathway_holdout_ids.csv: {len(pathway_holdout_ids)} metabolites (20% of the 645 with known pathway)")
