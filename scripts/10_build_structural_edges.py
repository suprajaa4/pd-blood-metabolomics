import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

entities = pd.read_csv("data/processed/entities_full.csv")
smiles_df = entities[entities["SMILES"].notna() & (entities["SMILES"] != "NA")][["entity_id", "SMILES"]]
print(f"Metabolites with usable SMILES: {len(smiles_df)}")

fps, ids = [], []
for _, row in smiles_df.iterrows():
    mol = Chem.MolFromSmiles(str(row["SMILES"]))
    if mol is not None:
        fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024))
        ids.append(row["entity_id"])
print(f"Fingerprinted: {len(fps)}")

THRESHOLD = 0.4
from_list, to_list, weight_list = [], [], []
for i in range(len(fps)):
    sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps)
    for j in range(i + 1, len(fps)):
        if sims[j] > THRESHOLD:
            from_list.append(ids[i]); to_list.append(ids[j]); weight_list.append(float(sims[j]))

edges = pd.DataFrame({"from": from_list, "to": to_list, "weight": weight_list})
edges.to_csv("data/processed/edges_structural.csv", index=False)
print(f"edges_structural.csv: {len(edges)} edges (Tanimoto > {THRESHOLD})")
