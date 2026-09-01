#!/usr/bin/env python3
"""
Missing-modality-aware multimodal semi-supervised VAE for metabolite annotation.

Evidence modalities (matched as closely as possible to the graph model):
  1) Morgan fingerprint (structural evidence; optional)
  2) Study membership (optional)
  3) Adduct membership (optional)
  4) Platform membership (optional)
  5) Standardized m/z (optional)
  6) Non-target annotation:
       - Class task: full pathway-membership multi-hot (optional)
       - Pathway task: Class one-hot (optional)

Key design points:
  * Metabolites are NOT dropped when SMILES/fingerprints are missing.
  * Every modality has an observation mask.
  * Reconstruction loss is calculated only where a modality is observed.
  * Binary modalities use BCE-with-logits; m/z uses MSE.
  * Strict holdout entities are excluded from ALL training losses by default.
  * Unlabeled, non-holdout entities still contribute reconstruction + KL loss.
  * Optional modality dropout improves robustness to missing evidence.
  * Full per-class probabilities are saved for ROC-AUC/PR-AUC/top-k metrics.

Example:
    python scripts/13_multimodal_vae.py --data-dir data/processed --workbook data/raw/pd_with_adducts.xlsx --out-dir results/vae --epochs 60
"""

from __future__ import annotations
import argparse
import json
import os
import random
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
from sklearn.metrics import (
    accuracy_score, precision_recall_fscore_support, roc_auc_score,
    average_precision_score
)
from sklearn.preprocessing import label_binarize

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAVE_RDKIT = True
except Exception:
    HAVE_RDKIT = False

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ---------------- defaults / reproducibility ----------------
SEED = 42
LATENT_DIM = 24
BATCH_SIZE = 64
LR = 1e-3
WEIGHT_DECAY = 1e-4
DROPOUT = 0.30
MODALITY_DROPOUT = 0.10

BIN_RECON_WEIGHT = 1.0
MZ_RECON_WEIGHT = 1.0
KL_BETA = 0.20
CLASSIFICATION_WEIGHT = 5.0
KL_WARMUP_EPOCHS = 15


def seed_everything(seed: int = SEED):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def resolve_file(data_dir: Path, candidates, required=True):
    for name in candidates:
        p = data_dir / name
        if p.exists():
            return p
    # permissive stem matching for uploaded '(1)' style names
    normalized = [c.lower().replace(".csv", "").replace(".xlsx", "") for c in candidates]
    for p in data_dir.iterdir():
        low = p.name.lower()
        if any(n in low for n in normalized):
            return p
    if required:
        raise FileNotFoundError(f"Could not find any of: {candidates} in {data_dir}")
    return None


def clean_cat(s):
    s = s.astype(object)
    return s.where(pd.notna(s), None)


def build_multihot(membership_df, entity_ids, group_col, allowed_groups=None):
    """Return matrix, modality-presence mask, and group columns."""
    df = membership_df[["entity_id", group_col]].dropna().drop_duplicates().copy()
    if allowed_groups is None:
        groups = sorted(df[group_col].astype(str).unique())
    else:
        groups = list(allowed_groups)
    g2i = {g: i for i, g in enumerate(groups)}
    e2i = {e: i for i, e in enumerate(entity_ids)}
    X = np.zeros((len(entity_ids), len(groups)), dtype=np.float32)
    present = np.zeros((len(entity_ids), 1), dtype=np.float32)
    for eid, grp in zip(df["entity_id"], df[group_col].astype(str)):
        if eid in e2i and grp in g2i:
            X[e2i[eid], g2i[grp]] = 1.0
            present[e2i[eid], 0] = 1.0
    return X, present, groups


def build_onehot_from_training(values, training_mask):
    """Build one-hot vocabulary only from training-visible values."""
    vals = pd.Series(values).astype(object)
    train_vals = vals[training_mask & vals.notna()].astype(str)
    groups = sorted(train_vals.unique())
    g2i = {g: i for i, g in enumerate(groups)}
    X = np.zeros((len(vals), len(groups)), dtype=np.float32)
    present = np.zeros((len(vals), 1), dtype=np.float32)
    for i, v in enumerate(vals):
        if pd.isna(v) or str(v) not in g2i:
            continue
        X[i, g2i[str(v)]] = 1.0
        present[i, 0] = 1.0
    return X, present, groups


def make_fingerprints(entities, n_bits=1024):
    X = np.zeros((len(entities), n_bits), dtype=np.float32)
    present = np.zeros((len(entities), 1), dtype=np.float32)
    if not HAVE_RDKIT:
        print("WARNING: RDKit unavailable; fingerprint modality will be absent for all entities.")
        return X, present
    valid = 0
    for i, smi in enumerate(entities["SMILES"] if "SMILES" in entities else [None]*len(entities)):
        if pd.isna(smi) or str(smi).strip() in {"", "NA", "nan"}:
            continue
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
        X[i] = np.asarray(fp, dtype=np.float32)
        present[i, 0] = 1.0
        valid += 1
    print(f"Fingerprint coverage: {valid}/{len(entities)} ({valid/len(entities):.1%})")
    return X, present


def extract_mz_from_workbook(workbook, entity_ids):
    X = np.zeros((len(entity_ids), 1), dtype=np.float32)
    present = np.zeros((len(entity_ids), 1), dtype=np.float32)
    if workbook is None:
        print("WARNING: no workbook found; m/z modality absent.")
        return X, present, np.nan, np.nan

    try:
        main = pd.read_excel(workbook, sheet_name="pd-with-adducts")
    except Exception as e:
        print(f"WARNING: could not read m/z workbook ({e}); m/z modality absent.")
        return X, present, np.nan, np.nan

    if "entity_id" not in main.columns:
        return X, present, np.nan, np.nan

    mz_cols = [c for c in main.columns if str(c).endswith("__mz")]
    values = {}
    for _, row in main[["entity_id"] + mz_cols].iterrows():
        eid = row["entity_id"]
        arr = []
        for c in mz_cols:
            v = row[c]
            if pd.isna(v):
                continue
            try:
                arr.extend(float(x.strip()) for x in str(v).split("|") if x.strip())
            except Exception:
                continue
        if arr:
            values[eid] = float(np.median(arr))

    raw = pd.Series(entity_ids).map(values).astype(float)
    obs = raw.notna().to_numpy()
    if obs.sum() == 0:
        return X, present, np.nan, np.nan
    mean = float(raw[obs].mean())
    std = float(raw[obs].std(ddof=0))
    if not np.isfinite(std) or std == 0:
        std = 1.0
    X[obs, 0] = ((raw[obs] - mean) / std).astype(np.float32)
    present[obs, 0] = 1.0
    print(f"m/z coverage: {obs.sum()}/{len(entity_ids)} ({obs.mean():.1%})")
    return X, present, mean, std


class MultiModalDataset(Dataset):
    def __init__(self, modalities, masks, labels, train_allowed, labeled_train):
        self.modalities = {k: torch.tensor(v, dtype=torch.float32) for k, v in modalities.items()}
        self.masks = {k: torch.tensor(v, dtype=torch.float32) for k, v in masks.items()}
        self.labels = torch.tensor(labels, dtype=torch.long)
        self.train_allowed = torch.tensor(train_allowed.astype(np.float32), dtype=torch.float32)
        self.labeled_train = torch.tensor(labeled_train.astype(np.float32), dtype=torch.float32)
        self.n = len(labels)
    def __len__(self): return self.n
    def __getitem__(self, i):
        return (
            {k: v[i] for k, v in self.modalities.items()},
            {k: v[i] for k, v in self.masks.items()},
            self.labels[i], self.train_allowed[i], self.labeled_train[i]
        )


class MultiModalVAE(nn.Module):
    def __init__(self, modality_dims, n_classes, latent_dim=LATENT_DIM, dropout=DROPOUT):
        super().__init__()
        self.modality_names = list(modality_dims)
        # Modality-specific embedding sizes prevent 1024 fingerprint bits from numerically swamping small modalities.
        emb_dims = {}
        for name, d in modality_dims.items():
            if name == "fingerprint": emb_dims[name] = 128
            elif name in {"study", "adduct", "cross"}: emb_dims[name] = min(48, max(12, d))
            elif name == "platform": emb_dims[name] = min(24, max(8, d))
            elif name == "mz": emb_dims[name] = 8
            else: emb_dims[name] = min(32, max(8, d))
        self.emb_dims = emb_dims

        self.encoders = nn.ModuleDict({
            name: nn.Sequential(
                nn.Linear(d, emb_dims[name]), nn.ReLU(), nn.Dropout(dropout)
            ) for name, d in modality_dims.items()
        })
        fusion_dim = sum(emb_dims.values()) + len(modality_dims)  # append modality-presence flags
        self.fusion = nn.Sequential(
            nn.Linear(fusion_dim, 256), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(256, 128), nn.ReLU()
        )
        self.mu = nn.Linear(128, latent_dim)
        self.logvar = nn.Linear(128, latent_dim)

        self.decoders = nn.ModuleDict({
            name: nn.Sequential(
                nn.Linear(latent_dim, emb_dims[name]), nn.ReLU(),
                nn.Linear(emb_dims[name], d)
            ) for name, d in modality_dims.items()
        })
        self.classifier = nn.Sequential(
            nn.Linear(latent_dim, 64), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(64, n_classes)
        )

    def encode(self, xs, masks, apply_modality_dropout=False, modality_dropout=0.0):
        embs, effective_masks = [], {}
        # Sample-level modality dropout only on observed modalities.
        for name in self.modality_names:
            m = masks[name]
            if apply_modality_dropout and modality_dropout > 0:
                keep = (torch.rand_like(m) > modality_dropout).float()
                m_eff = m * keep
            else:
                m_eff = m
            h = self.encoders[name](xs[name]) * m_eff
            embs.append(h)
            effective_masks[name] = m_eff
        mask_vector = torch.cat([effective_masks[n] for n in self.modality_names], dim=1)
        fused = torch.cat(embs + [mask_vector], dim=1)
        h = self.fusion(fused)
        return self.mu(h), self.logvar(h), effective_masks

    @staticmethod
    def reparameterize(mu, logvar):
        return mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)

    def forward(self, xs, masks, modality_dropout=0.0):
        mu, logvar, eff_masks = self.encode(xs, masks, self.training, modality_dropout)
        z = self.reparameterize(mu, logvar) if self.training else mu
        recon = {name: self.decoders[name](z) for name in self.modality_names}
        logits = self.classifier(z)
        return recon, logits, mu, logvar, eff_masks


def masked_binary_loss(logits, target, modality_mask, sample_train_mask):
    # Mean per sample, then include only observed modality + train-allowed samples.
    per_sample = F.binary_cross_entropy_with_logits(logits, target, reduction="none").mean(dim=1)
    w = modality_mask.squeeze(1) * sample_train_mask
    denom = w.sum().clamp_min(1.0)
    return (per_sample * w).sum() / denom


def masked_mse(pred, target, modality_mask, sample_train_mask):
    per_sample = ((pred - target) ** 2).mean(dim=1)
    w = modality_mask.squeeze(1) * sample_train_mask
    denom = w.sum().clamp_min(1.0)
    return (per_sample * w).sum() / denom


def topk_accuracy(y_true_idx, probs, ks=(1,2,3,5)):
    out = {}
    order = np.argsort(-probs, axis=1)
    for k in ks:
        kk = min(k, probs.shape[1])
        out[f"top_{k}_accuracy"] = float(np.mean([y_true_idx[i] in order[i,:kk] for i in range(len(y_true_idx))]))
    return out


def multiclass_metrics(y_true_idx, probs, classes):
    pred_idx = probs.argmax(axis=1)
    p, r, f1, _ = precision_recall_fscore_support(
        y_true_idx, pred_idx, labels=np.arange(len(classes)), average="macro", zero_division=0
    )
    metrics = {
        "accuracy": float(accuracy_score(y_true_idx, pred_idx)),
        "macro_precision": float(p), "macro_recall": float(r), "macro_f1": float(f1),
    }
    metrics.update(topk_accuracy(y_true_idx, probs))
    try:
        y_bin = label_binarize(y_true_idx, classes=np.arange(len(classes)))
        metrics["macro_roc_auc_ovr"] = float(roc_auc_score(y_bin, probs, average="macro", multi_class="ovr"))
        metrics["macro_pr_auc"] = float(average_precision_score(y_bin, probs, average="macro"))
    except Exception:
        metrics["macro_roc_auc_ovr"] = np.nan
        metrics["macro_pr_auc"] = np.nan
    return metrics


def build_base_modalities(data_dir: Path, entities: pd.DataFrame, workbook=None):
    eids = entities["entity_id"].tolist()
    modalities, masks, meta = {}, {}, {}

    modalities["fingerprint"], masks["fingerprint"] = make_fingerprints(entities, 1024)

    specs = [
        ("study", ["entity_study_membership.csv"], "study_id"),
        ("adduct", ["entity_adduct_membership.csv"], "matched_adduct"),
        ("platform", ["entity_platform_membership.csv"], "platform"),
    ]
    for name, candidates, col in specs:
        p = resolve_file(data_dir, candidates, required=False)
        if p is None:
            print(f"WARNING: {name} membership missing; modality absent.")
            modalities[name] = np.zeros((len(eids), 1), dtype=np.float32)
            masks[name] = np.zeros((len(eids), 1), dtype=np.float32)
            meta[name+"_columns"] = []
        else:
            df = pd.read_csv(p)
            X, M, cols = build_multihot(df, eids, col)
            # Keep a 1-column placeholder only if vocabulary is genuinely empty.
            if X.shape[1] == 0:
                X = np.zeros((len(eids), 1), dtype=np.float32)
            modalities[name], masks[name] = X, M
            meta[name+"_columns"] = cols
            print(f"{name.title()} coverage: {int(M.sum())}/{len(eids)} ({M.mean():.1%}); dims={X.shape[1]}")

    workbook = Path(workbook) if workbook else resolve_file(data_dir, ["pd_with_adducts.xlsx", "pdwithadduct.xlsx", "pd_data.xlsx"], required=False)
    X_mz, M_mz, mz_mean, mz_std = extract_mz_from_workbook(workbook, eids)
    modalities["mz"], masks["mz"] = X_mz, M_mz
    meta.update({"mz_mean": mz_mean, "mz_std": mz_std, "workbook": str(workbook) if workbook else None})
    return modalities, masks, meta


def run_task(task, entities, primary, base_modalities, base_masks, data_dir, holdout_ids,
             epochs, out_dir, device, strict_holdout=True, modality_dropout=MODALITY_DROPOUT):
    df = entities.merge(primary[["entity_id", "primary_pathway"]], on="entity_id", how="left")
    eids = df["entity_id"].tolist()
    is_holdout = df["entity_id"].isin(holdout_ids).to_numpy()

    if task == "class":
        target_col = "Class"
        # Match graph pathway evidence more closely: use FULL known pathway membership as cross modality.
        p = resolve_file(data_dir, ["metabolite_pathway_membership.csv"], required=True)
        path_mem = pd.read_csv(p)
        # Cross vocabulary learned from non-holdout entities to avoid test-only categories defining input dimensions.
        train_eids = set(df.loc[~is_holdout, "entity_id"])
        allowed_paths = sorted(path_mem.loc[path_mem["entity_id"].isin(train_eids), "Pathway"].dropna().astype(str).unique())
        X_cross, M_cross, cross_cols = build_multihot(path_mem, eids, "Pathway", allowed_paths)
        cross_description = "full pathway-membership multi-hot"
    elif task == "pathway":
        target_col = "primary_pathway"
        train_visible = ~is_holdout
        X_cross, M_cross, cross_cols = build_onehot_from_training(df["Class"], train_visible)
        cross_description = "Class one-hot"
    else:
        raise ValueError(task)

    modalities = {k: v.copy() for k, v in base_modalities.items()}
    masks = {k: v.copy() for k, v in base_masks.items()}
    modalities["cross"], masks["cross"] = X_cross, M_cross

    target = df[target_col]
    known = target.notna().to_numpy()
    training_label_mask = known & ~is_holdout
    classes = sorted(target[training_label_mask].astype(str).unique())
    class_to_idx = {c:i for i,c in enumerate(classes)}
    labels = np.full(len(df), -1, dtype=np.int64)
    for i, v in enumerate(target):
        if pd.notna(v) and str(v) in class_to_idx:
            labels[i] = class_to_idx[str(v)]

    eval_mask = is_holdout & known & (labels >= 0)
    # Strict: holdout entities contribute no reconstruction/KL/classification training signal.
    train_allowed = (~is_holdout) if strict_holdout else np.ones(len(df), dtype=bool)
    labeled_train = training_label_mask & (labels >= 0)

    # Ensure every entity has at least one modality. Entities with zero observed modalities cannot be predicted meaningfully.
    any_modality = np.zeros(len(df), dtype=bool)
    for m in masks.values(): any_modality |= (m[:,0] > 0)

    print(f"\n=== {task.upper()} ===")
    print(f"Target classes in training: {len(classes)}")
    print(f"Labeled training entities: {labeled_train.sum()}")
    print(f"Holdout evaluable: {eval_mask.sum()}/{is_holdout.sum()}")
    print(f"Entities with >=1 observed modality: {any_modality.sum()}/{len(df)} ({any_modality.mean():.1%})")
    print(f"Cross modality: {cross_description}; dims={X_cross.shape[1]}, coverage={M_cross.mean():.1%}")

    ds = MultiModalDataset(modalities, masks, labels, train_allowed, labeled_train)
    loader = DataLoader(ds, batch_size=BATCH_SIZE, shuffle=True)
    modality_dims = {k: v.shape[1] for k,v in modalities.items()}
    model = MultiModalVAE(modality_dims, len(classes)).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)

    history = []
    binary_names = [n for n in modalities if n != "mz"]
    for epoch in range(1, epochs+1):
        model.train()
        sums = {"total":0., "binary":0., "mz":0., "kl":0., "clf":0.}; nb=0
        beta = KL_BETA * min(1.0, epoch / max(1, KL_WARMUP_EPOCHS))
        for xs, ms, y, train_ok, labeled_ok in loader:
            xs = {k:v.to(device) for k,v in xs.items()}
            ms = {k:v.to(device) for k,v in ms.items()}
            y, train_ok, labeled_ok = y.to(device), train_ok.to(device), labeled_ok.to(device)
            opt.zero_grad()
            recon, logits, mu, logvar, eff_masks = model(xs, ms, modality_dropout=modality_dropout)

            bin_losses = [masked_binary_loss(recon[n], xs[n], eff_masks[n], train_ok) for n in binary_names]
            # Equal weight per modality (not per individual bit) prevents fingerprint dimensionality dominating.
            binary_loss = torch.stack(bin_losses).mean() if bin_losses else torch.tensor(0., device=device)
            mz_loss = masked_mse(recon["mz"], xs["mz"], eff_masks["mz"], train_ok)

            kl_per = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1) / mu.shape[1]
            kl_denom = train_ok.sum().clamp_min(1.0)
            kl_loss = (kl_per * train_ok).sum() / kl_denom

            clf_mask = (labeled_ok > 0.5) & (y >= 0)
            if clf_mask.sum() > 0:
                clf_loss = F.cross_entropy(logits[clf_mask], y[clf_mask])
            else:
                clf_loss = torch.tensor(0., device=device)

            total = (BIN_RECON_WEIGHT * binary_loss + MZ_RECON_WEIGHT * mz_loss +
                     beta * kl_loss + CLASSIFICATION_WEIGHT * clf_loss)
            total.backward()
            opt.step()
            for k,v in [("total",total),("binary",binary_loss),("mz",mz_loss),("kl",kl_loss),("clf",clf_loss)]:
                sums[k] += float(v.detach().cpu())
            nb += 1
        row = {"epoch":epoch, "kl_beta":beta, **{k:v/max(nb,1) for k,v in sums.items()}}
        history.append(row)
        if epoch == 1 or epoch % 10 == 0 or epoch == epochs:
            print(f"epoch {epoch:>3}/{epochs} total={row['total']:.4f} bin={row['binary']:.4f} mz={row['mz']:.4f} KL={row['kl']:.4f} clf={row['clf']:.4f} beta={beta:.3f}")

    # deterministic evaluation from posterior means
    model.eval()
    all_probs = np.full((len(df), len(classes)), np.nan, dtype=np.float32)
    latent = np.full((len(df), LATENT_DIM), np.nan, dtype=np.float32)
    valid_idx = np.where(any_modality)[0]
    chunk = 512
    with torch.no_grad():
        for start in range(0, len(valid_idx), chunk):
            inds = valid_idx[start:start+chunk]
            xs = {k:torch.tensor(v[inds], dtype=torch.float32, device=device) for k,v in modalities.items()}
            ms = {k:torch.tensor(v[inds], dtype=torch.float32, device=device) for k,v in masks.items()}
            mu, _, _ = model.encode(xs, ms, apply_modality_dropout=False)
            probs = F.softmax(model.classifier(mu), dim=1)
            all_probs[inds] = probs.cpu().numpy()
            latent[inds] = mu.cpu().numpy()

    eval_idx = np.where(eval_mask & any_modality)[0]
    y_eval = labels[eval_idx]
    p_eval = all_probs[eval_idx]
    metrics = multiclass_metrics(y_eval, p_eval, classes)
    metrics.update({
        "task": task,
        "model": "multimodal_missing_aware_vae",
        "n_total_entities": int(len(df)),
        "n_entities_with_any_modality": int(any_modality.sum()),
        "dataset_modality_coverage": float(any_modality.mean()),
        "n_holdout_requested": int(is_holdout.sum()),
        "n_holdout_evaluable": int(len(eval_idx)),
        "holdout_coverage": float(len(eval_idx) / max(1, is_holdout.sum())),
        "n_training_classes": int(len(classes)),
        "n_labeled_training": int(labeled_train.sum()),
        "strict_holdout": bool(strict_holdout),
        "latent_dim": LATENT_DIM,
        "binary_recon_weight": BIN_RECON_WEIGHT,
        "mz_recon_weight": MZ_RECON_WEIGHT,
        "kl_beta_final": KL_BETA,
        "classification_weight": CLASSIFICATION_WEIGHT,
        "modality_dropout": modality_dropout,
    })

    # save holdout predictions and all class probabilities
    pred_idx = np.argmax(p_eval, axis=1)
    holdout_df = pd.DataFrame({
        "entity_id": df.loc[eval_idx,"entity_id"].values,
        "true_label": [classes[i] for i in y_eval],
        "predicted": [classes[i] for i in pred_idx],
        "confidence": p_eval.max(axis=1),
    })
    for j,c in enumerate(classes): holdout_df[f"prob__{c}"] = p_eval[:,j]
    holdout_df.to_csv(out_dir/f"multimodal_vae_{task}_results.csv", index=False)

    # real unknown predictions: target missing, any modality available
    unknown_idx = np.where((~known) & any_modality)[0]
    p_u = all_probs[unknown_idx]
    pred_u = np.argmax(p_u, axis=1) if len(unknown_idx) else np.array([], dtype=int)
    unknown_df = pd.DataFrame({
        "entity_id": df.loc[unknown_idx,"entity_id"].values,
        "canonical_metabolite": df.loc[unknown_idx,"canonical_metabolite"].values,
        f"predicted_{target_col}": [classes[i] for i in pred_u],
        "confidence": p_u.max(axis=1) if len(unknown_idx) else [],
    })
    # availability audit per prediction
    for n,m in masks.items(): unknown_df[f"has_{n}"] = m[unknown_idx,0].astype(int)
    unknown_df.to_csv(out_dir/f"multimodal_vae_predicted_{task}_unknowns.csv", index=False)

    latent_df = pd.DataFrame(latent, columns=[f"z{i+1}" for i in range(LATENT_DIM)])
    latent_df.insert(0,"entity_id",df["entity_id"].values)
    latent_df.to_csv(out_dir/f"multimodal_vae_{task}_latent_means.csv", index=False)
    pd.DataFrame(history).to_csv(out_dir/f"multimodal_vae_{task}_loss_history.csv", index=False)
    with open(out_dir/f"multimodal_vae_{task}_metrics.json","w") as f: json.dump(metrics,f,indent=2)

    print("Metrics:", {k:round(v,4) if isinstance(v,float) and np.isfinite(v) else v for k,v in metrics.items() if k in ["accuracy","macro_precision","macro_recall","macro_f1","top_3_accuracy","macro_roc_auc_ovr","macro_pr_auc","holdout_coverage"]})
    return metrics


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=".")
    ap.add_argument("--out-dir", default="results/vae")
    ap.add_argument("--workbook", default="data/raw/pd_with_adducts.xlsx")
    ap.add_argument("--epochs", type=int, default=60)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--modality-dropout", type=float, default=MODALITY_DROPOUT)
    ap.add_argument("--transductive-holdout", action="store_true", help="Allow holdout features to contribute reconstruction/KL (not recommended for strict comparison).")
    args = ap.parse_args()

    seed_everything(args.seed)
    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Device:", device)

    entities_p = resolve_file(data_dir,["entities_full.csv","entities_full(1).csv"])
    primary_p = resolve_file(data_dir,["entity_primary_pathway.csv"])
    class_hold_p = resolve_file(data_dir,["holdout_ids_v2_stratified.csv","holdout_ids_v2_stratified(1).csv"])
    path_hold_p = resolve_file(data_dir,["pathway_holdout_ids_v2_stratified.csv","pathway_holdout_ids_v2_stratified(1).csv"])

    entities = pd.read_csv(entities_p).drop_duplicates("entity_id").reset_index(drop=True)
    primary = pd.read_csv(primary_p).drop_duplicates("entity_id")
    class_holdout = set(pd.read_csv(class_hold_p)["entity_id"])
    path_holdout = set(pd.read_csv(path_hold_p)["entity_id"])

    base_modalities, base_masks, meta = build_base_modalities(data_dir, entities, args.workbook)
    print("\nBase modality coverage:")
    for n,m in base_masks.items(): print(f"  {n:12s}: {int(m.sum()):4d}/{len(m)} = {m.mean():.1%}")

    metrics = []
    metrics.append(run_task("class", entities, primary, base_modalities, base_masks, data_dir,
                            class_holdout, args.epochs, out_dir, device,
                            strict_holdout=not args.transductive_holdout,
                            modality_dropout=args.modality_dropout))
    metrics.append(run_task("pathway", entities, primary, base_modalities, base_masks, data_dir,
                            path_holdout, args.epochs, out_dir, device,
                            strict_holdout=not args.transductive_holdout,
                            modality_dropout=args.modality_dropout))
    pd.DataFrame(metrics).to_csv(out_dir/"multimodal_vae_metrics_summary.csv", index=False)
    with open(out_dir/"multimodal_vae_feature_metadata.json","w") as f:
        json.dump(meta, f, indent=2, default=str)
    print("\nSaved outputs to", out_dir)


if __name__ == "__main__":
    main()
