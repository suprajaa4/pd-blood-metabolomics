# Parkinson's Disease Blood Metabolomics Integration

A reproducible analysis repository for integrating Parkinson's disease blood metabolomics studies, harmonizing metabolite identifiers, evaluating cross-study consistency, computing a PD evidence score, constructing a multilayer metabolite graph, and comparing graph-based annotation with a multimodal VAE.

## Repository structure

```text
PD-metabolomics-github/
├── notebooks/                 # literature extraction and metabolite annotation
├── scripts/                   # numbered analysis pipeline
├── data/
│   ├── raw/                   # supplied source/curated data
│   └── processed/             # model-ready tables and graph inputs
├── results/
│   ├── exploratory/           # descriptive tables and plots
│   ├── graph/                 # graph predictions and graph figures
│   ├── vae/                   # multimodal VAE outputs
│   └── model_comparison/      # graph-vs-VAE evaluation figures
├── archive/                   # older experimental model variants
├── requirements.txt
├── R-packages.txt
└── .gitignore
```

## Pipeline

Run commands from the repository root.

| Step | File | Purpose |
|---|---|---|
| 0 | `notebooks/00_pubmed_extraction_pipeline.ipynb` | Extract metabolite information from literature-derived inputs. |
| 1 | `notebooks/01_metabolite_annotation.ipynb` | Add and harmonize metabolite identifiers/annotations. |
| 2 | `scripts/01_group_by_inchikey.py` | Group repeated metabolites by InChIKey and clean names. |
| 3 | `scripts/02_exploratory_analysis.R` | Cross-study exploratory analysis, recurrence, direction consistency, PCA/MDS and report figures. |
| 4 | `scripts/03_pd_score_umap.Rmd` | Compute PD evidence scores and UMAP-based analyses. |
| 5 | `scripts/04_adduct_matching.R` | Assign candidate adducts from formula-derived monoisotopic mass and observed m/z. |
| 6 | `scripts/05_rebuild_entities_and_memberships.R` | Rebuild entity tables, study/adduct/platform memberships and graph layers. |
| 7 | `scripts/06_recompute_pd_score.R` | Recompute the entity-level PD score used by graph analyses. |
| 8 | `scripts/07_build_class_edges.R` | Build ClassyFire-class graph edges. |
| 9 | `scripts/08_build_mz_edges.R` | Build m/z-neighborhood graph edges. |
| 10 | `scripts/09_build_pathway_edges.R` | Build KEGG pathway edges and pathway membership. |
| 11 | `scripts/10_build_structural_edges.py` | Build Morgan-fingerprint/Tanimoto structural edges. |
| 12 | `scripts/11_build_pathway_holdout.py` | Build the pathway label table and basic pathway holdout. |
| 13 | `scripts/12_graph_propagation.py` | Predict missing class/pathway annotations using the multilayer graph and generate downstream graph analyses. |
| 14 | `scripts/13_multimodal_vae.py` | Train/evaluate the missing-modality-aware multimodal VAE. |
| 15 | `scripts/16_plot_train_test_representativeness.py` | Check training/holdout label representativeness. |

`14_compute_metrics.py` and `15_plot_pr_roc.py` are retained because they were used for an earlier probability-file comparison workflow. The current multimodal VAE already writes its main evaluation metrics directly to `results/vae/`.

## Setup

Python:

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
```

R packages are listed in `R-packages.txt`. The R Markdown score workflow also uses `openxlsx`, `scales`, and packages loaded directly in that file.

## Main commands

Exploratory analysis:

```bash
Rscript scripts/02_exploratory_analysis.R
```

Rebuild graph inputs:

```bash
Rscript scripts/05_rebuild_entities_and_memberships.R
Rscript scripts/06_recompute_pd_score.R
Rscript scripts/07_build_class_edges.R
Rscript scripts/08_build_mz_edges.R
Rscript scripts/09_build_pathway_edges.R
python scripts/10_build_structural_edges.py
```

Graph annotation:

```bash
python scripts/12_graph_propagation.py
```

Multimodal VAE:

```bash
python scripts/13_multimodal_vae.py \
  --data-dir data/processed \
  --workbook data/raw/pd_with_adducts.xlsx \
  --out-dir results/vae \
  --epochs 60
```

## Important reproducibility notes

The supplied ZIP did **not** contain the original ungrouped `pd_met_final(1).csv` expected by the grouping script. Therefore `scripts/01_group_by_inchikey.py` is preserved and cleaned, but requires that source file to be placed at `data/raw/pd_met_final.csv` if the full pipeline is rerun from the pre-grouping stage. The supplied grouped and direction-corrected dataset is already available as `data/raw/metabolites_grouped_direction_corrected.csv`, so all downstream exploratory work can start there.

The stratified class/pathway holdout files and corresponding training tables are retained in `data/processed/` because the ZIP did not contain the script that originally created the `*_v2_stratified.csv` split. These are treated as curated model inputs rather than silently regenerated with different sampling logic.

`edges_study.csv` and `edges_platform.csv` are not committed because the generated files are approximately 339 MB and 188 MB, above GitHub's standard 100 MB single-file limit. They are reproducible from `scripts/05_rebuild_entities_and_memberships.R` and are ignored by `.gitignore`.

## PD score

The project combines recurrence across studies, directional consistency, effect magnitude, statistical significance, identifier completeness, and demographic confounding penalties into an entity-level PD evidence score. The complete implementation and component outputs are in `scripts/03_pd_score_umap.Rmd` and `scripts/06_recompute_pd_score.R`.

## Archive

`archive/experimental_models/` contains older VAE, tiered-prediction and graph-propagation variants from development. They are kept for provenance but are not part of the canonical run order above.
