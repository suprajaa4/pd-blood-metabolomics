## ============================================================================
## Builds Metabolite -> Metabolite edges via shared KEGG Pathway membership,
## using the same inverse-group-size weighted projection pattern as every
## other layer (Class, Adduct, Platform, Structural). A shared pathway with
## few members counts more than a shared pathway with many - same principle
## throughout this whole project (the hub-domination fix from the Class
## enrichment work).
## ============================================================================
suppressPackageStartupMessages({library(readxl); library(dplyr); library(tidyr); library(stringr)})

data_path <- "data/raw/pd_with_adducts.xlsx"
pathway_path <- "data/raw/pathway_results.xlsx"

## ---- entity_id <-> KEGG ID is a direct column in the source sheet, no need
## to rebuild the whole PD-score pipeline just for this ----
main <- read_excel(data_path, sheet = "pd-with-adducts")
ent_kegg <- main %>%
  select(entity_id, canonical_metabolite, `KEGG ID`) %>%
  filter(!is.na(`KEGG ID`), `KEGG ID` != "NA") %>%
  distinct()
cat("Entities with a KEGG ID:", nrow(ent_kegg), "\n")

## ---- load KEGG pathway results and explode "one row per pathway" into
## "one row per (pathway, compound) pair" ----
pw <- read_excel(pathway_path, sheet = "pathway_results_kegg")
pw_long <- pw %>%
  select(Pathway, KEGG) %>%
  mutate(KEGG = str_split(KEGG, ",")) %>%
  unnest(KEGG) %>%
  mutate(KEGG = str_trim(KEGG)) %>%
  filter(KEGG != "")
cat("Pathway x compound pairs (KEGG-wide, before matching to our entities):", nrow(pw_long), "\n")
cat("Distinct pathways:", n_distinct(pw_long$Pathway), "\n")

## ---- join to OUR entities by KEGG ID (exact ID match - the whole reason
## KEGG was chosen over SMPDB for this step) ----
met_pathway <- pw_long %>%
  inner_join(ent_kegg, by = c("KEGG" = "KEGG ID"), relationship = "many-to-many")
cat("\nOf those, matched to our scored entities:", nrow(met_pathway), "rows,",
    n_distinct(met_pathway$entity_id), "distinct metabolites,",
    n_distinct(met_pathway$Pathway), "distinct pathways\n")

## ---- generic inverse-group-size weighted projection (same function used
## for every other layer in this project) ----
build_weighted_edges <- function(df, id_col, group_col) {
  d <- df %>% filter(!is.na(.data[[group_col]])) %>% distinct(.data[[id_col]], .data[[group_col]])
  colnames(d) <- c("met", "grp")
  grp_size <- d %>% count(grp, name = "group_size")
  d <- d %>% left_join(grp_size, by = "grp") %>% mutate(inv_weight = 1 / group_size)
  d %>%
    inner_join(d, by = "grp", suffix = c("_1", "_2"), relationship = "many-to-many") %>%
    filter(met_1 < met_2) %>%
    group_by(met_1, met_2) %>%
    summarise(weight = sum(inv_weight_1 * inv_weight_2), .groups = "drop") %>%
    rename(from = met_1, to = met_2)
}

edges_pathway <- build_weighted_edges(met_pathway, "entity_id", "Pathway")
write.csv(edges_pathway, "data/processed/edges_pathway.csv", row.names = FALSE)
cat("\nedges_pathway.csv:", nrow(edges_pathway), "edges among",
    n_distinct(c(edges_pathway$from, edges_pathway$to)), "metabolites\n")

## ---- also save the raw metabolite -> pathway membership table itself, in
## case you want the actual pathway NAMES later (e.g. for labeling plots),
## not just the derived metabolite-metabolite edges ----
write.csv(met_pathway %>% select(entity_id, canonical_metabolite, KEGG_ID = KEGG, Pathway),
          "data/processed/metabolite_pathway_membership.csv", row.names = FALSE)
cat("metabolite_pathway_membership.csv:", nrow(met_pathway), "rows\n")
