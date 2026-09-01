suppressPackageStartupMessages({library(readxl); library(dplyr); library(tidyr); library(stringr); library(purrr)})

data_path <- "data/raw/pd_with_adducts.xlsx"
main <- read_excel(data_path, sheet = "pd-with-adducts")
classyfire <- read_excel(data_path, sheet = "class of compounds")

studies <- c('S_24058461','S_24167579','S_25390405','S_28880465','S_29518718','S_33485385',
             'S_34577635','S_35829770','S_37335671','S_37589832','S_37755270','S_38516193',
             'S_38720721','ACS1517','MTBLS1219','MTBLS11904','ST001814','ST003159')

parse_cat_first <- function(x) map_chr(x, function(v) {
  if (is.na(v)) return(NA_character_)
  str_trim(str_split(as.character(v), "\\|")[[1]][1])
})

classyfire_dedup <- classyfire %>% select(InChIKey, Kingdom, Superclass, Class, Subclass) %>%
  distinct(InChIKey, .keep_all = TRUE)

entities <- main %>%
  select(entity_id, canonical_metabolite, SMILES, `KEGG ID`, InChIKey) %>%
  left_join(classyfire_dedup, by = "InChIKey")
write.csv(entities, "data/processed/entities_full.csv", row.names = FALSE)
cat("entities_full.csv:", nrow(entities), "rows,", sum(!is.na(entities$Class)), "with a known Class\n")


set.seed(123)
classified <- entities %>% filter(!is.na(Class))
holdout_ids <- classified %>% slice_sample(n = 200) %>% pull(entity_id)
write.csv(data.frame(entity_id = holdout_ids), "data/processed/holdout_ids.csv", row.names = FALSE)
cat("holdout_ids.csv:", length(holdout_ids), "metabolites\n")


build_detections <- function(s) {
  needed <- c('entity_id', paste0(s, c('__direction', '__p_value', '__fold_change',
                                        '__matched_adduct', '__matrix', '__mode')))
  if (!all(needed %in% colnames(main))) return(NULL)
  sub <- main[, needed]
  colnames(sub) <- c('entity_id', 'direction', 'p_value', 'fold_change', 'matched_adduct', 'matrix', 'mode')
  sub <- sub[!is.na(sub$direction) | !is.na(sub$p_value) | !is.na(sub$fold_change), ]
  if (nrow(sub) == 0) return(NULL)
  sub$study_id <- s
  sub$matched_adduct <- parse_cat_first(sub$matched_adduct)
  sub$matrix <- parse_cat_first(sub$matrix)
  sub$mode <- parse_cat_first(sub$mode)
  sub[, c('entity_id', 'study_id', 'matched_adduct', 'matrix', 'mode')]
}
detections <- bind_rows(map(studies, build_detections))
detections$platform <- paste(detections$matrix, detections$mode, sep = "_")
detections$platform[is.na(detections$matrix) | is.na(detections$mode)] <- NA

write.csv(detections %>% select(entity_id, study_id) %>% distinct(),
          "data/processed/entity_study_membership.csv", row.names = FALSE)
write.csv(detections %>% select(entity_id, matched_adduct) %>% filter(!is.na(matched_adduct)) %>% distinct(),
          "data/processed/entity_adduct_membership.csv", row.names = FALSE)
write.csv(detections %>% select(entity_id, platform) %>% filter(!is.na(platform)) %>% distinct(),
          "data/processed/entity_platform_membership.csv", row.names = FALSE)


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

write.csv(build_weighted_edges(detections, "entity_id", "study_id"), "data/processed/edges_study.csv", row.names = FALSE)
write.csv(build_weighted_edges(detections, "entity_id", "matched_adduct"), "data/processed/edges_adduct.csv", row.names = FALSE)
write.csv(build_weighted_edges(detections, "entity_id", "platform"), "data/processed/edges_platform.csv", row.names = FALSE)

