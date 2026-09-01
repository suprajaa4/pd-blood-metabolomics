suppressPackageStartupMessages({library(readxl); library(dplyr); library(tidyr); library(stringr); library(purrr)})

main <- read_excel("data/raw/pd_with_adducts.xlsx", sheet = "pd-with-adducts")
classyfire <- read_excel("data/raw/pd_with_adducts.xlsx", sheet = "class of compounds")
drug_sheet <- read_excel("data/raw/pd_with_adducts.xlsx", sheet = "drug")

studies <- c('S_24058461','S_24167579','S_25390405','S_28880465','S_29518718','S_33485385',
             'S_34577635','S_35829770','S_37335671','S_37589832','S_37755270','S_38516193',
             'S_38720721','ACS1517','MTBLS1219','MTBLS11904','ST001814','ST003159')
id_cols <- c('entity_group','entity_id','canonical_metabolite','all_metabolite_names',
             'n_studies','studies_present','InChIKey','formula','SMILES','HMDB ID','KEGG ID',
             'PubChem CID','annotation_score','has_formula','has_SMILES','has_HMDB','has_KEGG','has_PubChem')

parse_num_mean <- function(x) map_dbl(x, function(v) {
  if (is.na(v)) return(NA_real_)
  if (is.numeric(v)) return(as.numeric(v))
  nums <- suppressWarnings(as.numeric(str_trim(str_split(as.character(v), "\\|")[[1]])))
  nums <- nums[!is.na(nums)]
  if (length(nums) == 0) NA_real_ else mean(nums)
})
parse_num_min <- function(x) map_dbl(x, function(v) {
  if (is.na(v)) return(NA_real_)
  if (is.numeric(v)) return(as.numeric(v))
  nums <- suppressWarnings(as.numeric(str_trim(str_split(as.character(v), "\\|")[[1]])))
  nums <- nums[!is.na(nums)]
  if (length(nums) == 0) NA_real_ else min(nums)
})
parse_cat_first <- function(x) map_chr(x, function(v) {
  if (is.na(v)) return(NA_character_)
  str_trim(str_split(as.character(v), "\\|")[[1]][1])
})
parse_direction <- function(x) map_chr(x, function(v) {
  if (is.na(v)) return(NA_character_)
  parts <- str_trim(tolower(str_split(as.character(v), "\\|")[[1]]))
  norm <- case_when(str_detect(parts, "up") ~ "up", str_detect(parts, "down") ~ "down",
                     str_detect(parts, "no change") ~ "no_change", TRUE ~ parts)
  u <- unique(norm)
  if (setequal(u, "up")) return("up")
  if (setequal(u, "down")) return("down")
  if (setequal(u, "no_change")) return("no_change")
  if (setequal(u, c("up","no_change"))) return("up")
  if (setequal(u, c("down","no_change"))) return("down")
  if (all(c("up","down") %in% u)) return("conflicting_within_study")
  norm[1]
})

build_long_for_study <- function(s) {
  cols_needed <- c(id_cols, paste0(s, c('__repository','__phase','__mode','__matrix','__mz','__RT_min',
                    '__n_control','__n_PD','__p_value','__fold_change','__Log2FC','__direction',
                    '__matched_adduct','__control_female_ratio','__PD_female_ratio',
                    '__median_age_control','__median_age_PD')))
  if (!all(cols_needed %in% colnames(main))) return(NULL)
  sub <- main[, cols_needed]
  colnames(sub) <- c(id_cols, 'repository','phase','mode','matrix','mz','RT_min','n_control','n_PD',
                      'p_value','fold_change','Log2FC','direction','matched_adduct',
                      'control_female_ratio','PD_female_ratio','median_age_control','median_age_PD')
  keep <- !is.na(sub$direction) | !is.na(sub$p_value) | !is.na(sub$fold_change)
  sub <- sub[keep, ]
  if (nrow(sub) == 0) return(NULL)
  sub$study_id <- s
  sub$mode <- parse_cat_first(sub$mode); sub$matrix <- parse_cat_first(sub$matrix)
  sub$mz <- parse_num_mean(sub$mz); sub$RT_min <- parse_num_mean(sub$RT_min)
  sub$n_control <- parse_num_mean(sub$n_control); sub$n_PD <- parse_num_mean(sub$n_PD)
  sub$p_value <- parse_num_min(sub$p_value)
  sub$fold_change <- parse_num_mean(sub$fold_change); sub$Log2FC <- parse_num_mean(sub$Log2FC)
  sub$direction_raw <- sub$direction; sub$direction <- parse_direction(sub$direction)
  sub$matched_adduct <- parse_cat_first(sub$matched_adduct)
  sub$control_female_ratio <- parse_num_mean(sub$control_female_ratio)
  sub$PD_female_ratio <- parse_num_mean(sub$PD_female_ratio)
  sub$median_age_control <- parse_num_mean(sub$median_age_control)
  sub$median_age_PD <- parse_num_mean(sub$median_age_PD)
  sub
}
long_df <- bind_rows(map(studies, build_long_for_study))
need_fc <- is.na(long_df$Log2FC) & !is.na(long_df$fold_change) & long_df$fold_change > 0
long_df$Log2FC[need_fc] <- log2(long_df$fold_change[need_fc])
long_df <- long_df %>%
  mutate(age_delta_PD_minus_control = median_age_PD - median_age_control,
         female_ratio_delta_PD_minus_control = PD_female_ratio - control_female_ratio)
classyfire_dedup <- classyfire %>% select(InChIKey, Kingdom, Superclass, Class, Subclass) %>%
  distinct(InChIKey, .keep_all = TRUE)
long_df <- long_df %>% left_join(classyfire_dedup, by = "InChIKey")

long_df <- long_df %>%
  mutate(direction_final = case_when(
    !is.na(Log2FC) & Log2FC >  0.01 ~ "up",
    !is.na(Log2FC) & Log2FC < -0.01 ~ "down",
    !is.na(Log2FC) ~ "no_change",
    TRUE ~ direction)) %>%
  mutate(abs_log2fc = abs(Log2FC)) %>%
  group_by(study_id) %>%
  mutate(effect_percentile = percent_rank(abs_log2fc)) %>%
  ungroup() %>%
  mutate(neglog10p = -log10(pmax(p_value, 1e-300)))

minmax <- function(x) { rng <- range(x, na.rm=TRUE); if (diff(rng)==0) return(rep(0,length(x))); (x-rng[1])/diff(rng) }

ent <- long_df %>%
  group_by(entity_id) %>%
  summarise(
    canonical_metabolite = first(canonical_metabolite), InChIKey = first(InChIKey),
    formula = first(formula), SMILES = first(SMILES), `HMDB ID` = first(`HMDB ID`),
    `KEGG ID` = first(`KEGG ID`), `PubChem CID` = first(`PubChem CID`),
    Kingdom = first(Kingdom), Superclass = first(Superclass), Class = first(Class), Subclass = first(Subclass),
    annotation_score = first(annotation_score),
    has_formula = first(has_formula), has_SMILES = first(has_SMILES), has_HMDB = first(has_HMDB),
    has_KEGG = first(has_KEGG), has_PubChem = first(has_PubChem),
    n_studies_detected = n_distinct(study_id),
    n_dir_up = sum(direction_final == "up", na.rm = TRUE),
    n_dir_down = sum(direction_final == "down", na.rm = TRUE),
    mean_effect_percentile = mean(effect_percentile, na.rm = TRUE),
    mean_abs_log2fc = mean(abs_log2fc, na.rm = TRUE),
    mean_neglog10p = mean(neglog10p, na.rm = TRUE),
    n_studies_with_pvalue = sum(!is.na(p_value)),
    mean_age_delta_PD_minus_control = mean(age_delta_PD_minus_control, na.rm = TRUE),
    mean_female_ratio_delta_PD_minus_control = mean(female_ratio_delta_PD_minus_control, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(consistency_frac = ifelse((n_dir_up + n_dir_down) > 0, pmax(n_dir_up, n_dir_down) / (n_dir_up + n_dir_down), NA_real_),
         majority_direction = case_when(is.na(consistency_frac) ~ NA_character_, n_dir_up >= n_dir_down ~ "up", TRUE ~ "down"))

ent <- ent %>%
  mutate(comp_study_coverage = minmax(n_studies_detected),
         comp_consistency = ifelse(n_studies_detected <= 1, 0.5, ifelse(is.na(consistency_frac), 0.5, consistency_frac)),
         comp_effect_magnitude = ifelse(is.na(mean_effect_percentile), median(mean_effect_percentile, na.rm=TRUE), mean_effect_percentile),
         comp_identifier_completeness = (as.integer(has_formula)+as.integer(has_SMILES)+as.integer(has_HMDB)+as.integer(has_KEGG)+as.integer(has_PubChem))/5)
sig_norm <- minmax(ent$mean_neglog10p); fillval <- quantile(sig_norm, 0.25, na.rm=TRUE)
ent$comp_significance <- ifelse(is.na(sig_norm), fillval, sig_norm)
age_risk <- minmax(abs(ent$mean_age_delta_PD_minus_control))
ent$comp_age_confound_risk <- ifelse(is.na(age_risk), median(age_risk, na.rm=TRUE), age_risk)
gender_risk <- minmax(abs(ent$mean_female_ratio_delta_PD_minus_control))
ent$comp_gender_confound_risk <- ifelse(is.na(gender_risk), median(gender_risk, na.rm=TRUE), gender_risk)

## drug penalty
drug_lookup <- bind_rows(
  tibble(drug_name = colnames(drug_sheet)[1], item_type = "Parent drug / active compound"),
  drug_sheet %>% transmute(drug_name = .[[1]], item_type = item_type)
) %>% filter(!is.na(drug_name)) %>% mutate(drug_name_clean = str_trim(str_to_lower(drug_name)))
drug_severity <- c("Parent drug / active compound" = 1.00, "Drug, prodrug, ester, or salt/formulation" = 1.00,
                    "Veterinary / non-routine human-use drug" = 1.00, "Metabolite / degradation or conjugated form" = 0.50)
drug_lookup <- drug_lookup %>% mutate(severity = drug_severity[item_type]) %>%
  group_by(drug_name_clean) %>% summarise(severity = max(severity, na.rm=TRUE), .groups = "drop")
endogenous_dual_purpose <- str_to_lower(c("Hydrocortisone", "Cortisone", "Nicotinamide", "Pyridoxine"))
drug_lookup <- drug_lookup %>% mutate(severity = ifelse(drug_name_clean %in% endogenous_dual_purpose, 0.5, severity))
ent <- ent %>% mutate(.match_key = str_trim(str_to_lower(canonical_metabolite))) %>%
  left_join(drug_lookup, by = c(".match_key" = "drug_name_clean")) %>%
  mutate(is_known_drug = !is.na(severity), comp_drug_penalty = ifelse(is.na(severity), 0, severity)) %>%
  select(-.match_key, -severity)

W <- list(study_coverage=0.15, consistency=0.20, effect_magnitude=0.20, significance=0.20,
          identifier_completeness=0.10, age_confound_risk=0.075, gender_confound_risk=0.075, drug_penalty=0.15)
ent <- ent %>%
  mutate(PD_score_raw = W$study_coverage*comp_study_coverage + W$consistency*comp_consistency +
           W$effect_magnitude*comp_effect_magnitude + W$significance*comp_significance +
           W$identifier_completeness*comp_identifier_completeness -
           W$age_confound_risk*comp_age_confound_risk - W$gender_confound_risk*comp_gender_confound_risk -
           W$drug_penalty*comp_drug_penalty,
         PD_score = 100 * minmax(PD_score_raw),
         consistently_altered = (n_studies_detected >= 2) & (!is.na(consistency_frac) & consistency_frac >= 0.7)) %>%
  arrange(desc(PD_score))

export_df <- ent %>% select(entity_id, PD_score, consistently_altered, majority_direction,
                              n_studies_detected, is_known_drug,
                              has_formula, has_SMILES, has_HMDB, has_KEGG, has_PubChem, annotation_score)
write.csv(export_df, "data/processed/entity_pd_scores.csv", row.names = FALSE)
cat("Saved entity_pd_scores.csv:", nrow(export_df), "rows\n")
cat("Consistently altered:", sum(export_df$consistently_altered), "\n")
