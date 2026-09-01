#exploratory plots

library(tidyverse)
library(dplyr)

input_file <- "data/raw/metabolites_grouped_direction_corrected.csv"
out_dir    <- "results/exploratory"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots_png"), showWarnings = FALSE, recursive = TRUE)

## ---- 1. Helper functions ----------------------------------------------------
read_final_met <- function(path) {
  ## FINAL_MET.csv contains special characters from copied metabolite names.
  ## Windows-1252 is the safest default from this file inspection.
  readr::read_csv(
    file = path,
    locale = readr::locale(encoding = "Windows-1252"),
    show_col_types = FALSE,
    guess_max = 100000
  )
}

safe_num <- function(x) {
  x <- as.character(x)
  x <- ifelse(x %in% c("", "NA", "NaN", "NULL", "null"), NA_character_, x)
  suppressWarnings(as.numeric(x))
}

paste_unique <- function(x, sep = " | ") {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  paste(sort(unique(x)), collapse = sep)
}

mode_value <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

clean_name <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", " ") %>%
    str_squish()
}

save_plot <- function(plot_obj, filename, width = 11, height = 7) {
  ggsave(file.path(out_dir, "plots_png", paste0(filename, ".png")), plot_obj,
         width = width, height = height, dpi = 300)
  plot_obj
}

blank_plot <- function(title, subtitle = "No usable data available for this plot") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = subtitle, size = 5) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

plot_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F1F3F5", colour = "#D0D7DE")
  )

## Simple, non-complex colour palette used across the report.
## Teal/blue/orange/red/green are intentionally reused so the plots stay colourful
## without becoming visually noisy.
col_teal   <- "#2A9D8F"
col_blue   <- "#457B9D"
col_orange <- "#F4A261"
col_red    <- "#E76F51"
col_green  <- "#6A994E"
col_purple <- "#8E7CC3"
col_yellow <- "#E9C46A"
col_grey   <- "#8D99AE"
col_light  <- "#F7F7F7"

basic_palette <- c(col_teal, col_blue, col_orange, col_red, col_green, col_purple, col_yellow, col_grey)
direction_palette <- c("Up in PD" = col_red, "Down in PD" = col_blue, "No change" = col_grey, "Unknown" = "#B0B0B0")
present_palette   <- c("Present" = col_green, "Missing" = "#F2E8CF")
discord_palette   <- c("FALSE" = col_green, "TRUE" = col_red, "Consistent" = col_green, "Discordant" = col_red)
mode_palette      <- c("Positive" = col_red, "Negative" = col_blue, "Mixed" = col_purple, "Dual" = col_purple, "Not Reported" = col_grey, "Not reported" = col_grey)
matrix_palette    <- c("serum" = col_teal, "plasma" = col_orange, "csf" = col_purple, "Not reported" = col_grey, "not reported" = col_grey)

## ---- 2. Read wide dataset ---------------------------------------------------
wide <- read_final_met(input_file)

base_cols <- c(
  "entity_group", "entity_id", "canonical_metabolite", "all_metabolite_names",
  "n_input_rows", "n_studies", "studies_present", "InChIKey", "formula",
  "SMILES", "HMDB ID", "KEGG ID", "PubChem CID", "annotation_score",
  "has_formula", "has_SMILES", "has_HMDB", "has_KEGG", "has_PubChem"
)
base_cols <- intersect(base_cols, names(wide))

study_cols <- names(wide)[str_detect(names(wide), "__")]
study_prefixes <- sort(unique(str_replace(study_cols, "__.*$", "")))
study_ids <- str_remove(study_prefixes, "^S_")

message("Rows/entities: ", nrow(wide))
message("Study prefixes detected: ", length(study_prefixes), " -> ", paste(study_ids, collapse = ", "))

## ---- 3. Convert wide study columns to long format ---------------------------
## This is the key step. It preserves study IDs with letters, e.g. MTBLS11904,
## ST003159, ACS1517. Numeric PubMed IDs originally stored as S_12345678 are
## displayed without the S_ prefix in study_id.
## IMPORTANT: study-specific columns include both character fields, such as
## repository/mode/phase, and numeric fields, such as n_control/Log2FC.
## pivot_longer() must coerce all values to character first; otherwise tidyr
## tries to combine character and double columns and stops with a type error.
long <- wide %>%
  select(any_of(base_cols), all_of(study_cols)) %>%
  pivot_longer(
    cols = all_of(study_cols),
    names_to = c("study_prefix", "measure"),
    names_sep = "__",
    values_to = "value",
    values_drop_na = FALSE,
    values_transform = list(value = as.character)
  ) %>%
  mutate(
    study_id = str_remove(study_prefix, "^S_"),
    value = na_if(value, "")
  ) %>%
  group_by(across(any_of(base_cols)), study_prefix, study_id, measure) %>%
  summarise(value = paste_unique(value), .groups = "drop") %>%
  pivot_wider(names_from = measure, values_from = value)

## Drop blank study rows where no study-specific values are present.
study_measure_cols <- setdiff(names(long), c(base_cols, "study_prefix", "study_id"))
long <- long %>%
  filter(if_any(all_of(study_measure_cols), ~ !is.na(.x) & .x != ""))

## Numeric conversion for analysis.
numeric_cols <- intersect(c(
  "mz", "RT_min", "n_control", "int_mean_control", "int_median_control",
  "log_control_mean", "log_control_median", "z_control", "n_PD",
  "int_mean_PD", "int_median_PD", "log_PD_mean", "log_PD_median",
  "z_PD", "z_PD_minus_control", "p_value", "fold_change", "Log2FC",
  "annotation_score", "n_input_rows", "n_studies"
), names(long))

for (cc in numeric_cols) long[[cc]] <- safe_num(long[[cc]])

long <- long %>%
  mutate(
    entity_label = coalesce(canonical_metabolite, entity_id),
    entity_label_short = str_trunc(entity_label, 45),
    entity_has_inchikey = !is.na(InChIKey) & InChIKey != "",
    direction_from_log2fc = case_when(
      is.na(Log2FC) ~ NA_character_,
      Log2FC > 0 ~ "Up in PD",
      Log2FC < 0 ~ "Down in PD",
      TRUE ~ "No change"
    ),
    direction_final = coalesce(direction, direction_from_log2fc),
    neg_log10_p = if_else(!is.na(p_value) & p_value > 0, -log10(p_value), NA_real_),
    abs_log2fc = abs(Log2FC),
    raw_ratio_PD_control = if_else(!is.na(int_mean_PD) & !is.na(int_mean_control) & int_mean_control != 0,
                                   int_mean_PD / int_mean_control, NA_real_),
    mode_clean = case_when(
      is.na(mode) | mode == "" ~ "Not reported",
      TRUE ~ str_to_title(str_squish(mode))
    ),
    phase_clean = case_when(
      is.na(phase) | phase == "" ~ "Not reported",
      str_detect(phase, regex("hilic", ignore_case = TRUE)) ~ "HILIC",
      str_detect(phase, regex("c18|rp|reverse", ignore_case = TRUE)) ~ "Reverse phase/C18",
      str_detect(phase, regex("gc", ignore_case = TRUE)) ~ "GC-MS",
      str_detect(phase, regex("lc", ignore_case = TRUE)) ~ "LC-MS/unspecified",
      TRUE ~ phase
    ),
    matrix_clean = case_when(
      is.na(matrix) | matrix == "" ~ "Not reported",
      TRUE ~ str_to_lower(str_squish(matrix))
    )
  )

## ---- 4. Entity-level consistency tables ------------------------------------
fc_long <- long %>%
  filter(!is.na(Log2FC), direction_from_log2fc %in% c("Up in PD", "Down in PD"))

entity_consistency <- fc_long %>%
  group_by(entity_id) %>%
  summarise(
    canonical_metabolite = mode_value(canonical_metabolite),
    all_metabolite_names = mode_value(all_metabolite_names),
    InChIKey = mode_value(InChIKey),
    formula = mode_value(formula),
    n_studies_with_log2fc = n_distinct(study_id),
    n_effect_observations = n(),
    n_up = sum(Log2FC > 0, na.rm = TRUE),
    n_down = sum(Log2FC < 0, na.rm = TRUE),
    mean_log2fc = mean(Log2FC, na.rm = TRUE),
    median_log2fc = median(Log2FC, na.rm = TRUE),
    max_abs_log2fc = max(abs(Log2FC), na.rm = TRUE),
    min_p_value = suppressWarnings(min(p_value, na.rm = TRUE)),
    studies = paste_unique(study_id),
    .groups = "drop"
  ) %>%
  mutate(
    min_p_value = if_else(is.infinite(min_p_value), NA_real_, min_p_value),
    dominant_direction = case_when(
      n_up > n_down ~ "Up in PD",
      n_down > n_up ~ "Down in PD",
      TRUE ~ "Discordant/tie"
    ),
    consistency_score = pmax(n_up, n_down) / pmax(n_up + n_down, 1),
    discordant = n_up > 0 & n_down > 0,
    recurrent = n_studies_with_log2fc >= 2
  ) %>%
  arrange(desc(n_studies_with_log2fc), desc(consistency_score), desc(max_abs_log2fc))

study_summary <- long %>%
  group_by(study_id) %>%
  summarise(
    repository = mode_value(repository),
    phase = mode_value(phase_clean),
    mode = mode_value(mode_clean),
    matrix = mode_value(matrix_clean),
    n_entities_with_any_data = n_distinct(entity_id),
    n_log2fc = sum(!is.na(Log2FC)),
    n_up = sum(Log2FC > 0, na.rm = TRUE),
    n_down = sum(Log2FC < 0, na.rm = TRUE),
    n_p_lt_0_05 = sum(!is.na(p_value) & p_value < 0.05),
    median_abs_log2fc = median(abs(Log2FC), na.rm = TRUE),
    median_mz = median(mz, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_log2fc))

annotation_summary <- wide %>%
  summarise(
    rows_entities = n(),
    InChIKey = mean(!is.na(InChIKey) & InChIKey != ""),
    formula = mean(!is.na(formula) & formula != ""),
    SMILES = mean(!is.na(SMILES) & SMILES != ""),
    HMDB = mean(!is.na(`HMDB ID`) & `HMDB ID` != ""),
    KEGG = mean(!is.na(`KEGG ID`) & `KEGG ID` != ""),
    PubChem = mean(!is.na(`PubChem CID`) & `PubChem CID` != "")
  )

write_csv(long, file.path(out_dir, "long_study_metabolite_table.csv"))
write_csv(entity_consistency, file.path(out_dir, "inchikey_entity_up_down_consistency.csv"))
write_csv(study_summary, file.path(out_dir, "study_summary.csv"))
write_csv(annotation_summary, file.path(out_dir, "annotation_summary.csv"))
write_csv(entity_consistency %>% filter(recurrent, discordant), file.path(out_dir, "discordant_recurrent_entities.csv"))
write_csv(entity_consistency %>% filter(recurrent, !discordant), file.path(out_dir, "consistent_recurrent_entities.csv"))

## ---- 5. Plot data prep ------------------------------------------------------
top_studies <- study_summary %>% slice_max(n_log2fc, n = 12, with_ties = FALSE) %>% pull(study_id)

## For heatmaps, keep recurrent entities observed in >=2 studies and rank by recurrence and effect size.
heatmap_entities <- entity_consistency %>%
  filter(recurrent) %>%
  arrange(desc(n_studies_with_log2fc), desc(max_abs_log2fc)) %>%
  slice_head(n = 40) %>%
  pull(entity_id)

heatmap_df <- fc_long %>%
  filter(entity_id %in% heatmap_entities) %>%
  group_by(entity_id, study_id) %>%
  summarise(
    canonical_metabolite = mode_value(canonical_metabolite),
    Log2FC = mean(Log2FC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = str_trunc(canonical_metabolite, 38))

## ---- 6. Plots ---------------------------------------------------------------
plots <- list()

plots$entity_group_counts <- wide %>%
  count(entity_group, name = "n") %>%
  ggplot(aes(x = reorder(entity_group, n), y = n, fill = entity_group)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  labs(title = "Inchi", x = NULL, y = "Number of metabolite rows/entities") +
  plot_theme

plots$studies_per_entity <- wide %>%
  mutate(n_studies = safe_num(n_studies)) %>%
  ggplot(aes(x = n_studies)) +
  geom_histogram(binwidth = 1, boundary = 0.5, colour = "white", fill = col_blue) +
  scale_x_continuous(breaks = sort(unique(safe_num(wide$n_studies)))) +
  labs(title = "Number of studies represented per metabolite entity", x = "Studies per entity", y = "Count") +
  plot_theme

plots$metabolites_per_study <- study_summary %>%
  mutate(study_id = reorder(study_id, n_entities_with_any_data)) %>%
  ggplot(aes(x = study_id, y = n_entities_with_any_data)) +
  geom_col(fill = col_teal) +
  coord_flip() +
  labs(title = "Metabolite/entity coverage by study", x = "Study", y = "Entities with any study data") +
  plot_theme

plots$direction_by_study <- fc_long %>%
  count(study_id, direction_from_log2fc, name = "n") %>%
  ggplot(aes(x = reorder(study_id, n), y = n, fill = direction_from_log2fc)) +
  geom_col() +
  coord_flip() +
  labs(title = "Up/down regulation counts by study", x = "Study", y = "Number of Log2FC observations", fill = "Direction") +
  scale_fill_manual(values = direction_palette, na.value = col_grey) +
  plot_theme

plots$annotation_coverage <- annotation_summary %>%
  pivot_longer(cols = -rows_entities, names_to = "annotation", values_to = "coverage") %>%
  mutate(annotation = reorder(annotation, coverage)) %>%
  ggplot(aes(x = annotation, y = coverage)) +
  geom_col(fill = col_teal) +
  geom_text(aes(label = paste0(round(coverage * 100, 1), "%")), hjust = -0.1, size = 3.3) +
  coord_flip() +
  ylim(0, 1.05) +
  labs(title = "Structural annotation coverage", x = NULL, y = "Fraction of metabolite entities") +
  plot_theme

plots$annotation_score_distribution <- wide %>%
  mutate(annotation_score = safe_num(annotation_score)) %>%
  filter(!is.na(annotation_score)) %>%
  ggplot(aes(x = annotation_score)) +
  geom_histogram(binwidth = 1, colour = "white", fill = col_blue) +
  labs(title = "Annotation score distribution", x = "Annotation score", y = "Count") +
  plot_theme

plots$mode_by_study <- long %>%
  distinct(study_id, entity_id, mode_clean) %>%
  count(study_id, mode_clean, name = "n") %>%
  ggplot(aes(x = reorder(study_id, n), y = n, fill = mode_clean)) +
  geom_col() +
  coord_flip() +
  labs(title = "Ion mode composition by study", x = "Study", y = "Entity count", fill = "Mode") +
  scale_fill_brewer(palette = "Set2", na.value = col_grey) +
  plot_theme

plots$phase_by_study <- long %>%
  distinct(study_id, entity_id, phase_clean) %>%
  count(study_id, phase_clean, name = "n") %>%
  ggplot(aes(x = reorder(study_id, n), y = n, fill = phase_clean)) +
  geom_col() +
  coord_flip() +
  labs(title = "Analytical phase/platform category by study", x = "Study", y = "Entity count", fill = "Phase") +
  scale_fill_brewer(palette = "Set3", na.value = col_grey) +
  plot_theme

plots$matrix_by_study <- long %>%
  distinct(study_id, entity_id, matrix_clean) %>%
  count(study_id, matrix_clean, name = "n") %>%
  ggplot(aes(x = reorder(study_id, n), y = n, fill = matrix_clean)) +
  geom_col() +
  coord_flip() +
  labs(title = "Sample matrix by study", x = "Study", y = "Entity count", fill = "Matrix") +
  scale_fill_manual(values = matrix_palette, na.value = col_grey) +
  plot_theme

plots$mz_density_by_mode <- long %>%
  filter(!is.na(mz), mode_clean %in% c("Positive", "Negative")) %>%
  ggplot(aes(x = mz, colour = mode_clean, fill = mode_clean)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  labs(title = "m/z distribution by ion mode", x = "m/z", y = "Density", colour = "Mode", fill = "Mode") +
  scale_fill_manual(values = mode_palette, na.value = col_grey) +
  scale_colour_manual(values = mode_palette, na.value = col_grey) +
  plot_theme

plots$mz_vs_rt <- long %>%
  filter(!is.na(mz), !is.na(RT_min)) %>%
  ggplot(aes(x = RT_min, y = mz, colour = mode_clean)) +
  geom_point(alpha = 0.35, size = 1) +
  labs(title = "m/z versus retention time", x = "Retention time (min)", y = "m/z", colour = "Mode") +
  scale_colour_manual(values = mode_palette, na.value = col_grey) +
  plot_theme

plots$log2fc_density <- fc_long %>%
  filter(study_id %in% top_studies) %>%
  ggplot(aes(x = Log2FC)) +
  geom_density(fill = col_yellow, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ study_id, scales = "free_y") +
  labs(title = "Log2 fold-change distributions in the largest studies", x = "Log2FC (PD versus control)", y = "Density") +
  plot_theme + theme(axis.text.x = element_text(angle = 0))

plots$volcano_all <- long %>%
  filter(!is.na(Log2FC), !is.na(neg_log10_p)) %>%
  ggplot(aes(x = Log2FC, y = neg_log10_p, colour = direction_from_log2fc)) +
  geom_point(alpha = 0.45, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  labs(title = "Volcano-style overview across all study-metabolite observations",
       x = "Log2FC", y = "-log10(p-value)", colour = "Direction") +
  scale_colour_manual(values = direction_palette, na.value = col_grey) +
  plot_theme + theme(axis.text.x = element_text(angle = 0))

plots$log2fc_heatmap <- heatmap_df %>%
  ggplot(aes(x = study_id, y = reorder(label, Log2FC), fill = Log2FC)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  labs(title = "Recurrent metabolite Log2FC heatmap by study",
       x = "Study", y = "Metabolite/entity", fill = "Log2FC") +
  scale_fill_gradient2(low = col_blue, mid = "white", high = col_red, midpoint = 0, na.value = "grey95") +
  plot_theme + theme(axis.text.y = element_text(size = 7))

plots$consistency_score_distribution <- entity_consistency %>%
  filter(recurrent) %>%
  ggplot(aes(x = consistency_score)) +
  geom_histogram(binwidth = 0.1, boundary = 0, colour = "white", fill = col_blue) +
  labs(title = "Up/down consistency scores for recurrent entities", x = "Consistency score", y = "Count") +
  plot_theme + theme(axis.text.x = element_text(angle = 0))

plots$recurrent_direction_counts <- entity_consistency %>%
  filter(recurrent) %>%
  count(dominant_direction, discordant, name = "n") %>%
  ggplot(aes(x = dominant_direction, y = n, fill = discordant)) +
  geom_col() +
  labs(title = "Direction calls among recurrent entities", x = "Dominant direction", y = "Entity count", fill = "Discordant?") +
  scale_fill_manual(values = discord_palette, na.value = col_grey) +
  plot_theme + theme(axis.text.x = element_text(angle = 0))

plots$top_consistent_up <- entity_consistency %>%
  filter(recurrent, !discordant, dominant_direction == "Up in PD") %>%
  slice_max(mean_log2fc, n = 20, with_ties = FALSE) %>%
  mutate(label = reorder(str_trunc(canonical_metabolite, 42), mean_log2fc)) %>%
  ggplot(aes(x = label, y = mean_log2fc)) +
  geom_col(fill = col_red) +
  coord_flip() +
  labs(title = "Most consistently up-regulated recurrent entities", x = NULL, y = "Mean Log2FC") +
  plot_theme

plots$top_consistent_down <- entity_consistency %>%
  filter(recurrent, !discordant, dominant_direction == "Down in PD") %>%
  slice_min(mean_log2fc, n = 20, with_ties = FALSE) %>%
  mutate(label = reorder(str_trunc(canonical_metabolite, 42), -mean_log2fc)) %>%
  ggplot(aes(x = label, y = mean_log2fc)) +
  geom_col(fill = col_blue) +
  coord_flip() +
  labs(title = "Most consistently down-regulated recurrent entities", x = NULL, y = "Mean Log2FC") +
  plot_theme

plots$discordant_entities <- entity_consistency %>%
  filter(recurrent, discordant) %>%
  arrange(desc(n_studies_with_log2fc), consistency_score) %>%
  slice_head(n = 25) %>%
  mutate(label = reorder(str_trunc(canonical_metabolite, 42), n_studies_with_log2fc)) %>%
  ggplot(aes(x = label, y = n_studies_with_log2fc, fill = consistency_score)) +
  geom_col() +
  coord_flip() +
  labs(title = "Discordant recurrent entities", x = NULL, y = "Studies with Log2FC", fill = "Consistency") +
  scale_fill_gradient(low = col_yellow, high = col_red, na.value = "grey90") +
  plot_theme

plots$z_difference_vs_log2fc <- long %>%
  filter(!is.na(Log2FC), !is.na(z_PD_minus_control)) %>%
  ggplot(aes(x = Log2FC, y = z_PD_minus_control, colour = study_id)) +
  geom_point(alpha = 0.45, size = 1.1, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Study-scaled Z difference versus Log2FC",
       x = "Log2FC", y = "Z-scaled PD - control") +
  plot_theme + theme(axis.text.x = element_text(angle = 0))

plots$raw_intensity_pd_vs_control <- long %>%
  filter(!is.na(int_mean_control), !is.na(int_mean_PD), int_mean_control > 0, int_mean_PD > 0) %>%
  ggplot(aes(x = int_mean_control, y = int_mean_PD, colour = study_id)) +
  geom_point(alpha = 0.35, size = 1, show.legend = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Raw/  group mean intensities: PD versus control",
        x = "Control mean intensity", y = "PD mean intensity") +
  plot_theme + theme(axis.text.x = element_text(angle = 0))

plots$raw_ratio_vs_log2fc <- long %>%
  filter(!is.na(raw_ratio_PD_control), raw_ratio_PD_control > 0, !is.na(Log2FC)) %>%
  ggplot(aes(x = log2(raw_ratio_PD_control), y = Log2FC, colour = study_id)) +
  geom_point(alpha = 0.35, size = 1, show.legend = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "raw-intensity ratio versus Log2FC",
       x = "log2(PD mean / control mean)", y = "Log2FC") +
  plot_theme + theme(axis.text.x = element_text(angle = 0))

plots$pvalue_distribution <- long %>%
  filter(!is.na(p_value), p_value >= 0, p_value <= 1, study_id %in% top_studies) %>%
  ggplot(aes(x = p_value)) +
  geom_histogram(bins = 30, colour = "white", fill = col_blue) +
  facet_wrap(~ study_id, scales = "free_y") +
  labs(title = "P-value distributions in largest studies", x = "p-value", y = "Count") +
  plot_theme + theme(axis.text.x = element_text(angle = 0))

plots$annotation_score_vs_abs_effect <- long %>%
  filter(!is.na(annotation_score), !is.na(abs_log2fc)) %>%
  ggplot(aes(x = factor(annotation_score), y = abs_log2fc)) +
  geom_boxplot(outlier.alpha = 0.15, fill = col_yellow, colour = col_blue) +
  labs(title = "Annotation score versus absolute effect size", x = "Annotation score", y = "|Log2FC|") +
  plot_theme + theme(axis.text.x = element_text(angle = 0))



## ---- 6B. Additional requested overview, recurrence, discordance, and ordination plots ----

## A. What information is available in each study?
## This creates a simple study-by-information heatmap similar to the study summary table.
study_data_availability <- long %>%
  group_by(study_id) %>%
  summarise(
    `Metabolite/entity` = n_distinct(entity_id) > 0,
    `Sample sizes` = any(!is.na(n_control) | !is.na(n_PD)),
    `Raw intensity` = any(!is.na(int_mean_control) | !is.na(int_mean_PD) |
                            !is.na(int_median_control) | !is.na(int_median_PD)),
    `Log intensity` = any(!is.na(log_control_mean) | !is.na(log_PD_mean) |
                            !is.na(log_control_median) | !is.na(log_PD_median)),
    `Z-scaled intensity` = any(!is.na(z_control) | !is.na(z_PD) | !is.na(z_PD_minus_control)),
    `P-value` = any(!is.na(p_value)),
    `Fold-change` = any(!is.na(fold_change)),
    `Curated Log2FC` = any(!is.na(Log2FC)),
    `m/z` = any(!is.na(mz)),
    `Retention time` = any(!is.na(RT_min)),
    `Ion mode` = any(!is.na(mode_clean) & mode_clean != "Not reported"),
    `Phase/platform` = any(!is.na(phase_clean) & phase_clean != "Not reported"),
    `Sample matrix` = any(!is.na(matrix_clean) & matrix_clean != "not reported"),
    `Annotation IDs` = any((!is.na(InChIKey) & InChIKey != "") |
                           (!is.na(formula) & formula != "") |
                           (!is.na(`HMDB ID`) & `HMDB ID` != "") |
                           (!is.na(`KEGG ID`) & `KEGG ID` != "") |
                           (!is.na(`PubChem CID`) & `PubChem CID` != "")),
    .groups = "drop"
  )

study_data_availability_long <- study_data_availability %>%
  pivot_longer(cols = -study_id, names_to = "info_type", values_to = "present") %>%
  mutate(
    present_label = if_else(present, "Present", "Missing"),
    info_type = factor(info_type, levels = rev(c(
      "Metabolite/entity", "Annotation IDs", "Sample sizes", "Sample matrix",
      "Ion mode", "Phase/platform", "m/z", "Retention time", "Raw intensity",
      "Log intensity", "Z-scaled intensity", "P-value", "Fold-change", "Curated Log2FC"
    )))
  )

study_completeness_score <- study_data_availability_long %>%
  group_by(study_id) %>%
  summarise(
    n_available = sum(present, na.rm = TRUE),
    n_possible = n(),
    completeness_fraction = n_available / n_possible,
    .groups = "drop"
  ) %>%
  arrange(desc(completeness_fraction))

write_csv(study_data_availability, file.path(out_dir, "study_data_availability_flags.csv"))
write_csv(study_completeness_score, file.path(out_dir, "study_data_completeness_score.csv"))

plots$study_data_availability_heatmap <- study_data_availability_long %>%
  ggplot(aes(x = study_id, y = info_type, fill = present_label)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  labs(title = "What information is available for each study?",
       x = "Study", y = NULL, fill = NULL) +
  scale_fill_manual(values = present_palette, na.value = col_grey) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9))

plots$study_data_completeness_score <- study_completeness_score %>%
  mutate(study_id = reorder(study_id, completeness_fraction)) %>%
  ggplot(aes(x = study_id, y = completeness_fraction)) +
  geom_col(fill = col_teal) +
  geom_text(aes(label = paste0(round(100 * completeness_fraction), "%")),
            hjust = -0.1, size = 3.2) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.05)) +
  labs(title = "Study-level data completeness score",
       x = "Study", y = "Information present") +
  plot_theme

## B. Top recurrent metabolites/entities across all studies.
entity_recurrence_overall <- wide %>%
  mutate(
    n_studies_num = safe_num(n_studies),
    annotation_score_num = safe_num(annotation_score),
    label = str_trunc(coalesce(canonical_metabolite, entity_id), 45),
    recurrent_status = if_else(!is.na(n_studies_num) & n_studies_num >= 2,
                               "Recurrent (>=2 studies)", "Not recurrent")
  ) %>%
  arrange(desc(n_studies_num), desc(annotation_score_num), label)

write_csv(entity_recurrence_overall, file.path(out_dir, "entity_recurrence_overall.csv"))

plots$top_recurrent_metabolites_all <- entity_recurrence_overall %>%
  filter(!is.na(n_studies_num), n_studies_num >= 2) %>%
  slice_max(n_studies_num, n = 25, with_ties = FALSE) %>%
  mutate(label = reorder(label, n_studies_num)) %>%
  ggplot(aes(x = label, y = n_studies_num)) +
  geom_col(fill = col_teal) +
  geom_text(aes(label = n_studies_num), hjust = -0.1, size = 3.1) +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, max(entity_recurrence_overall$n_studies_num, na.rm = TRUE), by = 1),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Top recurrent metabolite entities across studies",
      x = NULL, y = "Number of studies") +
  plot_theme

## C. Ratio of recurrent versus non-recurrent entities.
recurrent_ratio <- entity_recurrence_overall %>%
  count(recurrent_status, name = "n") %>%
  mutate(percent = n / sum(n))
write_csv(recurrent_ratio, file.path(out_dir, "recurrent_vs_not_recurrent_ratio.csv"))

plots$recurrent_vs_not_recurrent_ratio <- recurrent_ratio %>%
  mutate(recurrent_status = reorder(recurrent_status, n)) %>%
  ggplot(aes(x = recurrent_status, y = n)) +
  geom_col(fill = col_teal) +
  geom_text(aes(label = paste0(n, " (", round(100 * percent, 1), "%)")),
            hjust = -0.05, size = 3.5) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Ratio of recurrent to non-recurrent metabolite entities",
       x = NULL, y = "Number of entities") +
  plot_theme

## D. Recurrent up/down entities ranked by recurrence first, then effect size.
top_recurrent_up <- entity_consistency %>%
  filter(recurrent, dominant_direction == "Up in PD") %>%
  arrange(desc(n_studies_with_log2fc), desc(mean_log2fc)) %>%
  slice_head(n = 25) %>%
  mutate(label = reorder(str_trunc(canonical_metabolite, 45), mean_log2fc))

top_recurrent_down <- entity_consistency %>%
  filter(recurrent, dominant_direction == "Down in PD") %>%
  arrange(desc(n_studies_with_log2fc), mean_log2fc) %>%
  slice_head(n = 25) %>%
  mutate(label = reorder(str_trunc(canonical_metabolite, 45), -mean_log2fc))

write_csv(top_recurrent_up, file.path(out_dir, "top_recurrent_up_entities.csv"))
write_csv(top_recurrent_down, file.path(out_dir, "top_recurrent_down_entities.csv"))

plots$top_recurrent_up_by_mean_log2fc <- if (nrow(top_recurrent_up) > 0) {
  top_recurrent_up %>%
    ggplot(aes(x = label, y = mean_log2fc)) +
    geom_col(fill = col_red) +
    geom_text(aes(label = paste0("n=", n_studies_with_log2fc)), hjust = -0.1, size = 3) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Top recurrent up-regulated entities",
          x = NULL, y = "Mean Log2FC") +
    plot_theme
} else blank_plot("Top recurrent up-regulated entities")

plots$top_recurrent_down_by_mean_log2fc <- if (nrow(top_recurrent_down) > 0) {
  top_recurrent_down %>%
    ggplot(aes(x = label, y = mean_log2fc)) +
    geom_col(fill = col_blue) +
    geom_text(aes(label = paste0("n=", n_studies_with_log2fc)), hjust = 1.1, size = 3) +
    coord_flip() +
    labs(title = "Top recurrent down-regulated entities",
          x = NULL, y = "Mean Log2FC") +
    plot_theme
} else blank_plot("Top recurrent down-regulated entities")

## E. Ratio of non-recurrent, consistent recurrent, and discordant recurrent entities.
consistency_ratio_all <- entity_consistency %>%
  mutate(consistency_category = case_when(
    !recurrent ~ "Not recurrent with Log2FC",
    recurrent & discordant ~ "Recurrent discordant",
    recurrent & !discordant ~ "Recurrent consistent",
    TRUE ~ "Unclassified"
  )) %>%
  count(consistency_category, name = "n") %>%
  mutate(percent = n / sum(n))
write_csv(consistency_ratio_all, file.path(out_dir, "consistency_discordance_ratio.csv"))

plots$consistency_discordance_ratio <- consistency_ratio_all %>%
  mutate(consistency_category = reorder(consistency_category, n)) %>%
  ggplot(aes(x = consistency_category, y = n)) +
  geom_col(fill = col_teal) +
  geom_text(aes(label = paste0(n, " (", round(100 * percent, 1), "%)")),
            hjust = -0.05, size = 3.5) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Consistency and discordance among entities with Log2FC",
       x = NULL, y = "Number of entities") +
  plot_theme

## F. Top discordant entities: show up-study and down-study counts directly.
top_discordant_counts <- entity_consistency %>%
  filter(recurrent, discordant) %>%
  arrange(desc(n_studies_with_log2fc), consistency_score, desc(max_abs_log2fc)) %>%
  slice_head(n = 25) %>%
  transmute(
    entity_id,
    label = str_trunc(canonical_metabolite, 45),
    `Up in PD` = n_up,
    `Down in PD` = n_down,
    consistency_score,
    n_studies_with_log2fc
  )
write_csv(top_discordant_counts, file.path(out_dir, "top_discordant_entities_up_down_counts.csv"))

top_discordant_long <- top_discordant_counts %>%
  pivot_longer(cols = c(`Up in PD`, `Down in PD`), names_to = "direction", values_to = "n_direction")

plots$top_discordant_up_down_counts <- if (nrow(top_discordant_long) > 0) {
  top_discordant_long %>%
    mutate(label = reorder(label, n_studies_with_log2fc)) %>%
    ggplot(aes(x = label, y = n_direction, fill = direction)) +
    geom_col(position = "stack") +
    coord_flip() +
    labs(title = "Top discordant recurrent entities",
         x = NULL, y = "Number of studies", fill = "Direction") +
    scale_fill_manual(values = direction_palette, na.value = col_grey) +
    plot_theme
} else blank_plot("Top discordant recurrent entities", "No recurrent discordant entities found")

## G. PCA and MDS of metabolite-study observations.
## Points are entity-study observations. Features are curated effect/measurement variables.
## -------------------------------
## PCA/MDS ordination
## -------------------------------
## Main biological/effect ordination:
## Uses only effect-size/statistical-strength variables.
## m/z and RT are excluded here because they mainly reflect analytical
## feature properties and platform/study-method differences.

effect_features <- c(
  "Log2FC",
  "z_PD_minus_control",
  "neg_log10_p"
)

effect_features <- intersect(effect_features, names(long))

ordination_df <- long %>%
  select(
    entity_id, canonical_metabolite, study_id,
    matrix_clean, mode_clean, phase_clean,
    all_of(effect_features)
  ) %>%
  mutate(across(all_of(effect_features), safe_num)) %>%
  mutate(non_missing_features = rowSums(!is.na(across(all_of(effect_features))))) %>%
  filter(non_missing_features >= 2)

if (nrow(ordination_df) >= 10 && length(effect_features) >= 2) {
  
  num_df <- ordination_df %>%
    select(all_of(effect_features))
  
  ## Keep only numeric features with enough values and variation.
  keep_feature <- vapply(
    num_df,
    function(x) sum(!is.na(x)) >= 10 && sd(x, na.rm = TRUE) > 0,
    logical(1)
  )
  
  num_df <- num_df[, keep_feature, drop = FALSE]
  
  if (ncol(num_df) >= 2) {
    
    ## Median-impute missing values for exploratory ordination.
    for (jj in seq_along(num_df)) {
      med <- median(num_df[[jj]], na.rm = TRUE)
      num_df[[jj]][is.na(num_df[[jj]])] <- med
    }
    
    ## Scale features so PCA/MDS are not dominated by variable ranges.
    X <- scale(as.matrix(num_df))
    
    keep_rows <- apply(X, 1, function(z) all(is.finite(z)))
    X <- X[keep_rows, , drop = FALSE]
    ord_meta <- ordination_df[keep_rows, , drop = FALSE]
    
    ## -------------------------------
    ## PCA
    ## -------------------------------
    pca_fit <- prcomp(X, center = FALSE, scale. = FALSE)
    
    pca_scores <- as_tibble(pca_fit$x[, 1:2, drop = FALSE]) %>%
      setNames(c("PC1", "PC2")) %>%
      bind_cols(
        ord_meta %>%
          select(entity_id, canonical_metabolite, study_id,
                 matrix_clean, mode_clean, phase_clean)
      )
    
    pca_var <- round(
      100 * (pca_fit$sdev^2 / sum(pca_fit$sdev^2))[1:2],
      1
    )
    
    write_csv(
      pca_scores,
      file.path(out_dir, "pca_effect_metabolite_study_scores.csv")
    )
    
    plots$pca_effect_by_study <- pca_scores %>%
      ggplot(aes(x = PC1, y = PC2, colour = study_id)) +
      geom_point(alpha = 0.45, size = 1.2) +
      labs(
        title = "PCA of reported PD-control effect patterns by study",
        x = paste0("PC1 (", pca_var[1], "%)"),
        y = paste0("PC2 (", pca_var[2], "%)"),
        colour = "Study"
      ) +
      plot_theme +
      theme(axis.text.x = element_text(angle = 0))
    
    plots$pca_effect_by_matrix <- pca_scores %>%
      ggplot(aes(x = PC1, y = PC2, colour = matrix_clean)) +
      geom_point(alpha = 0.45, size = 1.2) +
      labs(
        title = "PCA of reported PD-control effect patterns by sample matrix",
        x = paste0("PC1 (", pca_var[1], "%)"),
        y = paste0("PC2 (", pca_var[2], "%)"),
        colour = "Matrix"
      ) +
      scale_colour_manual(values = matrix_palette, na.value = col_grey) +
      plot_theme +
      theme(axis.text.x = element_text(angle = 0))
    
    ## -------------------------------
    ## MDS
    ## -------------------------------
    ## MDS uses Euclidean distance on the scaled effect matrix.
    ## Sampling is capped to avoid very large distance matrices.
    
    set.seed(7)
    mds_max_n <- 2500
    
    mds_idx <- if (nrow(X) > mds_max_n) {
      sample(seq_len(nrow(X)), mds_max_n)
    } else {
      seq_len(nrow(X))
    }
    
    mds_fit <- cmdscale(
      dist(X[mds_idx, , drop = FALSE]),
      k = 2
    )
    
    mds_scores <- as_tibble(mds_fit) %>%
      setNames(c("MDS1", "MDS2")) %>%
      bind_cols(
        ord_meta[mds_idx, ] %>%
          select(entity_id, canonical_metabolite, study_id,
                 matrix_clean, mode_clean, phase_clean)
      )
    
    write_csv(
      mds_scores,
      file.path(out_dir, "mds_effect_metabolite_study_scores.csv")
    )
    
    plots$mds_effect_by_study <- mds_scores %>%
      ggplot(aes(x = MDS1, y = MDS2, colour = study_id)) +
      geom_point(alpha = 0.45, size = 1.2) +
      labs(
        title = "MDS of reported PD-control effect patterns by study",
        x = "MDS1",
        y = "MDS2",
        colour = "Study"
      ) +
      plot_theme +
      theme(axis.text.x = element_text(angle = 0))
    
    plots$mds_effect_by_matrix <- mds_scores %>%
      ggplot(aes(x = MDS1, y = MDS2, colour = matrix_clean)) +
      geom_point(alpha = 0.45, size = 1.2) +
      labs(
        title = "MDS of reported PD-control effect patterns by sample matrix",
        x = "MDS1",
        y = "MDS2",
        colour = "Matrix"
      ) +
      scale_colour_manual(values = matrix_palette, na.value = col_grey) +
      plot_theme +
      theme(axis.text.x = element_text(angle = 0))
    
  } else {
    
    plots$pca_effect_by_study <- blank_plot(
      "PCA of reported PD-control effect patterns by study",
      "Not enough variable numeric effect features"
    )
    
    plots$pca_effect_by_matrix <- blank_plot(
      "PCA of reported PD-control effect patterns by matrix",
      "Not enough variable numeric effect features"
    )
    
    plots$mds_effect_by_study <- blank_plot(
      "MDS of reported PD-control effect patterns by study",
      "Not enough variable numeric effect features"
    )
    
    plots$mds_effect_by_matrix <- blank_plot(
      "MDS of reported PD-control effect patterns by matrix",
      "Not enough variable numeric effect features"
    )
  }
  
} else {
  
  plots$pca_effect_by_study <- blank_plot(
    "PCA of reported PD-control effect patterns by study",
    "Not enough observations with at least two numeric effect features"
  )
  
  plots$pca_effect_by_matrix <- blank_plot(
    "PCA of reported PD-control effect patterns by matrix",
    "Not enough observations with at least two numeric effect features"
  )
  
  plots$mds_effect_by_study <- blank_plot(
    "MDS of reported PD-control effect patterns by study",
    "Not enough observations with at least two numeric effect features"
  )
  
  plots$mds_effect_by_matrix <- blank_plot(
    "MDS of reported PD-control effect patterns by matrix",
    "Not enough observations with at least two numeric effect features"
  )
}


## -------------------------------
## Technical/analytical QC plots
## -------------------------------
## m/z and RT are useful, but better treated as analytical coverage/QC
## variables instead of being mixed into the main biological ordination.

technical_features <- c("mz", "RT_min", "annotation_score")
technical_features <- intersect(technical_features, names(long))

tech_df <- long %>%
  select(
    entity_id, canonical_metabolite, study_id,
    matrix_clean, mode_clean, phase_clean,
    all_of(technical_features)
  ) %>%
  mutate(across(all_of(technical_features), safe_num))

if ("mz" %in% names(tech_df)) {
  
  mz_df <- tech_df %>%
    filter(!is.na(mz))
  
  if (nrow(mz_df) >= 10) {
    
    plots$mz_distribution_by_study <- mz_df %>%
      ggplot(aes(x = mz, fill = study_id)) +
      geom_density(alpha = 0.30) +
      labs(
        title = "m/z distribution across studies",
        x = "m/z",
        y = "Density",
        fill = "Study"
      ) +
      plot_theme +
      theme(axis.text.x = element_text(angle = 0))
    
    plots$mz_distribution_by_matrix <- mz_df %>%
      ggplot(aes(x = mz, fill = matrix_clean)) +
      geom_density(alpha = 0.35) +
      labs(
        title = "m/z distribution across sample matrices",
        x = "m/z",
        y = "Density",
        fill = "Matrix"
      ) +
      scale_fill_manual(values = matrix_palette, na.value = col_grey) +
      plot_theme +
      theme(axis.text.x = element_text(angle = 0))
    
  } else {
    
    plots$mz_distribution_by_study <- blank_plot(
      "m/z distribution across studies",
      "Not enough m/z values"
    )
    
    plots$mz_distribution_by_matrix <- blank_plot(
      "m/z distribution across sample matrices",
      "Not enough m/z values"
    )
  }
}

if ("RT_min" %in% names(tech_df)) {
  
  rt_df <- tech_df %>%
    filter(!is.na(RT_min))
  
  if (nrow(rt_df) >= 10) {
    
    plots$rt_distribution_by_study <- rt_df %>%
      ggplot(aes(x = RT_min, fill = study_id)) +
      geom_density(alpha = 0.30) +
      labs(
        title = "Retention-time distribution across studies",
        x = "Retention time, min",
        y = "Density",
        fill = "Study"
      ) +
      plot_theme +
      theme(axis.text.x = element_text(angle = 0))
    
    plots$rt_distribution_by_matrix <- rt_df %>%
      ggplot(aes(x = RT_min, fill = matrix_clean)) +
      geom_density(alpha = 0.35) +
      labs(
        title = "Retention-time distribution across sample matrices",
        x = "Retention time, min",
        y = "Density",
        fill = "Matrix"
      ) +
      scale_fill_manual(values = matrix_palette, na.value = col_grey) +
      plot_theme +
      theme(axis.text.x = element_text(angle = 0))
    
  } else {
    
    plots$rt_distribution_by_study <- blank_plot(
      "Retention-time distribution across studies",
      "Not enough RT values"
    )
    
    plots$rt_distribution_by_matrix <- blank_plot(
      "Retention-time distribution across sample matrices",
      "Not enough RT values"
    )
  }
}

if (all(c("mz", "RT_min") %in% names(tech_df))) {
  
  mz_rt_df <- tech_df %>%
    filter(!is.na(mz), !is.na(RT_min))
  
  if (nrow(mz_rt_df) >= 10) {
    
    plots$mz_rt_space_by_study <- mz_rt_df %>%
      ggplot(aes(x = RT_min, y = mz, colour = study_id)) +
      geom_point(alpha = 0.45, size = 1.2) +
      labs(
        title = "Analytical feature space by study",
        x = "Retention time, min",
        y = "m/z",
        colour = "Study"
      ) +
      plot_theme +
      theme(axis.text.x = element_text(angle = 0))
    
    plots$mz_rt_space_by_matrix <- mz_rt_df %>%
      ggplot(aes(x = RT_min, y = mz, colour = matrix_clean)) +
      geom_point(alpha = 0.45, size = 1.2) +
      labs(
        title = "Analytical feature space by sample matrix",
        x = "Retention time, min",
        y = "m/z",
        colour = "Matrix"
      ) +
      scale_colour_manual(values = matrix_palette, na.value = col_grey) +
      plot_theme +
      theme(axis.text.x = element_text(angle = 0))
    
  } else {
    
    plots$mz_rt_space_by_study <- blank_plot(
      "Analytical feature space by study",
      "Not enough observations with both m/z and RT"
    )
    
    plots$mz_rt_space_by_matrix <- blank_plot(
      "Analytical feature space by sample matrix",
      "Not enough observations with both m/z and RT"
    )
  }
}
## ---- 7. Save plots ----------------------------------------------------------
pdf(file.path(out_dir, "pd_met_plos_exploratory_plots.pdf"), width = 12, height = 8)
for (nm in names(plots)) {
  print(plots[[nm]])
  save_plot(plots[[nm]], nm)
}
dev.off()

## ---- 8. Console summary -----------------------------------------------------
message("\nDone.")
message("Output directory: ", normalizePath(out_dir))
message("PDF: ", file.path(out_dir, "pd_met_plos_exploratory_plots.pdf"))
message("Core tables:")
message(" - long_study_metabolite_table.csv")
message(" - inchikey_entity_up_down_consistency.csv")
message(" - study_summary.csv")
message(" - discordant_recurrent_entities.csv")
message(" - consistent_recurrent_entities.csv")
