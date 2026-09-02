suppressPackageStartupMessages({library(dplyr); library(stringr); library(tidyr); library(readr)})
TOLERANCE_DA <- 0.02
input_file <- "data/raw/metabolites_grouped_direction_corrected.csv"
out_dir <- "data/processed"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
adducts <- tribble(
  ~matched_adduct, ~mass_shift, ~charge, ~multiplier,
  "[M-H]-", -1.007276, -1, 1,
  "[M+H]+", 1.007276, 1, 1,
  "[M+FA-H]- / [M+HCOO]-", 44.998201, -1, 1,
  "[M-H2O-H]-", -19.017841, -1, 1,
  "[M+H-H2O]+", -17.003289, 1, 1,
  "[M+NH4]+", 18.033823, 1, 1,
  "[M+Na]+", 22.989218, 1, 1,
  "[M+K]+", 38.963158, 1, 1,
  "[M+CH3COO]-", 59.013851, -1, 1,
  "[2M-H]-", -1.007276, -1, 2,
  "[2M+H]+", 1.007276, 1, 2,
  "[2M+Na]+", 22.989218, 1, 2,
  "[M+2H]2+", 2.014552, 2, 1,
  "[M-2H]2-", -2.014552, -2, 1
)

atomic_mass <- c(C=12.000000, H=1.007825, N=14.003074, O=15.994915,
                 P=30.973762, S=31.972071, F=18.998403, Cl=34.968853,
                 Br=78.918338, I=126.904473, Na=22.989770, K=38.963707,
                 Si=27.976927)

formula_mass <- function(formula) {
  if (is.na(formula) || formula == "") return(NA_real_)
  parts <- str_match_all(formula, "([A-Z][a-z]?)(\\d*)")[[1]]
  if (nrow(parts) == 0) return(NA_real_)
  vals <- mapply(function(el, n) {
    if (!el %in% names(atomic_mass)) return(NA_real_)
    atomic_mass[[el]] * ifelse(n == "", 1, as.numeric(n))
  }, parts[,2], parts[,3])
  if (any(is.na(vals))) NA_real_ else sum(vals)
}
parse_numeric_cell <- function(x) {
  if (is.na(x)) return(numeric(0))
  vals <- suppressWarnings(as.numeric(str_trim(str_split(as.character(x), "\\|")[[1]])))
  vals[!is.na(vals)]
}
df <- read_csv(input_file, show_col_types = FALSE, guess_max = 100000)
mz_cols <- names(df)[str_detect(names(df), "__mz$")]
df <- df %>%
  rowwise() %>%
  mutate(
    neutral_mass = formula_mass(formula),
    observed_mz = {
      vals <- unlist(lapply(c_across(all_of(mz_cols)), parse_numeric_cell))
      if (length(vals) == 0) NA_real_ else median(vals)
    }
  ) %>%
  ungroup()
matches <- df %>%
  select(entity_id, canonical_metabolite, formula, neutral_mass, observed_mz) %>%
  filter(!is.na(neutral_mass), !is.na(observed_mz)) %>%
  crossing(adducts) %>%
  mutate(
    theoretical_mz = (multiplier * neutral_mass + mass_shift) / abs(charge),
    diff_da = observed_mz - theoretical_mz,
    abs_diff_da = abs(diff_da)
  ) %>%
  filter(abs_diff_da <= TOLERANCE_DA) %>%
  arrange(entity_id, abs_diff_da)

best <- matches %>% group_by(entity_id) %>% slice_min(abs_diff_da, n = 1, with_ties = FALSE) %>% ungroup()
write_csv(matches, file.path(out_dir, "adduct_matches_all.csv"))
write_csv(best, file.path(out_dir, "adduct_matches_best.csv"))

