suppressPackageStartupMessages({library(readxl); library(dplyr)})
main <- read_excel("data/raw/pd_with_adducts.xlsx", sheet = "pd-with-adducts")
classyfire <- read_excel("data/raw/pd_with_adducts.xlsx", sheet = "class of compounds")
classyfire_dedup <- classyfire %>% select(InChIKey, Class) %>% distinct(InChIKey, .keep_all = TRUE)

met_class <- main %>% select(entity_id, InChIKey) %>%
  left_join(classyfire_dedup, by = "InChIKey") %>% filter(!is.na(Class))

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

edges_class <- build_weighted_edges(met_class, "entity_id", "Class")
write.csv(edges_class, "data/processed/edges_class.csv", row.names = FALSE)
cat("edges_class.csv:", nrow(edges_class), "edges\n")
