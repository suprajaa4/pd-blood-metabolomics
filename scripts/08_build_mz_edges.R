suppressPackageStartupMessages({library(readxl); library(dplyr); library(purrr)})

data_path <- "data/raw/pd_with_adducts.xlsx"
main <- read_excel(data_path, sheet = "pd-with-adducts")

studies <- c('S_24058461','S_24167579','S_25390405','S_28880465','S_29518718','S_33485385',
             'S_34577635','S_35829770','S_37335671','S_37589832','S_37755270','S_38516193',
             'S_38720721','ACS1517','MTBLS1219','MTBLS11904','ST001814','ST003159')

#get from cols
get_mz_values <- function(s) {
  col <- paste0(s, "__mz")
  if (!(col %in% colnames(main))) return(NULL)
  vals <- main[[col]]
  parsed <- sapply(vals, function(v) {
    if (is.na(v)) return(NA_real_)
    nums <- suppressWarnings(as.numeric(trimws(strsplit(as.character(v), "\\|")[[1]])))
    nums <- nums[!is.na(nums)]
    if (length(nums) == 0) NA_real_ else mean(nums)
  })
  data.frame(entity_id = main$entity_id, mz = parsed)
}
mz_long <- bind_rows(map(studies, get_mz_values)) %>% filter(!is.na(mz))
mz_rep <- mz_long %>% group_by(entity_id) %>% summarise(mz = median(mz), .groups = "drop")
cat("Metabolites with a representative m/z:", nrow(mz_rep), "\n")


mz_rep$mz_bin <- round(mz_rep$mz)

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

edges_mz <- build_weighted_edges(mz_rep, "entity_id", "mz_bin")
write.csv(edges_mz, "data/processed/edges_mz.csv", row.names = FALSE)
cat("edges_mz.csv:", nrow(edges_mz), "edges among",
    n_distinct(c(edges_mz$from, edges_mz$to)), "metabolites\n")
