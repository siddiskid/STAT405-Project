set.seed(405)
options(stringsAsFactors = FALSE)

pkgs <- c("gssr", "dplyr", "tidyr", "readr", "tibble")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need) > 0) install.packages(need, repos = "https://cloud.r-project.org")
for (p in pkgs) library(p, character.only = TRUE)

g <- gssr::gss_get_yr(2024)

keep <- c("id", "confinan", "polviews", "degree", "age", "sex", "income", "region")

d <- g %>% dplyr::select(dplyr::any_of(keep))

dplyr::glimpse(d)

missingness <- d %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~sum(is.na(.)))) %>%
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "variable", values_to = "n_missing") %>%
  dplyr::mutate(pct_missing = n_missing / nrow(d)) %>%
  dplyr::arrange(dplyr::desc(pct_missing))

missingness

v1 <- d %>% dplyr::count(confinan, name = "n") %>% dplyr::mutate(variable = "confinan", level = as.character(confinan))
v2 <- d %>% dplyr::count(polviews, name = "n") %>% dplyr::mutate(variable = "polviews", level = as.character(polviews))
v3 <- d %>% dplyr::count(degree, name = "n") %>% dplyr::mutate(variable = "degree", level = as.character(degree))
v4 <- d %>% dplyr::count(sex, name = "n") %>% dplyr::mutate(variable = "sex", level = as.character(sex))
v5 <- d %>% dplyr::count(income, name = "n") %>% dplyr::mutate(variable = "income", level = as.character(income))
v6 <- d %>% dplyr::count(region, name = "n") %>% dplyr::mutate(variable = "region", level = as.character(region))

value_counts <- dplyr::bind_rows(v1, v2, v3, v4, v5, v6)
value_counts <- value_counts %>% dplyr::select(variable, level, n) %>% dplyr::arrange(variable, dplyr::desc(n))

value_counts

audit_summary <- tibble::tibble(
  metric = c("n_rows_raw", "n_cols_raw", "n_rows_core", "n_complete_core"),
  value = c(nrow(g), ncol(g), nrow(d), sum(stats::complete.cases(d)))
)

audit_summary

dir.create("../output/audit", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(missingness, "../output/audit/missingness_2024_core.csv")
readr::write_csv(value_counts, "../output/audit/value_counts_2024_core.csv")
readr::write_csv(audit_summary, "../output/audit/audit_summary_2024_core.csv")
