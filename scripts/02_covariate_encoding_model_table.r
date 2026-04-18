set.seed(405)
options(stringsAsFactors = FALSE)

pkgs <- c("gssr", "dplyr", "forcats", "readr", "tibble")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need) > 0) install.packages(need, repos = "https://cloud.r-project.org")
for (p in pkgs) library(p, character.only = TRUE)
g <- gssr::gss_get_yr(2024)
keep <- c("id", "confinan", "polviews", "degree", "age", "sex", "income", "region")

raw <- g %>% dplyr::select(dplyr::any_of(keep))

tmp <- raw %>%

  dplyr::mutate(
    confinan = as.character(confinan),
    polviews = suppressWarnings(as.numeric(as.character(polviews))),
    degree = forcats::fct_infreq(as.factor(degree)),
    sex = as.factor(sex),
    income = suppressWarnings(as.numeric(as.character(income))),
    region = as.factor(region)
  )

tmp$confinan_ord <- ifelse(tmp$confinan %in% c("1", "A GREAT DEAL"), 1L,
  ifelse(tmp$confinan %in% c("2", "ONLY SOME"), 2L,
    ifelse(tmp$confinan %in% c("3", "HARDLY ANY"), 3L, NA_integer_)
  )
)

tmp$polviews_bin <- ifelse(!is.na(tmp$polviews) & tmp$polviews >= 1 & tmp$polviews <= 3, "Liberal",
  ifelse(!is.na(tmp$polviews) & tmp$polviews == 4, "Moderate",
    ifelse(!is.na(tmp$polviews) & tmp$polviews >= 5 & tmp$polviews <= 7, "Conservative", NA_character_)
  )
)

tmp$polviews_bin <- factor(tmp$polviews_bin, levels = c("Liberal", "Moderate", "Conservative"))
tmp$age_std <- as.numeric(scale(tmp$age))
tmp$income_std <- as.numeric(scale(tmp$income))

model_tbl <- tmp %>%
  dplyr::select(
    id, confinan_ord, polviews, polviews_bin, degree, age, age_std, sex, income, income_std, region
  ) %>%
  dplyr::filter(!is.na(confinan_ord), !is.na(polviews_bin), !is.na(age_std), !is.na(income_std))

encoding_audit <- tibble::tibble(
  metric = c("n_raw", "n_model_tbl", "n_removed", "pct_removed"),
  value = c(
    nrow(raw),
    nrow(model_tbl),
    nrow(raw) - nrow(model_tbl),
    (nrow(raw) - nrow(model_tbl)) / nrow(raw)
  )
)

encoding_audit

dplyr::count(model_tbl, polviews_bin)
dplyr::count(model_tbl, confinan_ord)
dir.create("../output/derived", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(model_tbl, "../output/derived/gss_2024_model_table.csv")
readr::write_csv(encoding_audit, "../output/derived/encoding_audit.csv")

