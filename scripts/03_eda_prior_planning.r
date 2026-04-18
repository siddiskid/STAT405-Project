set.seed(405)
options(stringsAsFactors = FALSE)

pkgs <- c("dplyr", "ggplot2", "readr", "tibble", "tidyr", "forcats")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need) > 0) install.packages(need, repos = "https://cloud.r-project.org")
for (p in pkgs) library(p, character.only = TRUE)

d <- readr::read_csv("../output/derived/gss_2024_model_table.csv", show_col_types = FALSE)

d <- d %>%

  dplyr::mutate(
    confinan_ord = as.integer(confinan_ord),
    polviews_bin = factor(polviews_bin, levels = c("Liberal", "Moderate", "Conservative")),
    degree = as.factor(degree),
    sex = as.factor(sex),
    region = as.factor(region)
  )

sample_overview <- tibble::tibble(
  metric = c("n_rows", "n_cols", "age_mean", "age_sd", "income_mean", "income_sd"),
  value = c(
    nrow(d),
    ncol(d),
    mean(d$age, na.rm = TRUE),
    sd(d$age, na.rm = TRUE),
    mean(d$income, na.rm = TRUE),
    sd(d$income, na.rm = TRUE)
  )
)

sample_overview

outcome_dist <- d %>%
  dplyr::count(confinan_ord, name = "n") %>%
  dplyr::mutate(prop = n / sum(n))

outcome_dist

outcome_by_ideo <- d %>%
  dplyr::count(polviews_bin, confinan_ord, name = "n") %>%
  dplyr::group_by(polviews_bin) %>%
  dplyr::mutate(prop_within_ideo = n / sum(n)) %>%
  dplyr::ungroup()

outcome_by_ideo

tmp_age <- d %>%
  dplyr::group_by(confinan_ord) %>%
  dplyr::summarise(
    n = dplyr::n(),
    age_mean = mean(age, na.rm = TRUE),
    age_sd = sd(age, na.rm = TRUE),
    income_mean = mean(income, na.rm = TRUE),
    income_sd = sd(income, na.rm = TRUE),
    .groups = "drop"
  )

age_income_by_outcome <- tmp_age

age_income_by_outcome

p_outcome <- ggplot2::ggplot(outcome_dist, ggplot2::aes(x = factor(confinan_ord), y = prop)) +
  ggplot2::geom_col() +
  ggplot2::labs(x = "confinan_ord", y = "Proportion")

p_outcome

p_ideo <- ggplot2::ggplot(outcome_by_ideo, ggplot2::aes(x = polviews_bin, y = prop_within_ideo, fill = factor(confinan_ord))) +
  ggplot2::geom_col(position = "stack") +
  ggplot2::labs(x = "Ideology group", y = "Within-group proportion", fill = "confinan_ord")

p_ideo

prior_plan <- data.frame(
  parameter_block = c(
    "cutpoints",
    "beta_std_predictors",
    "beta_categorical",
    "hier_alpha_mu",
    "hier_alpha_sigma"
  ),
  baseline_prior = c(
    "normal(0, 2.5)",
    "normal(0, 1)",
    "normal(0, 1)",
    "normal(0, 1)",
    "normal+(0, 1)"
  ),
  sensitivity_prior_1 = c(
    "normal(0, 1.5)",
    "normal(0, 0.5)",
    "normal(0, 0.5)",
    "normal(0, 0.5)",
    "normal+(0, 0.5)"
  ),
  sensitivity_prior_2 = c(
    "normal(0, 5)",
    "normal(0, 2)",
    "normal(0, 2)",
    "normal(0, 2)",
    "normal+(0, 2)"
  ),
  rationale = c(
    "Outcome has 3 levels; moderate scale on latent logit axis",
    "age_std and income_std are standardized",
    "Regularize dummy coefficients to avoid overfit",
    "Centered ideology-level baseline on logit scale",
    "Controls between-ideology dispersion"
  ),
  stringsAsFactors = FALSE
)

prior_plan
dir.create("../output/eda", recursive = TRUE, showWarnings = FALSE)
dir.create("../output/prior", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(sample_overview, "../output/eda/sample_overview.csv")
readr::write_csv(outcome_dist, "../output/eda/outcome_distribution.csv")
readr::write_csv(outcome_by_ideo, "../output/eda/outcome_by_ideology.csv")
readr::write_csv(age_income_by_outcome, "../output/eda/age_income_by_outcome.csv")
readr::write_csv(prior_plan, "../output/prior/prior_plan.csv")
ggplot2::ggsave(filename = "../output/eda/outcome_distribution.png", plot = p_outcome, width = 7, height = 4, dpi = 150)
ggplot2::ggsave(filename = "../output/eda/outcome_by_ideology.png", plot = p_ideo, width = 8, height = 4.5, dpi = 150)

