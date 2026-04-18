set.seed(405)

options(stringsAsFactors = FALSE)

pkgs <- c("rstan", "loo", "dplyr", "readr", "tibble", "ggplot2", "bayesplot", "tidyr")

need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(need) > 0) install.packages(need, repos = "https://cloud.r-project.org")

for (p in pkgs) library(p, character.only = TRUE)

if (!file.exists("../output/model1_nuts/fit_m1_nuts.rds")) stop("Missing model1 fit")

if (!file.exists("../output/model2_nuts/fit_m2_nuts.rds")) stop("Missing model2 fit")

fit_m1 <- readRDS("../output/model1_nuts/fit_m1_nuts.rds")

fit_m2 <- readRDS("../output/model2_nuts/fit_m2_nuts.rds")

d <- readr::read_csv("../output/derived/gss_2024_model_table.csv", show_col_types = FALSE) %>%

  dplyr::mutate(

    degree = as.factor(degree),

    sex = as.factor(sex),

    region = as.factor(region),

    polviews_bin = factor(polviews_bin, levels = c("Liberal", "Moderate", "Conservative"))

  )

m1 <- d %>%

  dplyr::filter(!is.na(confinan_ord), !is.na(polviews), !is.na(age_std), !is.na(income_std), !is.na(degree), !is.na(sex), !is.na(region))

m2 <- d %>%

  dplyr::filter(!is.na(confinan_ord), !is.na(polviews_bin), !is.na(age_std), !is.na(income_std), !is.na(degree), !is.na(sex), !is.na(region))

y_m1 <- as.integer(m1$confinan_ord)

y_m2 <- as.integer(m2$confinan_ord)

ideology_m2 <- m2$polviews_bin

log_lik_m1 <- rstan::extract(fit_m1, pars = "log_lik", permuted = TRUE)$log_lik

log_lik_m2 <- rstan::extract(fit_m2, pars = "log_lik", permuted = TRUE)$log_lik

if (is.null(log_lik_m1) || nrow(log_lik_m1) == 0) stop("Missing log_lik m1")

if (is.null(log_lik_m2) || nrow(log_lik_m2) == 0) stop("Missing log_lik m2")

loo_m1 <- loo::loo(log_lik_m1)

loo_m2 <- loo::loo(log_lik_m2)

lc <- loo::loo_compare(loo_m1, loo_m2)

lc

loo_summary <- tibble::tibble(

  model = c("Model1", "Model2"),

  elpd_loo = c(loo_m1$estimates["elpd_loo", "Estimate"], loo_m2$estimates["elpd_loo", "Estimate"]),

  se_elpd_loo = c(loo_m1$estimates["elpd_loo", "SE"], loo_m2$estimates["elpd_loo", "SE"]),

  p_loo = c(loo_m1$estimates["p_loo", "Estimate"], loo_m2$estimates["p_loo", "Estimate"]),

  looic = c(loo_m1$estimates["looic", "Estimate"], loo_m2$estimates["looic", "Estimate"])

) %>%

  dplyr::arrange(dplyr::desc(elpd_loo))

elpd_diff_m2_minus_m1 <- loo_m2$estimates["elpd_loo", "Estimate"] - loo_m1$estimates["elpd_loo", "Estimate"]

se_diff_approx <- sqrt(loo_m1$estimates["elpd_loo", "SE"]^2 + loo_m2$estimates["elpd_loo", "SE"]^2)

loo_diff <- tibble::tibble(metric = c("elpd_diff_m2_minus_m1", "se_diff_approx"), value = c(elpd_diff_m2_minus_m1, se_diff_approx))

loo_summary

yrep_m1 <- rstan::extract(fit_m1, pars = "y_rep", permuted = TRUE)$y_rep

yrep_m2 <- rstan::extract(fit_m2, pars = "y_rep", permuted = TRUE)$y_rep

if (is.null(yrep_m1) || nrow(yrep_m1) == 0) stop("Missing y_rep m1")

if (is.null(yrep_m2) || nrow(yrep_m2) == 0) stop("Missing y_rep m2")

yrep_m1_subset <- yrep_m1[sample.int(nrow(yrep_m1), min(200, nrow(yrep_m1))), , drop = FALSE]

yrep_m2_subset <- yrep_m2[sample.int(nrow(yrep_m2), min(200, nrow(yrep_m2))), , drop = FALSE]

ppc_m1 <- bayesplot::ppc_bars(y = y_m1, yrep = yrep_m1_subset) + ggplot2::ggtitle("Model 1 PPC")

ppc_m2 <- bayesplot::ppc_bars(y = y_m2, yrep = yrep_m2_subset) + ggplot2::ggtitle("Model 2 PPC")

ppc_m1

category_props <- function(y_vec) {

  tab <- table(factor(y_vec, levels = c(1, 2, 3)))

  as.numeric(tab) / sum(tab)

}

obs_prop_m1 <- category_props(y_m1)

obs_prop_m2 <- category_props(y_m2)

pred_prop_m1 <- t(apply(yrep_m1, 1, category_props))

pred_prop_m2 <- t(apply(yrep_m2, 1, category_props))

ppc_table_m1 <- tibble::tibble(

  category = c(1, 2, 3),

  observed = obs_prop_m1,

  pred_mean = apply(pred_prop_m1, 2, mean),

  pred_q025 = apply(pred_prop_m1, 2, stats::quantile, probs = 0.025),

  pred_q975 = apply(pred_prop_m1, 2, stats::quantile, probs = 0.975)

) %>%

  dplyr::mutate(model = "Model1")

ppc_table_m2 <- tibble::tibble(

  category = c(1, 2, 3),

  observed = obs_prop_m2,

  pred_mean = apply(pred_prop_m2, 2, mean),

  pred_q025 = apply(pred_prop_m2, 2, stats::quantile, probs = 0.025),

  pred_q975 = apply(pred_prop_m2, 2, stats::quantile, probs = 0.975)

) %>%

  dplyr::mutate(model = "Model2")

ppc_category_table <- dplyr::bind_rows(ppc_table_m1, ppc_table_m2)

ppc_category_table

ideology_levels <- levels(ideology_m2)

ppc_by_ideology <- dplyr::bind_rows(lapply(ideology_levels, function(g) {

  idx <- which(ideology_m2 == g)

  obs <- category_props(y_m2[idx])

  pred_sub <- yrep_m2[, idx, drop = FALSE]

  pred_props <- t(apply(pred_sub, 1, category_props))

  tibble::tibble(

    ideology = g,

    category = c(1, 2, 3),

    observed = obs,

    pred_mean = apply(pred_props, 2, mean),

    pred_q025 = apply(pred_props, 2, stats::quantile, probs = 0.025),

    pred_q975 = apply(pred_props, 2, stats::quantile, probs = 0.975)

  )

}))

ppc_by_ideology

p_loo <- ggplot2::ggplot(loo_summary, ggplot2::aes(x = model, y = elpd_loo)) +

  ggplot2::geom_col() +

  ggplot2::geom_errorbar(ggplot2::aes(ymin = elpd_loo - se_elpd_loo, ymax = elpd_loo + se_elpd_loo), width = 0.2) +

  ggplot2::labs(x = "Model", y = "elpd_loo", title = "PSIS-LOO comparison")

p_ppc_cat <- ggplot2::ggplot(ppc_category_table, ggplot2::aes(x = factor(category), y = observed, color = model, group = model)) +

  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.25), size = 2) +

  ggplot2::geom_point(ggplot2::aes(y = pred_mean), position = ggplot2::position_dodge(width = 0.25), shape = 1, size = 2) +

  ggplot2::labs(x = "Outcome category", y = "Proportion", title = "Observed vs predicted category proportions")

p_ppc_ideo <- ggplot2::ggplot(ppc_by_ideology, ggplot2::aes(x = factor(category), y = observed)) +

  ggplot2::geom_point(size = 2) +

  ggplot2::geom_point(ggplot2::aes(y = pred_mean), shape = 1, size = 2) +

  ggplot2::facet_wrap(~ ideology) +

  ggplot2::labs(x = "Outcome category", y = "Proportion", title = "Model 2 PPC by ideology")

p_loo

dir.create("../output/model_compare", recursive = TRUE, showWarnings = FALSE)

readr::write_csv(loo_summary, "../output/model_compare/loo_summary.csv")

readr::write_csv(loo_diff, "../output/model_compare/loo_diff.csv")

readr::write_csv(ppc_category_table, "../output/model_compare/ppc_category_table.csv")

readr::write_csv(ppc_by_ideology, "../output/model_compare/ppc_by_ideology.csv")

ggplot2::ggsave(filename = "../output/model_compare/loo_comparison.png", plot = p_loo, width = 6.5, height = 4.5, dpi = 150)

ggplot2::ggsave(filename = "../output/model_compare/ppc_model1_bars.png", plot = ppc_m1, width = 6.5, height = 4.5, dpi = 150)

ggplot2::ggsave(filename = "../output/model_compare/ppc_model2_bars.png", plot = ppc_m2, width = 6.5, height = 4.5, dpi = 150)

ggplot2::ggsave(filename = "../output/model_compare/ppc_category_compare.png", plot = p_ppc_cat, width = 7, height = 4.5, dpi = 150)

ggplot2::ggsave(filename = "../output/model_compare/ppc_model2_by_ideology.png", plot = p_ppc_ideo, width = 8, height = 4.5, dpi = 150)

