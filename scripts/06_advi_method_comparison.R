set.seed(405)
options(stringsAsFactors = FALSE)

pkgs <- c("rstan", "dplyr", "readr", "tibble", "stringr", "tidyr", "ggplot2")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need) > 0) install.packages(need, repos = "https://cloud.r-project.org")
for (p in pkgs) library(p, character.only = TRUE)

rstan_options(auto_write = TRUE)

d <- readr::read_csv("../output/derived/gss_2024_model_table.csv", show_col_types = FALSE) %>%
  dplyr::mutate(
    degree = as.factor(degree),
    sex = as.factor(sex),
    region = as.factor(region),
    polviews_bin = factor(polviews_bin, levels = c("Liberal", "Moderate", "Conservative"))
  )

m1 <- d %>%
  dplyr::filter(
    !is.na(confinan_ord), !is.na(polviews), !is.na(age_std), !is.na(income_std),
    !is.na(degree), !is.na(sex), !is.na(region)
  )

X_m1 <- stats::model.matrix(~ polviews + age_std + income_std + degree + sex + region, data = m1)
X_m1 <- X_m1[, colnames(X_m1) != "(Intercept)", drop = FALSE]
y_m1 <- as.integer(m1$confinan_ord)

stan_data_m1 <- list(
  N = nrow(X_m1),
  K = ncol(X_m1),
  y = y_m1,
  X = X_m1
)

tibble::tibble(metric = c("m1_N", "m1_K"), value = c(stan_data_m1$N, stan_data_m1$K))

m2 <- d %>%
  dplyr::filter(
    !is.na(confinan_ord), !is.na(polviews_bin), !is.na(age_std), !is.na(income_std),
    !is.na(degree), !is.na(sex), !is.na(region)
  )

X_m2 <- stats::model.matrix(~ age_std + income_std + degree + sex + region, data = m2)
X_m2 <- X_m2[, colnames(X_m2) != "(Intercept)", drop = FALSE]
y_m2 <- as.integer(m2$confinan_ord)

ideology_levels <- levels(m2$polviews_bin)
ideology_id <- as.integer(m2$polviews_bin)

stan_data_m2 <- list(
  N = nrow(X_m2),
  K = ncol(X_m2),
  G = length(ideology_levels),
  y = y_m2,
  X = X_m2,
  ideology_id = ideology_id
)

tibble::tibble(metric = c("m2_N", "m2_K", "m2_G"), value = c(stan_data_m2$N, stan_data_m2$K, stan_data_m2$G))
sm_m1 <- rstan::stan_model(file = "../stan/model1_pooled_ordinal.stan", model_name = "model1_pooled_ordinal")
sm_m2 <- rstan::stan_model(file = "../stan/model2_hierarchical_ordinal.stan", model_name = "model2_hierarchical_ordinal")

fit_m1_vi <- rstan::vb(
  object = sm_m1,
  data = stan_data_m1,
  algorithm = "meanfield",
  output_samples = 4000,
  seed = 405,
  tol_rel_obj = 0.001,
  eval_elbo = 100,
  iter = 20000
)

lp_m1 <- rstan::extract(fit_m1_vi, pars = "lp__", permuted = TRUE)$lp__

length(lp_m1)

fit_m2_vi <- rstan::vb(
  object = sm_m2,
  data = stan_data_m2,
  algorithm = "meanfield",
  output_samples = 4000,
  seed = 406,
  tol_rel_obj = 0.001,
  eval_elbo = 100,
  iter = 20000
)

lp_m2 <- rstan::extract(fit_m2_vi, pars = "lp__", permuted = TRUE)$lp__

length(lp_m2)

m1_vi_summary <- as.data.frame(rstan::summary(fit_m1_vi)$summary, check.names = FALSE)
m1_vi_summary$parameter <- rownames(m1_vi_summary)
m1_vi_summary <- m1_vi_summary %>% dplyr::relocate(parameter)

m2_vi_summary <- as.data.frame(rstan::summary(fit_m2_vi)$summary, check.names = FALSE)
m2_vi_summary$parameter <- rownames(m2_vi_summary)
m2_vi_summary <- m2_vi_summary %>% dplyr::relocate(parameter)

m1_vi_key <- m1_vi_summary %>%
  dplyr::filter(stringr::str_detect(parameter, "^beta\\[|^c\\[")) %>%
  dplyr::select(parameter, mean_vi = mean, sd_vi = sd, q025_vi = `2.5%`, q50_vi = `50%`, q975_vi = `97.5%`)

m2_vi_key <- m2_vi_summary %>%
  dplyr::filter(stringr::str_detect(parameter, "^beta\\[|^c\\[|^mu_alpha$|^sigma_alpha$|^alpha\\[")) %>%
  dplyr::select(parameter, mean_vi = mean, sd_vi = sd, q025_vi = `2.5%`, q50_vi = `50%`, q975_vi = `97.5%`)

list(nrow(m1_vi_key), nrow(m2_vi_key))

m1_nuts_summary <- readr::read_csv("../output/model1_nuts/fit_summary.csv", show_col_types = FALSE) %>%
  dplyr::filter(stringr::str_detect(parameter, "^beta\\[|^c\\[")) %>%
  dplyr::select(parameter, mean_nuts = mean, sd_nuts = sd, q025_nuts = `2.5%`, q50_nuts = `50%`, q975_nuts = `97.5%`)

m2_nuts_summary <- readr::read_csv("../output/model2_nuts/fit_summary.csv", show_col_types = FALSE) %>%
  dplyr::filter(stringr::str_detect(parameter, "^beta\\[|^c\\[|^mu_alpha$|^sigma_alpha$|^alpha\\[")) %>%
  dplyr::select(parameter, mean_nuts = mean, sd_nuts = sd, q025_nuts = `2.5%`, q50_nuts = `50%`, q975_nuts = `97.5%`)

m1_compare <- m1_nuts_summary %>%
  dplyr::inner_join(m1_vi_key, by = "parameter") %>%
  dplyr::mutate(
    abs_mean_diff = abs(mean_nuts - mean_vi),
    sd_ratio_vi_to_nuts = sd_vi / sd_nuts
  ) %>%
  dplyr::arrange(dplyr::desc(abs_mean_diff))

m2_compare <- m2_nuts_summary %>%
  dplyr::inner_join(m2_vi_key, by = "parameter") %>%
  dplyr::mutate(
    abs_mean_diff = abs(mean_nuts - mean_vi),
    sd_ratio_vi_to_nuts = sd_vi / sd_nuts
  ) %>%
  dplyr::arrange(dplyr::desc(abs_mean_diff))

list(nrow(m1_compare), nrow(m2_compare))

compare_metrics <- tibble::tibble(
  model = c("Model1", "Model2"),
  n_parameters = c(nrow(m1_compare), nrow(m2_compare)),
  median_abs_mean_diff = c(stats::median(m1_compare$abs_mean_diff), stats::median(m2_compare$abs_mean_diff)),
  median_sd_ratio_vi_to_nuts = c(stats::median(m1_compare$sd_ratio_vi_to_nuts), stats::median(m2_compare$sd_ratio_vi_to_nuts)),
  p10_sd_ratio_vi_to_nuts = c(stats::quantile(m1_compare$sd_ratio_vi_to_nuts, 0.10), stats::quantile(m2_compare$sd_ratio_vi_to_nuts, 0.10)),
  p90_sd_ratio_vi_to_nuts = c(stats::quantile(m1_compare$sd_ratio_vi_to_nuts, 0.90), stats::quantile(m2_compare$sd_ratio_vi_to_nuts, 0.90))
)

compare_metrics

p_m1_mean <- ggplot2::ggplot(m1_compare, ggplot2::aes(x = mean_nuts, y = mean_vi)) +
  ggplot2::geom_point(alpha = 0.7) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::labs(x = "NUTS mean", y = "ADVI mean", title = "Model 1: posterior means")

p_m2_mean <- ggplot2::ggplot(m2_compare, ggplot2::aes(x = mean_nuts, y = mean_vi)) +
  ggplot2::geom_point(alpha = 0.7) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::labs(x = "NUTS mean", y = "ADVI mean", title = "Model 2: posterior means")

p_m1_sd <- ggplot2::ggplot(m1_compare, ggplot2::aes(x = sd_nuts, y = sd_vi)) +
  ggplot2::geom_point(alpha = 0.7) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::labs(x = "NUTS sd", y = "ADVI sd", title = "Model 1: posterior sd")

p_m2_sd <- ggplot2::ggplot(m2_compare, ggplot2::aes(x = sd_nuts, y = sd_vi)) +
  ggplot2::geom_point(alpha = 0.7) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::labs(x = "NUTS sd", y = "ADVI sd", title = "Model 2: posterior sd")

p_m2_alpha <- ggplot2::ggplot(m2_compare %>% dplyr::filter(stringr::str_detect(parameter, "^alpha\\[")), ggplot2::aes(x = parameter, y = sd_ratio_vi_to_nuts)) +
  ggplot2::geom_col() +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
  ggplot2::labs(x = "parameter", y = "sd ratio (VI/NUTS)", title = "Model 2: ideology intercept uncertainty ratio")

p_m1_mean

dir.create("../output/vi", recursive = TRUE, showWarnings = FALSE)
saveRDS(fit_m1_vi, file = "../output/vi/fit_m1_vi.rds")
saveRDS(fit_m2_vi, file = "../output/vi/fit_m2_vi.rds")
readr::write_csv(m1_vi_summary, "../output/vi/m1_vi_summary.csv")
readr::write_csv(m2_vi_summary, "../output/vi/m2_vi_summary.csv")
readr::write_csv(m1_compare, "../output/vi/m1_nuts_vs_vi.csv")
readr::write_csv(m2_compare, "../output/vi/m2_nuts_vs_vi.csv")
readr::write_csv(compare_metrics, "../output/vi/compare_metrics.csv")
ggplot2::ggsave(filename = "../output/vi/m1_mean_compare.png", plot = p_m1_mean, width = 6.5, height = 4.5, dpi = 150)
ggplot2::ggsave(filename = "../output/vi/m2_mean_compare.png", plot = p_m2_mean, width = 6.5, height = 4.5, dpi = 150)
ggplot2::ggsave(filename = "../output/vi/m1_sd_compare.png", plot = p_m1_sd, width = 6.5, height = 4.5, dpi = 150)
ggplot2::ggsave(filename = "../output/vi/m2_sd_compare.png", plot = p_m2_sd, width = 6.5, height = 4.5, dpi = 150)
ggplot2::ggsave(filename = "../output/vi/m2_alpha_sd_ratio.png", plot = p_m2_alpha, width = 7, height = 4.5, dpi = 150)
