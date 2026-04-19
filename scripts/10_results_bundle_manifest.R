set.seed(405)

options(stringsAsFactors = FALSE)
pkgs <- c("dplyr", "readr", "tibble", "tidyr", "stringr")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need) > 0) install.packages(need, repos = "https://cloud.r-project.org")
for (p in pkgs) library(p, character.only = TRUE)

diag_m1 <- readr::read_csv("../output/model1_nuts/diagnostics_summary.csv", show_col_types = FALSE) %>%
  tidyr::pivot_wider(names_from = metric, values_from = value) %>%
  dplyr::mutate(model = "Model1")

diag_m2 <- readr::read_csv("../output/model2_nuts/diagnostics_summary.csv", show_col_types = FALSE) %>%
  tidyr::pivot_wider(names_from = metric, values_from = value) %>%
  dplyr::mutate(model = "Model2")

nuts_diag <- dplyr::bind_rows(diag_m1, diag_m2) %>%
  dplyr::select(model, dplyr::everything())

nuts_diag

vi_compare <- readr::read_csv("../output/vi/compare_metrics.csv", show_col_types = FALSE)
loo_summary <- readr::read_csv("../output/model_compare/loo_summary.csv", show_col_types = FALSE)
loo_diff <- readr::read_csv("../output/model_compare/loo_diff.csv", show_col_types = FALSE)
sens_metrics <- readr::read_csv("../output/sensitivity/sensitivity_metrics.csv", show_col_types = FALSE)
recovery_metrics <- readr::read_csv("../output/validation/recovery_metrics.csv", show_col_types = FALSE)

vi_compare

f1 <- nuts_diag %>%
  dplyr::transmute(section = "NUTS diagnostics", key = model, value = paste0("max_rhat=", round(max_rhat, 4), ", min_ess=", round(min_ess, 1), ", divergences=", ifelse(is.na(n_divergences), 0, n_divergences)))

f2 <- vi_compare %>%
  dplyr::transmute(section = "VI vs NUTS", key = model, value = paste0("median_sd_ratio=", round(median_sd_ratio_vi_to_nuts, 3), ", median_abs_mean_diff=", round(median_abs_mean_diff, 3)))

f3 <- loo_summary %>%
  dplyr::transmute(section = "Model comparison", key = model, value = paste0("elpd_loo=", round(elpd_loo, 2), ", se=", round(se_elpd_loo, 2)))

f4 <- loo_diff %>%
  dplyr::transmute(section = "Model comparison", key = metric, value = round(value, 3) %>% as.character())

f5 <- sens_metrics %>%
  dplyr::transmute(section = "Prior sensitivity", key = paste(model, prior, sep = "_"), value = paste0("median_abs_mean_shift=", round(median_abs_mean_shift, 3), ", median_sd_ratio=", round(median_sd_ratio_vs_base, 3)))

f6 <- recovery_metrics %>%
  dplyr::transmute(section = "Synthetic recovery", key = model, value = paste0("coverage_rate=", round(coverage_rate, 3), ", median_abs_error=", round(median_abs_error, 3), ", max_rhat=", round(max_rhat, 3)))

final_summary <- dplyr::bind_rows(
  f1,
  f2,
  f3,
  f4,
  f5,
  f6
)

final_summary

artifact_manifest <- data.frame(
  path = c(
    "output/model1_nuts/diagnostics_summary.csv",
    "output/model2_nuts/diagnostics_summary.csv",
    "output/vi/compare_metrics.csv",
    "output/model_compare/loo_summary.csv",
    "output/model_compare/loo_diff.csv",
    "output/model_compare/ppc_model1_bars.png",
    "output/model_compare/ppc_model2_bars.png",
    "output/sensitivity/sensitivity_metrics.csv",
    "output/validation/recovery_metrics.csv",
    "output/validation/recovery_scatter.png"
  ),
  stage = c(
    "Model 1 NUTS",
    "Model 2 NUTS",
    "VI comparison",
    "Model comparison",
    "Model comparison",
    "Model comparison",
    "Model comparison",
    "Prior sensitivity",
    "Validation",
    "Validation"
  ),
  description = c(
    "Convergence and ESS diagnostics",
    "Convergence, ESS, divergences, treedepth",
    "Aggregate VI vs NUTS discrepancy metrics",
    "PSIS-LOO estimates for both models",
    "elpd difference and uncertainty",
    "PPC bars for Model 1",
    "PPC bars for Model 2",
    "Sensitivity metrics by model and prior",
    "Synthetic recovery metrics",
    "True vs posterior mean recovery plot"
  ),
  stringsAsFactors = FALSE
) %>%

  dplyr::mutate(exists = file.exists(file.path("..", path)))

artifact_manifest

dir.create("../output/final", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(final_summary, "../output/final/final_summary_table.csv")
readr::write_csv(artifact_manifest, "../output/final/artifact_manifest.csv")
