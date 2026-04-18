set.seed(405)
options(stringsAsFactors = FALSE)

pkgs <- c("rstan", "dplyr", "readr", "tibble", "bayesplot", "ggplot2")

need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(need) > 0) install.packages(need, repos = "https://cloud.r-project.org")
for (p in pkgs) library(p, character.only = TRUE)
rstan_options(auto_write = TRUE)
safe_cores <- if (tolower(Sys.info()[["sysname"]]) == "darwin") 1L else max(1L, parallel::detectCores(logical = FALSE) - 1L)
options(mc.cores = safe_cores)

d <- readr::read_csv("../output/derived/gss_2024_model_table.csv", show_col_types = FALSE) %>%

  dplyr::mutate(
    degree = as.factor(degree),
    sex = as.factor(sex),
    region = as.factor(region)
  )

df <- d %>%
  dplyr::filter(
    !is.na(confinan_ord), !is.na(polviews), !is.na(age_std), !is.na(income_std),
    !is.na(degree), !is.na(sex), !is.na(region)
  )

f <- stats::as.formula(~ polviews + age_std + income_std + degree + sex + region)
X <- stats::model.matrix(f, data = df)
X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
y <- as.integer(df$confinan_ord)

if (!all(y %in% c(1L, 2L, 3L))) stop("Bad y values")

if (length(y) != nrow(X)) stop("y and X mismatch")

stan_data <- list(
  N = nrow(X),
  K = ncol(X),
  y = y,
  X = X
)

tibble::tibble(metric = c("N", "K"), value = c(stan_data$N, stan_data$K))

stan_code <- "

data {

  int<lower=1> N;

  int<lower=1> K;

  int<lower=1,upper=3> y[N];

  matrix[N, K] X;

}

parameters {

  vector[K] beta;

  ordered[2] c;

}

model {

  beta ~ normal(0, 1);

  c ~ normal(0, 2.5);

  y ~ ordered_logistic(X * beta, c);

}

generated quantities {

  vector[N] log_lik;

  int<lower=1,upper=3> y_rep[N];

  for (n in 1:N) {

    log_lik[n] = ordered_logistic_lpmf(y[n] | (X[n] * beta), c);

    y_rep[n] = ordered_logistic_rng((X[n] * beta), c);

  }

}

"

dir.create("../stan", recursive = TRUE, showWarnings = FALSE)

writeLines(stan_code, con = "../stan/model1_pooled_ordinal.stan")

if (exists("fit_m1_nuts")) rm(fit_m1_nuts)

gc()

sm_m1 <- rstan::stan_model(file = "../stan/model1_pooled_ordinal.stan", model_name = "model1_pooled_ordinal")

fit_m1_nuts <- rstan::sampling(
  object = sm_m1,
  data = stan_data,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 405,
  refresh = 100,
  control = list(adapt_delta = 0.9)
 )

lp <- rstan::extract(fit_m1_nuts, pars = "lp__", permuted = TRUE)$lp__
if (length(lp) == 0) stop("No draws")
length(lp)

sum_raw <- as.data.frame(rstan::summary(fit_m1_nuts)$summary, check.names = FALSE)

if (nrow(sum_raw) == 0) stop("Summary empty")

sum_raw$parameter <- rownames(sum_raw)

fit_summary <- sum_raw %>% dplyr::relocate(parameter)

show_cols <- c("parameter")
if ("mean" %in% colnames(fit_summary)) show_cols <- c(show_cols, "mean")
if ("Mean" %in% colnames(fit_summary)) show_cols <- c(show_cols, "Mean")
if ("sd" %in% colnames(fit_summary)) show_cols <- c(show_cols, "sd")
if ("SD" %in% colnames(fit_summary)) show_cols <- c(show_cols, "SD")
if ("2.5%" %in% colnames(fit_summary)) show_cols <- c(show_cols, "2.5%")
if ("50%" %in% colnames(fit_summary)) show_cols <- c(show_cols, "50%")
if ("97.5%" %in% colnames(fit_summary)) show_cols <- c(show_cols, "97.5%")
if ("n_eff" %in% colnames(fit_summary)) show_cols <- c(show_cols, "n_eff")
if ("Bulk_ESS" %in% colnames(fit_summary)) show_cols <- c(show_cols, "Bulk_ESS")
if ("Rhat" %in% colnames(fit_summary)) show_cols <- c(show_cols, "Rhat")

fit_summary %>%
  dplyr::select(dplyr::all_of(unique(show_cols))) %>%
  head(20)

colnames(fit_summary)

ess_col <- if ("n_eff" %in% colnames(fit_summary)) "n_eff" else if ("Bulk_ESS" %in% colnames(fit_summary)) "Bulk_ESS" else NA_character_
rhat_col <- if ("Rhat" %in% colnames(fit_summary)) "Rhat" else NA_character_

diag_summary <- tibble::tibble(

  metric = c("max_rhat", "min_ess", "n_params_rhat_gt_1.01"),

  value = c(
    if (!is.na(rhat_col)) max(fit_summary[[rhat_col]], na.rm = TRUE) else NA_real_,
    if (!is.na(ess_col)) min(fit_summary[[ess_col]], na.rm = TRUE) else NA_real_,
    if (!is.na(rhat_col)) sum(fit_summary[[rhat_col]] > 1.01, na.rm = TRUE) else NA_real_
  )
)

diag_summary

yrep <- rstan::extract(fit_m1_nuts, pars = "y_rep", permuted = TRUE)$y_rep

if (is.null(yrep) || nrow(yrep) == 0) stop("No y_rep")

yrep_subset <- yrep[sample.int(nrow(yrep), size = min(200, nrow(yrep))), , drop = FALSE]

ppc_plot <- bayesplot::ppc_bars(y = y, yrep = yrep_subset)
ppc_plot

dir.create("../output/model1_nuts", recursive = TRUE, showWarnings = FALSE)
saveRDS(fit_m1_nuts, file = "../output/model1_nuts/fit_m1_nuts.rds")
readr::write_csv(fit_summary, "../output/model1_nuts/fit_summary.csv")
readr::write_csv(diag_summary, "../output/model1_nuts/diagnostics_summary.csv")
ggplot2::ggsave(filename = "../output/model1_nuts/ppc_bars_model1.png", plot = ppc_plot, width = 7, height = 4, dpi = 150)

