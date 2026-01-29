#!/usr/bin/env Rscript
# =============================================================================
# Fit pyrethroid resistance over time (DRC provinces) using a Beta-Binomial Stan
#
# Inputs:
#   - dt_new.csv with (at minimum): year(s), tested, resistant, and a province field
#
# Outputs:
#   - ggplot object showing:
#       * predictive uncertainty band (beta-binomial)
#       * parameter (model) credible interval band
#       * median fitted curve
#       * observed points (coloured by province) + boxplots by year
#   - printed summary statistics + posterior parameter summaries
#
# Notes:
#   - Province column is detected robustly (by name heuristics + value overlap).
#   - Column names are normalized to be accent/space/punctuation tolerant.
# =============================================================================

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringi)     # accent/ASCII normalization
  library(extraDistr)  # rbbinom()
  library(scales)      # label_percent()
})

# -----------------------------
# Global options
# -----------------------------
set.seed(42)
rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))

# -----------------------------
# User configuration
# -----------------------------
cfg <- list(
  csv_path = "D:/Documents/Excel files/dt_new.csv",
  chains   = 8,
  iter     = 6000,
  warmup   = 2000,
  seed     = 123,
  n_pred   = 100,     # pseudo-sample size for predictive beta-binomial draws
  fine_by  = 0.05     # time grid resolution (in year-index units)
)

# =============================================================================
# Helpers
# =============================================================================

normalize_colnames <- function(nms) {
  nms2 <- trimws(nms)
  nms2 <- stringi::stri_trans_general(nms2, "Latin-ASCII")
  nms2 <- tolower(nms2)
  nms2 <- gsub("[^a-z0-9]+", "_", nms2)
  nms2
}

normalize_key <- function(x) {
  x <- as.character(x)
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- tolower(x)
  gsub("[^a-z0-9]+", "", x)
}

infer_province_column <- function(df) {
  # Candidate names (already normalized)
  candidate_names <- c(
    "adm1_fr","region","province","name_1","adm1fr","adm1","prov",
    "province_name","adm1_name","adm_1","adm1_francais","adm1francais",
    "adm1label","adm1_label","adm_1_label","provincefrancais","province_fr"
  )
  
  hit <- intersect(candidate_names, names(df))
  if (length(hit) > 0) return(hit[1])
  
  rx <- "(^|_)(adm1|prov|province|region|name_?1|adm_?1)($|_)"
  hit2 <- names(df)[grepl(rx, names(df), ignore.case = TRUE)]
  if (length(hit2) > 0) return(hit2[1])
  
  # Fallback: choose best character column by overlap with known provinces
  known_provs <- c(
    "Equateur","Haut-Katanga","Ituri","Kasaï-Central","Kasaï-Oriental",
    "Kinshasa","Kongo-Central","Kwilu","Maniema","Maï-Ndombe",
    "Nord-Kivu","Nord-Ubangi","Sankuru","Sud-Kivu","Tanganyika","Tshopo"
  )
  known_key <- normalize_key(known_provs)
  
  char_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
  best_col <- NULL
  best_overlap <- -1
  
  for (cn in char_cols) {
    vals <- unique(df[[cn]])
    vals_key <- normalize_key(vals)
    overlap <- sum(vals_key %in% known_key, na.rm = TRUE)
    if (overlap > best_overlap) {
      best_overlap <- overlap
      best_col <- cn
    }
  }
  
  if (!is.null(best_col) && best_overlap >= 5) {
    message(sprintf("Heuristically selected '%s' as province column (overlap=%d).",
                    best_col, best_overlap))
    return(best_col)
  }
  
  stop(
    "Could not infer a province column.\n",
    "Tip: rename your province field to something like 'province' or 'adm1_fr'."
  )
}

make_fine_grid <- function(yearnum, years, by = 0.05) {
  tmin <- min(yearnum, na.rm = TRUE)
  tmax <- max(yearnum, na.rm = TRUE)
  ymin <- min(years, na.rm = TRUE)
  
  fine_df <- data.frame(
    time = seq(tmin - 0.5, tmax + 0.05, by = by)
  )
  fine_df$year <- fine_df$time + ymin - 1
  fine_df
}

summarise_numeric <- function(df) {
  df %>%
    summarise(
      across(
        where(is.numeric),
        list(
          mean = ~mean(.x, na.rm = TRUE),
          sd   = ~sd(.x, na.rm = TRUE),
          min  = ~min(.x, na.rm = TRUE),
          max  = ~max(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    )
}

# =============================================================================
# 1) Load and standardize data
# =============================================================================

if (!file.exists(cfg$csv_path)) {
  stop("CSV not found at: ", cfg$csv_path)
}

res_data <- read_csv(cfg$csv_path, show_col_types = FALSE)
message("Original column names:")
print(names(res_data))

names(res_data) <- normalize_colnames(names(res_data))
message("Normalized column names:")
print(names(res_data))

# Ensure year column named `years`
if (!"years" %in% names(res_data) && "year" %in% names(res_data)) {
  res_data <- res_data %>% rename(years = year)
}
if (!"years" %in% names(res_data)) {
  stop("No 'years' (or 'year') column found after normalization.")
}
res_data <- res_data %>% mutate(years = suppressWarnings(as.numeric(years)))

# Ensure tested/resistant columns
if (!("msq_tested" %in% names(res_data)) && "tested" %in% names(res_data)) {
  res_data <- res_data %>% rename(msq_tested = tested)
}
if (!("msq_resistant" %in% names(res_data)) && "resistant" %in% names(res_data)) {
  res_data <- res_data %>% rename(msq_resistant = resistant)
}

# Compute probresistance if missing
if (!"probresistance" %in% names(res_data)) {
  if (all(c("msq_resistant", "msq_tested") %in% names(res_data))) {
    res_data <- res_data %>%
      mutate(probresistance = ifelse(msq_tested > 0, msq_resistant / msq_tested, NA_real_))
  } else {
    stop("Need 'probresistance' OR both 'msq_resistant' and 'msq_tested'.")
  }
}

# Province column detection + standard fields
prov_col <- infer_province_column(res_data)
message(sprintf("Using '%s' as province column.", prov_col))

res_data <- res_data %>%
  mutate(
    Province_raw = .data[[prov_col]],
    Province     = as.character(Province_raw),
    Province_key = normalize_key(Province_raw)
  )

# Time index for Stan
res_data <- res_data %>%
  mutate(yearnum = years - min(years, na.rm = TRUE) + 1)

# Sanity checks
needed <- c("msq_tested", "msq_resistant", "probresistance", "years", "yearnum", "Province")
missing <- setdiff(needed, names(res_data))
if (length(missing) > 0) stop("Missing columns after processing: ", paste(missing, collapse = ", "))

# Fine grid for smooth curve
fine_df <- make_fine_grid(res_data$yearnum, res_data$years, by = cfg$fine_by)

# Optional summary
message("Numeric summary (quick check):")
print(summarise_numeric(res_data))

# =============================================================================
# 2) Stan model (Beta-Binomial with a flexible sigmoid curve)
# =============================================================================

stan_code <- "
data {
  int<lower = 1> N;
  int<lower = 0> n[N];
  int<lower = 0> r[N];
  vector<lower = 0>[N] t;
  int<lower = 1> T;
  vector<lower = 0>[T] tau;
}
parameters {
  real<lower = 0> a;
  real<lower = 0> b;
  real<lower = 0> c;
  real<lower = 0> k;
}
transformed parameters {
  vector[N] m_t = (1.0 ./ (1.0 + exp(-a * (t - b))) ^ (1.0 / c));
}
model {
  a ~ lognormal(0.75, 0.75);
  b ~ lognormal(1.0, 1.0);
  c ~ lognormal(0.5, 0.5);
  k ~ pareto(1.0, 1.5);

  r ~ beta_binomial(n, k * m_t, k * (1.0 - m_t));
}
generated quantities {
  vector[T] p = (1.0 ./ (1.0 + exp(-a * (tau - b))) ^ (1.0 / c));
}
"

sm <- rstan::stan_model(model_code = stan_code)

stan_data <- list(
  N   = nrow(res_data),
  n   = tidyr::replace_na(res_data$msq_tested, 0L),
  r   = tidyr::replace_na(res_data$msq_resistant, 0L),
  t   = res_data$yearnum,
  T   = nrow(fine_df),
  tau = fine_df$time
)

fit <- rstan::sampling(
  object = sm,
  data   = stan_data,
  chains = cfg$chains,
  iter   = cfg$iter,
  warmup = cfg$warmup,
  seed   = cfg$seed
)

post <- rstan::extract(fit)

# =============================================================================
# 3) Derived summaries (median curve, model CI, predictive CI)
# =============================================================================

# Median parameter curve
median_a <- median(post$a)
median_b <- median(post$b)
median_c <- median(post$c)

median_curve <- 1 / (1 + exp(-median_a * (fine_df$time - median_b)))^(1 / median_c)

# Model CI from generated quantities p
model_ci <- apply(post$p, 2, quantile, probs = c(0.025, 0.975))

# Predictive uncertainty band: simulate beta-binomial outcomes then convert to proportions
n_samples <- length(post$a)
pred_draws <- matrix(NA_real_, nrow = n_samples, ncol = nrow(fine_df))

for (i in seq_len(n_samples)) {
  a_i <- post$a[i]; b_i <- post$b[i]; c_i <- post$c[i]; k_i <- post$k[i]
  m_pred <- 1 / (1 + exp(-a_i * (fine_df$time - b_i)))^(1 / c_i)
  
  alpha <- k_i * m_pred
  beta  <- k_i * (1 - m_pred)
  
  # rbbinom(n, size, alpha, beta): returns counts
  pred_counts <- extraDistr::rbbinom(n = length(m_pred), size = cfg$n_pred, alpha = alpha, beta = beta)
  pred_draws[i, ] <- pred_counts / cfg$n_pred
}

pred_ci <- apply(pred_draws, 2, quantile, probs = c(0.025, 0.975))

# Posterior parameter summaries
cat("\nPosterior means and 95% credible intervals:\n")
for (nm in c("a", "b", "c", "k")) {
  mu <- mean(post[[nm]])
  ci <- quantile(post[[nm]], probs = c(0.025, 0.975))
  cat(sprintf("  %s: %.3f [%.3f, %.3f]\n", nm, mu, ci[1], ci[2]))
}

# =============================================================================
# 4) Plot
# =============================================================================

year_min <- floor(min(res_data$years, na.rm = TRUE))
year_max <- ceiling(max(res_data$years, na.rm = TRUE))

p <- ggplot() +
  # Predictive CI (wider band)
  geom_ribbon(
    data = fine_df,
    aes(x = year, ymin = pred_ci[1, ], ymax = pred_ci[2, ]),
    alpha = 0.30
  ) +
  # Model CI (narrower band)
  geom_ribbon(
    data = fine_df,
    aes(x = year, ymin = model_ci[1, ], ymax = model_ci[2, ]),
    alpha = 0.50
  ) +
  # Boxplots by year
  geom_boxplot(
    data = res_data,
    aes(x = years, y = probresistance, group = factor(years)),
    width = 0.5, alpha = 0.6, outlier.shape = NA
  ) +
  # Median curve
  geom_line(
    data = fine_df,
    aes(x = year, y = median_curve),
    linewidth = 0.9
  ) +
  # Observations (coloured by province)
  geom_point(
    data = res_data,
    aes(x = years, y = probresistance, color = Province),
    size = 2,
    position = position_jitter(width = 0.15, height = 0, seed = 42)
  ) +
  scale_x_continuous(name = "Year", breaks = seq(year_min, year_max, by = 1)) +
  scale_y_continuous(
    name   = "Probability of resistance",
    limits = c(0, 1),
    labels = scales::label_percent(accuracy = 1)
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    plot.background  = element_rect(fill = "white", colour = NA),
    axis.text        = element_text(colour = "black"),
    axis.title       = element_text(colour = "black"),
    legend.position  = "right"
  ) +
  labs(title = NULL, subtitle = NULL)

print(p)

# =============================================================================
# 5) Optional: mean fitted resistance per calendar year (from post$p)
# =============================================================================

mean_p_time <- colMeans(post$p)
fine_df$mean_resistance <- mean_p_time
fine_df$year_rounded <- round(fine_df$year)

mean_res_by_year <- fine_df %>%
  group_by(year_rounded) %>%
  summarise(mean_resistance = mean(mean_resistance), .groups = "drop") %>%
  arrange(year_rounded)

message("\nMean fitted resistance per calendar year (posterior mean):")
print(mean_res_by_year)