
# =============================================================================
# ITN usage trajectories (posterior draws) + DHS points
# - 4 provinces (facet 2x2)
# - National summary (weighted)
#
# Author: Chadrac Ntembue Tshiakuisha
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readxl)
  library(readr)
  library(janitor)
  library(purrr)
  library(ggplot2)
  library(future)
  library(furrr)
  library(scales)
})

# -----------------------------------------------------------------------------
# USER SETTINGS (EDIT PATHS)
# -----------------------------------------------------------------------------
paths <- list(
  dist_file       = "D:/Exel file/bednet_distribution_by_province_2008_2022.xlsx",
  usage_draws_rds  = "D:/file rds/usage_draws_df.rds",
  excel_file       = "D:/Exel file/DRC_Data.xlsx",
  obs2324_path     = "D:/Desktop/DHS_data_2023_2024/itn_usage_by_province_with_counts_and_time.csv",
  obs1314_csv      = "D:/fichier csv/dhs_data_survey_2013_2014.csv",
  obs1314_xlsx     = "D:/Desktop/DHS_data_2013_2014/dhs_data_survey_2013_2014.xlsx"
)

stopifnot(
  file.exists(paths$dist_file),
  file.exists(paths$usage_draws_rds),
  file.exists(paths$excel_file),
  file.exists(paths$obs2324_path)
)

# Provinces to show in 2x2 facet
selected_provinces <- c("Haut-Katanga", "Kasaï-Central", "Kinshasa", "Sud-Kivu")

# Posterior draws to sample (for visualisation only)
n_draws_to_plot <- 100
set.seed(123)

# Parallel
workers <- max(1, parallel::detectCores() - 1)
future::plan(future::multisession, workers = workers)

# Plot window
x_min <- 2009
x_max <- 2025
x_breaks <- seq(x_min, x_max, by = 3)

# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------
date_to_CMC <- function(year, month) ((year - 1900) * 12) + month

CMC_to_decimal_year <- function(cmc) {
  yr <- floor((cmc - 1) / 12) + 1900
  mn <- cmc - ((yr - 1900) * 12)
  yr + (mn - 1) / 12
}

safe_wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

pick_col <- function(df, candidates, required = TRUE) {
  nm <- names(df)
  for (cand in candidates) {
    idx <- which(tolower(nm) == tolower(cand))
    if (length(idx)) return(nm[idx[1]])
  }
  if (required) stop("Expected column not found. Tried: ", paste(candidates, collapse = ", "))
  NA_character_
}

extract_year <- function(x) {
  y <- str_extract(x, "\\d{4}")
  suppressWarnings(as.integer(y))
}

coerce_flags01 <- function(v) {
  if (is.logical(v)) return(as.integer(v))
  if (is.numeric(v)) {
    out <- ifelse(is.finite(v) & v > 0, 1L, 0L)
    out[is.na(out)] <- 0L
    return(out)
  }
  vv <- tolower(trimws(as.character(v)))
  ones <- c("1","true","t","yes","y","x","✓")
  out <- ifelse(vv %in% ones, 1L, suppressWarnings(as.integer(vv)))
  out[is.na(out)] <- 0L
  out
}

make_campaign_tbl <- function(df) {
  stopifnot("province" %in% names(df))
  yr_cols <- setdiff(names(df), "province")
  
  yr_map <- tibble(
    col  = yr_cols,
    year = vapply(yr_cols, extract_year, integer(1))
  ) %>% filter(!is.na(year))
  
  if (nrow(yr_map) == 0) {
    warning("No 4-digit year columns detected in distribution sheet.")
    return(tibble(province = character(0), year = integer(0), flag = integer(0)))
  }
  
  df %>%
    select(province, all_of(yr_map$col)) %>%
    pivot_longer(-province, names_to = "col", values_to = "raw_flag") %>%
    left_join(yr_map, by = "col") %>%
    mutate(flag = coerce_flags01(raw_flag)) %>%
    select(province, year, flag)
}

normalize_prov <- function(x) {
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x <- tolower(gsub("[^a-z]", "", x))
  x
}

# campaign months by year (fallback = June if missing)
campaign_month_lookup <- setNames(
  c(6,6,7,6,7,9,6,8,8,6,6,7,8,7,8,6),
  as.character(2009:2024)
)

compute_for_province <- function(prov, draws_df, campaign_tbl) {
  yrs <- campaign_tbl %>%
    filter(province == prov, flag == 1L, !is.na(year)) %>%
    pull(year) %>%
    unique() %>%
    sort()
  
  if (length(yrs) == 0) {
    draws_df %>%
      mutate(
        province = prov,
        m = Inf,
        C = 0,
        P = 0
      ) -> out
    
    vlines <- tibble(province = character(0), year = integer(0))
    return(list(df = out, v = vlines))
  }
  
  months <- vapply(as.character(yrs), function(y) {
    if (y %in% names(campaign_month_lookup)) campaign_month_lookup[[y]] else 6
  }, numeric(1))
  
  campaign_cmc <- mapply(date_to_CMC, yrs, months)
  
  idx <- findInterval(draws_df$cmc, campaign_cmc)
  no_prior <- idx == 0
  idx[no_prior] <- NA_integer_
  
  last_cmc <- campaign_cmc[idx]
  m <- draws_df$cmc - last_cmc
  m[no_prior] <- 240
  
  C <- draws_df$c0 * exp(-m / draws_df$invlam)
  P <- pmin(C + draws_df$d, 1)
  
  out <- draws_df %>%
    mutate(province = prov, m = m, C = C, P = P)
  
  vlines <- tibble(province = prov, year = yrs)
  
  list(df = out, v = vlines)
}

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------
# Distribution flags
dist_df <- read_excel(paths$dist_file) %>% clean_names()
if (!"province" %in% names(dist_df)) {
  prov_col <- pick_col(dist_df, c("province","name_1","adm1_fr","region"))
  dist_df <- dist_df %>% rename(province = !!prov_col)
}
dist_df <- dist_df %>% mutate(province = str_squish(str_to_title(province)))
campaign_tbl <- make_campaign_tbl(dist_df)

# Usage draws
usage_draws <- readRDS(paths$usage_draws_rds) %>% clean_names()
req <- c("sample_id", "cmc", "c0", "d", "invlam")
if (!all(req %in% names(usage_draws))) {
  stop("usage_draws_df must contain: ", paste(req, collapse = ", "))
}

draw_ids <- sample(unique(usage_draws$sample_id),
                   size = min(n_draws_to_plot, length(unique(usage_draws$sample_id))))
draws_sample <- usage_draws %>%
  filter(sample_id %in% draw_ids) %>%
  mutate(decimal_year = CMC_to_decimal_year(cmc),
         year = floor(decimal_year))

# DHS 2023–2024
obs2324_raw <- read_csv(paths$obs2324_path, show_col_types = FALSE) %>% clean_names()
prov_2324 <- pick_col(obs2324_raw, c("province","adm1_fr","region"))
time_2324 <- pick_col(obs2324_raw, c("survey_mid_decimal_year","decimal_year","mid_decimal_year"))
use_2324  <- pick_col(obs2324_raw, c("usage_pct","overall_use_pct","use_pct"))
lcl_2324  <- pick_col(obs2324_raw, c("lcl","lower_ci","lower"), required = FALSE)
ucl_2324  <- pick_col(obs2324_raw, c("ucl","upper_ci","upper"), required = FALSE)

obs2324 <- obs2324_raw %>%
  transmute(
    province = str_squish(str_to_title(.data[[prov_2324]])),
    decimal_year = as.numeric(.data[[time_2324]]),
    overall_use  = pmin(pmax(as.numeric(.data[[use_2324]]) / 100, 0), 1),
    lcl_prop     = if (!is.na(lcl_2324)) pmin(pmax(as.numeric(.data[[lcl_2324]]) / 100, 0), 1) else NA_real_,
    ucl_prop     = if (!is.na(ucl_2324)) pmin(pmax(as.numeric(.data[[ucl_2324]]) / 100, 0), 1) else NA_real_,
    source = "DHS 2023–2024"
  ) %>%
  filter(is.finite(decimal_year), is.finite(overall_use)) %>%
  group_by(province) %>%
  slice_max(decimal_year, n = 1, with_ties = FALSE) %>%
  ungroup()

# DHS 2013–2014
if (file.exists(paths$obs1314_csv)) {
  obs1314_raw <- read_csv(paths$obs1314_csv, show_col_types = FALSE) %>% clean_names()
} else if (file.exists(paths$obs1314_xlsx)) {
  obs1314_raw <- read_excel(paths$obs1314_xlsx) %>% clean_names()
} else {
  stop("No DHS 2013–2014 file found. Provide CSV or XLSX.")
}

prov_1314 <- pick_col(obs1314_raw, c("province","adm1_fr","region"))
time_1314 <- pick_col(obs1314_raw, c("survey_mid_decimal_year","decimal_year","mid_decimal_year"))
use_1314  <- pick_col(obs1314_raw, c("usage_pct","overall_use_pct","use_pct"))
lcl_1314  <- pick_col(obs1314_raw, c("lcl","lower_ci","lower"), required = FALSE)
ucl_1314  <- pick_col(obs1314_raw, c("ucl","upper_ci","upper"), required = FALSE)

obs1314 <- obs1314_raw %>%
  transmute(
    province = str_squish(str_to_title(.data[[prov_1314]])),
    decimal_year = as.numeric(.data[[time_1314]]),
    overall_use  = pmin(pmax(as.numeric(.data[[use_1314]]) / 100, 0), 1),
    lcl_prop     = if (!is.na(lcl_1314)) pmin(pmax(as.numeric(.data[[lcl_1314]]) / 100, 0), 1) else NA_real_,
    ucl_prop     = if (!is.na(ucl_1314)) pmin(pmax(as.numeric(.data[[ucl_1314]]) / 100, 0), 1) else NA_real_,
    source = "DHS 2013–2014"
  ) %>%
  filter(is.finite(decimal_year), is.finite(overall_use))

obs_points <- bind_rows(obs1314, obs2324) %>%
  mutate(year = floor(decimal_year))

# Population (for national weighting)
pop_raw <- read_excel(paths$excel_file, sheet = "Population") %>% clean_names()
col_prov <- pick_col(pop_raw, c("name_1","province","adm1_fr","region"))
col_year <- pick_col(pop_raw, c("year","annee"))
col_pop  <- pick_col(pop_raw, c("pop","population","rural_pop"))

population <- pop_raw %>%
  transmute(
    province = str_squish(str_to_title(.data[[col_prov]])),
    year     = as.integer(.data[[col_year]]),
    pop      = as.numeric(.data[[col_pop]])
  ) %>%
  filter(!is.na(province), !is.na(year), is.finite(pop))

# -----------------------------------------------------------------------------
# MATCH PROVINCE NAMES SAFELY (requested 4 provinces)
# -----------------------------------------------------------------------------
all_provs <- dist_df %>% distinct(province) %>% pull(province)
norm_all <- normalize_prov(all_provs)

picked <- vapply(selected_provinces, function(tgt) {
  idx <- which(norm_all == normalize_prov(tgt))
  if (length(idx)) all_provs[idx[1]] else NA_character_
}, character(1))

if (anyNA(picked)) {
  warning("Some provinces were not found in distribution file: ",
          paste(selected_provinces[is.na(picked)], collapse = ", "))
}
selected_provs_found <- na.omit(unique(picked))
if (!length(selected_provs_found)) stop("None of the requested provinces were found.")

# -----------------------------------------------------------------------------
# COMPUTE PROVINCE TRAJECTORIES (4 provinces)
# -----------------------------------------------------------------------------
prov_results <- future_map(
  selected_provs_found,
  ~ compute_for_province(.x, draws_sample, campaign_tbl),
  .progress = TRUE
)

final_df_4  <- bind_rows(map(prov_results, "df"))
vlines_df_4 <- bind_rows(map(prov_results, "v"))

# -----------------------------------------------------------------------------
# NATIONAL SERIES (all provinces with campaigns & population)
# -----------------------------------------------------------------------------
provinces_nat <- intersect(
  campaign_tbl %>% distinct(province) %>% pull(province),
  population  %>% distinct(province) %>% pull(province)
)
if (!length(provinces_nat)) stop("No overlap between campaign_tbl provinces and population provinces.")

nat_results <- future_map(
  provinces_nat,
  ~ compute_for_province(.x, draws_sample, campaign_tbl),
  .progress = TRUE
)

final_df_nat <- bind_rows(map(nat_results, "df")) %>%
  mutate(year = floor(decimal_year))

final_df_nat_pop <- final_df_nat %>%
  inner_join(population, by = c("province","year"))

nat_by_sample <- final_df_nat_pop %>%
  group_by(sample_id, decimal_year, year) %>%
  summarise(P_nat_w = safe_wmean(P, pop), .groups = "drop")

nat_summary <- nat_by_sample %>%
  group_by(decimal_year, year) %>%
  summarise(
    P_med = median(P_nat_w, na.rm = TRUE),
    P_lo  = quantile(P_nat_w, 0.025, na.rm = TRUE),
    P_hi  = quantile(P_nat_w, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# National DHS points (weighted) — keep circles (same as provincial)
obs_points_w <- obs_points %>%
  inner_join(population, by = c("province","year"))

nat_dhs_points <- obs_points_w %>%
  group_by(source) %>%
  summarise(
    decimal_year = safe_wmean(decimal_year, pop),
    overall_use  = safe_wmean(overall_use,  pop),
    lcl_prop     = safe_wmean(lcl_prop,     pop),
    ucl_prop     = safe_wmean(ucl_prop,     pop),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# PLOTS
# -----------------------------------------------------------------------------
pal_src <- c("DHS 2013–2014" = "#0072B2", "DHS 2023–2024" = "#D55E00")

# ---- Plot 4 provinces (no national diamonds; circles only) ----
p_prov <- ggplot(
  final_df_4,
  aes(x = decimal_year, y = P, group = sample_id)
) +
  geom_line(alpha = 0.15, linewidth = 0.4, color = "#66BB6A") +
  facet_wrap(~ province, ncol = 2) +
  geom_vline(
    data = vlines_df_4,
    aes(xintercept = year),
    inherit.aes = FALSE,
    color = "grey65", linetype = "22", alpha = 0.7, linewidth = 0.5
  ) +
  geom_errorbar(
    data = obs_points %>%
      filter(province %in% selected_provs_found,
             is.finite(lcl_prop), is.finite(ucl_prop)),
    aes(x = decimal_year, ymin = lcl_prop, ymax = ucl_prop, color = source),
    inherit.aes = FALSE,
    width = 0.06,
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  geom_point(
    data = obs_points %>% filter(province %in% selected_provs_found),
    aes(x = decimal_year, y = overall_use, color = source),
    inherit.aes = FALSE,
    shape = 16,
    size = 3.2
  ) +
  scale_color_manual(values = pal_src, name = "DHS survey") +
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = x_breaks,
    labels = function(x) sprintf("%d", x),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(x = "Year", y = "ITN usage") +
  theme_classic(base_size = 18) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    
    # ✅ Province label box (strip) in skyblue
    strip.background = element_rect(fill = "skyblue", color = "black", linewidth = 0.8),
    strip.text       = element_text(face = "bold", color = "black"),
    
    # Legend inside plot
    legend.position      = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background    = element_rect(fill = scales::alpha("white", 0.75), color = NA),
    legend.key           = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 3.4)))


# ---- Plot national ----
p_nat <- ggplot(nat_summary, aes(x = decimal_year)) +
  geom_ribbon(aes(ymin = P_lo, ymax = P_hi), fill = "#A5D6A7", alpha = 0.25) +
  geom_line(aes(y = P_med), color = "#1B5E20", linewidth = 1.2) +
  geom_errorbar(
    data = nat_dhs_points %>% filter(is.finite(lcl_prop), is.finite(ucl_prop)),
    aes(x = decimal_year, ymin = lcl_prop, ymax = ucl_prop, color = source),
    width = 0.06, linewidth = 0.7
  ) +
  geom_point(
    data = nat_dhs_points,
    aes(x = decimal_year, y = overall_use, color = source),
    shape = 16, size = 3.6
  ) +
  scale_color_manual(values = pal_src, name = "Survey") +
  scale_x_continuous(limits = c(x_min, x_max), breaks = x_breaks,
                     labels = function(x) sprintf("%d", x),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1),
                     labels = percent_format(accuracy = 1)) +
  labs(x = "Year", y = "DRC ITN use") +
  theme_classic(base_size = 18) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position      = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background    = element_rect(fill = alpha("white", 0.75), color = NA),
    legend.key           = element_blank()
  )

print(p_prov)
print(p_nat)

# Optional: save
# ggsave("fig_itn_usage_4provinces.png", p_prov, width = 11, height = 8, dpi = 300)
# ggsave("fig_itn_usage_national.png",   p_nat,  width = 10, height = 6, dpi = 300)
#----------------------------------------------------------------------

