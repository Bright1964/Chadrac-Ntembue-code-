#!/usr/bin/env Rscript
# =============================================================================
# DRC malaria prevalence simulations (2008–2024)
#
# - Province-level simulations using malariasimulation + malariaEquilibrium
# - Two scenarios:
#     1) With resistance: time-varying dn0, rn0, gamman from posterior draws
#     2) Without resistance: fixed dn0, rn, gamman
# - ITN usage (coverage) from posterior samples per province/year
# - ITN campaign years from distribution calendar (wide -> long)
# - Produces:
#     (i) Top 10 provinces (faceted)
#     (ii) Four target provinces (faceted)
#     (iii) National mean (population-weighted)
# - Overlays DHS prevalence points (2013–2014 and 2023–2024) with 95% CI
#
# Outputs (PNG):
#   outputs/prev_top10_with_DHS_points.png
#   outputs/prev_four_with_DHS_points.png
#   outputs/prev_nat_with_DHS_points.png
# =============================================================================

suppressPackageStartupMessages({
  library(malariasimulation)
  library(malariaEquilibrium)
  library(ggplot2)
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(future)
  library(furrr)
  library(stringi)
  library(stringr)
  library(janitor)
  library(purrr)
  library(tibble)
})

# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------
YEAR_MIN <- 2008L
YEAR_MAX <- 2024L

years            <- 25L
year_in_days     <- 365L
sim_length       <- years * year_in_days
human_population <- 10000L

# Number of posterior draws to run (increase for final runs)
n_draws_use <- 20L

# Parallel
workers <- min(8L, max(1L, parallel::detectCores() - 1L))
future::plan(future::multisession, workers = workers)
options(future.rng.onMisuse = "ignore")

# Output folder
out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("Running with ", workers, " parallel workers; n_draws_use = ", n_draws_use)

# -----------------------------------------------------------------------------
# PATHS (EDIT IF NEEDED)
# -----------------------------------------------------------------------------
paths <- list(
  efficacy_draw_file = "D:/fichier csv/itn_efficacy_posterior_draws_100samples.csv",
  pbo_file           = "D:/fichier csv/itn_efficacy_pbo_posterior_draws_100samples.csv",
  gamma_draw_file    = "D:/fichier csv/gamma_p_first100_samples.csv",
  usage_file         = "D:/fichier csv/Fully_Corrected_ITN_Usage_Coverage_Data.csv",
  dist_file          = "D:/Exel file/bednet_distribution_by_province_2008_2022.xlsx",
  excel_file_user    = "D:/Exel file/DRC_Data.xlsx",
  excel_file_fallback = "/mnt/data/DRC_Data.xlsx",
  
  # DHS prevalence point files
  dhs_prev_2324_csv_user  = "D:/Desktop/DHS_data_2023_2024/DHS_2023_2024_Malaria_Prevalence_Provincial_And_National.csv",
  dhs_prev_1314_xlsx_user = "D:/Exel file/2013-2014_prevalence_DHS.xlsx",
  
  # fallbacks (if uploaded to this chat)
  dhs_prev_2324_csv_fallback  = "/mnt/data/DHS_2023_2024_Malaria_Prevalence_Provincial_And_National.csv",
  dhs_prev_1314_xlsx_fallback = "/mnt/data/2013-2014_prevalence_DHS.xlsx"
)

excel_file <- if (file.exists(paths$excel_file_user)) paths$excel_file_user else paths$excel_file_fallback

required_files <- c(
  paths$efficacy_draw_file, paths$pbo_file, paths$gamma_draw_file,
  paths$usage_file, paths$dist_file, excel_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n- ", paste(missing_files, collapse = "\n- "))
}

# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------
norm_key <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- tolower(x)
  gsub("[^a-z0-9]", "", x)
}

qL <- function(x) as.numeric(quantile(x, 0.025, na.rm = TRUE, names = FALSE))
qU <- function(x) as.numeric(quantile(x, 0.975, na.rm = TRUE, names = FALSE))

extract_year <- function(x) {
  y <- stringr::str_extract(x, "\\d{4}")
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
  yr_cols <- setdiff(names(df), c("province", "provkey"))
  
  yr_map <- tibble(
    col  = yr_cols,
    year = vapply(yr_cols, extract_year, integer(1))
  ) %>% filter(!is.na(year))
  
  if (nrow(yr_map) == 0) {
    warning("No 4-digit year columns detected in distribution calendar.")
    return(tibble(province = character(0), year = integer(0), flag = integer(0)))
  }
  
  df %>%
    select(province, all_of(yr_map$col)) %>%
    pivot_longer(-province, names_to = "col", values_to = "raw_flag") %>%
    left_join(yr_map, by = "col") %>%
    mutate(flag = coerce_flags01(raw_flag)) %>%
    select(province, year, flag)
}

# Robust wrapper for run_simulation() across package versions
run_sim <- function(sim_length, params) {
  fmls <- names(formals(malariasimulation::run_simulation))
  if (all(c("timesteps", "parameters") %in% fmls)) {
    return(malariasimulation::run_simulation(timesteps = sim_length, parameters = params))
  }
  if (all(c("timesteps", "params") %in% fmls)) {
    return(malariasimulation::run_simulation(timesteps = sim_length, params = params))
  }
  if (all(c("tmax", "parameters") %in% fmls)) {
    return(malariasimulation::run_simulation(tmax = sim_length, parameters = params))
  }
  malariasimulation::run_simulation(sim_length, params)
}

# DHS loader
load_dhs_prev <- function(path, point_year, source_label, norm_key_fun) {
  df <- if (grepl("\\.xlsx$", tolower(path))) {
    readxl::read_excel(path)
  } else {
    readr::read_csv(path, show_col_types = FALSE, locale = readr::locale(encoding = "latin1"))
  }
  
  df <- df %>% janitor::clean_names()
  required <- c("province", "prevalence_pct", "lower_ci_pct", "upper_ci_pct")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("DHS file missing columns: ", paste(missing, collapse = ", "), " in ", path)
  }
  
  df %>%
    transmute(
      Region  = as.character(province),
      ProvKey = norm_key_fun(province),
      Year    = as.numeric(point_year),
      Prev    = as.numeric(prevalence_pct),
      Lower   = as.numeric(lower_ci_pct),
      Upper   = as.numeric(upper_ci_pct),
      Source  = source_label
    ) %>%
    filter(is.finite(Year), is.finite(Prev))
}

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------
message("Loading inputs...")

efficacy_draws <- read.csv(paths$efficacy_draw_file, stringsAsFactors = FALSE)
pbo_draws      <- read.csv(paths$pbo_file,           stringsAsFactors = FALSE)
gamma_draws    <- read.csv(paths$gamma_draw_file,    stringsAsFactors = FALSE)
usage_draws    <- read.csv(paths$usage_file,         stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1")

interventions_data <- read_excel(excel_file, sheet = "Interventions") %>% clean_names()
eir_data           <- read_excel(excel_file, sheet = "EIR")           %>% clean_names()
pop_raw            <- read_excel(excel_file, sheet = "Population")    %>% clean_names()

dist_calendar_raw  <- read_excel(paths$dist_file) %>% clean_names()
if (!"province" %in% names(dist_calendar_raw)) names(dist_calendar_raw)[1] <- "province"

# Restrict to rural if present
if ("urban_rural" %in% names(interventions_data)) {
  interventions_data <- interventions_data %>% filter(urban_rural == "rural")
}
if ("urban_rural" %in% names(pop_raw) && any(pop_raw$urban_rural == "rural", na.rm = TRUE)) {
  pop_raw <- pop_raw %>% filter(urban_rural == "rural")
}

# -----------------------------------------------------------------------------
# NORMALISE KEYS + CHECKS
# -----------------------------------------------------------------------------
if (!"name_1" %in% names(pop_raw)) stop("Population sheet must contain column 'name_1'.")
if (!"year"   %in% names(pop_raw)) stop("Population sheet must contain column 'year'.")
if (!"pop"    %in% names(pop_raw)) stop("Population sheet must contain column 'pop'.")

name_dict <- pop_raw %>%
  transmute(ProvKey = norm_key(name_1), Region = as.character(name_1)) %>%
  distinct()

pop_df <- pop_raw %>%
  transmute(
    ProvKey = norm_key(name_1),
    Year    = as.integer(year),
    Pop     = as.numeric(pop)
  ) %>%
  filter(!is.na(ProvKey), !is.na(Year), is.finite(Pop)) %>%
  group_by(ProvKey, Year) %>%
  summarise(Pop = sum(Pop, na.rm = TRUE), .groups = "drop")

if (!"name_1" %in% names(interventions_data)) stop("Interventions sheet must contain column 'name_1'.")
if (!"year"   %in% names(interventions_data)) stop("Interventions sheet must contain column 'year'.")
if (!"eir_drc" %in% names(eir_data))          stop("EIR sheet must contain column 'eir_drc'.")

interventions_data <- interventions_data %>% mutate(ProvKey = norm_key(name_1), year = as.integer(year))
eir_data           <- eir_data           %>% mutate(ProvKey = norm_key(name_1))

# usage file columns
if (!all(c("province","year","sample_id","p") %in% names(usage_draws))) {
  stop("usage_file must contain columns: province, year, sample_id, p")
}
usage_draws <- usage_draws %>% mutate(ProvKey = norm_key(province), year = as.integer(year))

dist_calendar_raw <- dist_calendar_raw %>% mutate(ProvKey = norm_key(province))
campaign_tbl <- make_campaign_tbl(dist_calendar_raw) %>%
  mutate(ProvKey = norm_key(province)) %>%
  filter(year >= YEAR_MIN, year <= YEAR_MAX)

region_keys_to_sim <- sort(unique(campaign_tbl$ProvKey))

get_region_label <- function(pk) {
  out <- name_dict %>% filter(ProvKey == pk) %>% pull(Region)
  if (length(out) == 0 || is.na(out[1])) pk else out[1]
}

# draws
if (!"draw_id" %in% names(efficacy_draws)) stop("efficacy_draws must contain column 'draw_id'.")
if (!"draw_id" %in% names(pbo_draws))      stop("pbo_draws must contain column 'draw_id'.")
if (!"draw_id" %in% names(gamma_draws))    stop("gamma_draws must contain column 'draw_id'.")

draw_ids <- sort(unique(efficacy_draws$draw_id))
draw_ids <- draw_ids[seq_len(min(n_draws_use, length(draw_ids)))]

sample_ids <- sort(unique(usage_draws$sample_id))
if (length(sample_ids) == 0) stop("usage_draws has no sample_id values.")

# -----------------------------------------------------------------------------
# CORE SIMULATION PER PROVINCE
# -----------------------------------------------------------------------------
simulate_one_province <- function(region_key, region_name, draw, usage_sample,
                                  eff_py, eff_pbo, gam_df_orig,
                                  interventions_data, eir_data, campaign_tbl,
                                  year_in_days, sim_length, human_population) {
  
  region_interventions <- interventions_data %>% filter(ProvKey == region_key)
  if (nrow(region_interventions) == 0) return(NULL)
  
  starting_EIR <- eir_data %>% filter(ProvKey == region_key) %>% pull(eir_drc) %>% .[1]
  if (!is.finite(starting_EIR)) return(NULL)
  
  campaign_years <- campaign_tbl %>%
    filter(ProvKey == region_key, flag == 1L) %>%
    pull(year) %>% unique() %>% sort()
  if (length(campaign_years) == 0) return(NULL)
  
  # Efficacy: py-only (<2019) + PBO (>=2019)
  eff_py_  <- eff_py  %>% filter(year < 2019)
  eff_pbo_ <- eff_pbo %>% filter(year >= 2019)
  
  # Extend PBO params: hold 2021 for 2022–2024
  last_pbo <- eff_pbo_ %>% filter(year == 2021)
  if (nrow(last_pbo) > 0) {
    extend_pbo <- last_pbo[rep(1, 3), ]
    extend_pbo$year <- 2022:2024
    eff_pbo_ <- bind_rows(eff_pbo_, extend_pbo)
  }
  
  eff_combined  <- bind_rows(eff_py_, eff_pbo_)
  region_params <- eff_combined %>% filter(year %in% campaign_years)
  
  # Fix 2009 if campaign exists but params start 2010
  if (2009 %in% campaign_years && !(2009 %in% region_params$year)) {
    params_2010 <- eff_combined %>% filter(year == 2010)
    if (nrow(params_2010) > 0) {
      params_2009 <- params_2010
      params_2009$year <- 2009
      region_params <- bind_rows(params_2009, region_params)
    }
  }
  region_params <- arrange(region_params, year)
  
  # gamma_p: extend 2021 to 2022–2024 and fix 2009 similarly
  gam_df <- gam_df_orig
  if (any(c(2022, 2023, 2024) %in% campaign_years)) {
    gam_2021 <- gam_df %>% filter(year == 2021)
    if (nrow(gam_2021) > 0) {
      gam_ext <- gam_2021[rep(1, 3), ]
      gam_ext$year <- 2022:2024
      gam_df <- bind_rows(gam_df, gam_ext)
    }
  }
  if (2009 %in% campaign_years && !(2009 %in% gam_df$year)) {
    gam_2010 <- gam_df %>% filter(year == 2010)
    if (nrow(gam_2010) > 0) {
      gam_2009 <- gam_2010
      gam_2009$year <- 2009
      gam_df <- bind_rows(gam_2009, gam_df)
    }
  }
  gam_df <- arrange(gam_df, year)
  region_gamman <- gam_df %>% filter(year %in% campaign_years)
  
  # ITN usage (coverage): yearly mean for province & sample
  region_usage <- usage_sample %>%
    filter(ProvKey == region_key, year %in% campaign_years) %>%
    group_by(year) %>%
    summarise(coverages = mean(p, na.rm = TRUE), .groups = "drop")
  
  # Fill missing campaign years (last observed; fallback 0.6)
  missing_years <- setdiff(campaign_years, region_usage$year)
  if (length(missing_years) > 0) {
    last_cov <- if (nrow(region_usage) > 0) tail(region_usage$coverages, 1) else 0.6
    region_usage <- bind_rows(
      region_usage,
      tibble(year = missing_years, coverages = rep(last_cov, length(missing_years)))
    )
  }
  region_usage <- arrange(region_usage, year)
  
  # Common valid years
  valid_years <- Reduce(intersect, list(region_usage$year, region_params$year, region_gamman$year))
  if (length(valid_years) == 0) return(NULL)
  
  region_usage  <- region_usage  %>% filter(year %in% valid_years)
  region_params <- region_params %>% filter(year %in% valid_years)
  region_gamman <- region_gamman %>% filter(year %in% valid_years)
  
  start_year <- min(region_interventions$year, na.rm = TRUE)
  bednet_timesteps <- (valid_years - start_year) * year_in_days
  
  # Base equilibrium params
  simparams <- get_parameters(list(human_population = human_population))
  simparams <- set_equilibrium(parameters = simparams, init_EIR = starting_EIR)
  
  # Scenario 1: With resistance (time-varying)
  params_resist <- set_bednets(
    parameters = simparams,
    timesteps  = bednet_timesteps,
    coverages  = matrix(region_usage$coverages, ncol = 1),
    retention  = 19.8 / 12 * year_in_days,
    dn0        = matrix(region_params$dn0, ncol = 1),
    rn         = matrix(region_params$rn0, ncol = 1),
    rnm        = matrix(rep(0.24, length(valid_years)), ncol = 1),
    gamman     = matrix(region_gamman$gamma_p * year_in_days, ncol = 1)
  )
  
  out_resist <- run_sim(sim_length, params_resist)
  prev_resist <- 100 * (out_resist$n_detect_730_3650 / out_resist$n_730_3650)
  
  # Scenario 2: Without resistance (fixed)
  fixed_dn0    <- 0.533
  fixed_rn     <- 0.56
  fixed_gamman <- 2.64 * year_in_days
  
  params_noresist <- set_bednets(
    parameters = simparams,
    timesteps  = bednet_timesteps,
    coverages  = matrix(region_usage$coverages, ncol = 1),
    retention  = 19.8 / 12 * year_in_days,
    dn0        = matrix(rep(fixed_dn0,    length(valid_years)), ncol = 1),
    rn         = matrix(rep(fixed_rn,     length(valid_years)), ncol = 1),
    rnm        = matrix(rep(0.24,         length(valid_years)), ncol = 1),
    gamman     = matrix(rep(fixed_gamman, length(valid_years)), ncol = 1)
  )
  
  out_noresist <- run_sim(sim_length, params_noresist)
  prev_noresist <- 100 * (out_noresist$n_detect_730_3650 / out_noresist$n_730_3650)
  
  time_years <- out_resist$timestep / year_in_days + start_year
  
  df <- bind_rows(
    tibble(Year = time_years, Prevalence = prev_resist,   Scenario = "With resistance",
           Region = region_name, Draw = draw),
    tibble(Year = time_years, Prevalence = prev_noresist, Scenario = "Without resistance",
           Region = region_name, Draw = draw)
  )
  
  lines <- tibble(Year = valid_years, Region = region_name, Draw = draw)
  list(data = df, lines = lines)
}

# -----------------------------------------------------------------------------
# RUN SIMULATIONS (PARALLEL OVER DRAWS)
# -----------------------------------------------------------------------------
message("Simulating...")

results_list <- future_map(
  seq_along(draw_ids),
  function(i) {
    draw <- draw_ids[i]
    sid  <- sample_ids[(i - 1) %% length(sample_ids) + 1]
    
    eff_py  <- efficacy_draws %>% filter(draw_id == draw)
    eff_pbo <- pbo_draws      %>% filter(draw_id == draw)
    gam_df  <- gamma_draws    %>% filter(draw_id == draw)
    
    usage_sample <- usage_draws %>% filter(sample_id == sid)
    
    per_region <- lapply(region_keys_to_sim, function(pk) {
      region_name <- get_region_label(pk)
      simulate_one_province(
        region_key = pk,
        region_name = region_name,
        draw = draw,
        usage_sample = usage_sample,
        eff_py = eff_py,
        eff_pbo = eff_pbo,
        gam_df_orig = gam_df,
        interventions_data = interventions_data,
        eir_data = eir_data,
        campaign_tbl = campaign_tbl,
        year_in_days = year_in_days,
        sim_length = sim_length,
        human_population = human_population
      )
    })
    
    per_region <- Filter(Negate(is.null), per_region)
    if (!length(per_region)) return(NULL)
    
    list(
      data  = bind_rows(lapply(per_region, `[[`, "data")),
      lines = bind_rows(lapply(per_region, `[[`, "lines"))
    )
  },
  .options  = furrr_options(seed = TRUE),
  .progress = TRUE
)

results_list <- Filter(Negate(is.null), results_list)
if (!length(results_list)) stop("No simulation results produced. Check inputs/keys.")

all_data  <- bind_rows(lapply(results_list, `[[`, "data"))
all_lines <- bind_rows(lapply(results_list, `[[`, "lines")) %>%
  filter(Year >= YEAR_MIN, Year <= YEAR_MAX)

# Campaign line subsets
bednet_lines_pyrethroid <- all_lines %>% filter(Year < 2019)
bednet_lines_pbo        <- all_lines %>% filter(Year >= 2019)

# -----------------------------------------------------------------------------
# NATIONAL PREVALENCE (POPULATION-WEIGHTED)
# -----------------------------------------------------------------------------
message("Computing national mean...")

all_prev <- all_data %>%
  mutate(
    ProvKey  = norm_key(Region),
    Year_int = as.integer(floor(Year))
  )

all_prev_pop <- all_prev %>%
  inner_join(pop_df, by = c("ProvKey" = "ProvKey", "Year_int" = "Year")) %>%
  mutate(Prev_abs = (Prevalence / 100) * Pop)

nat_prev_year_draw <- all_prev_pop %>%
  group_by(Year_int, Scenario, Draw) %>%
  summarise(
    Prev_abs_tot = sum(Prev_abs, na.rm = TRUE),
    Pop_tot      = sum(Pop, na.rm = TRUE),
    Prev_nat     = if_else(Pop_tot > 0, 100 * Prev_abs_tot / Pop_tot, NA_real_),
    .groups = "drop"
  ) %>%
  transmute(
    Year       = Year_int,
    Prevalence = Prev_nat,
    Scenario   = Scenario,
    Region     = "DRC (mean)",
    Draw       = Draw
  )

all_prev_plus <- bind_rows(all_data, nat_prev_year_draw)

summary_data <- all_prev_plus %>%
  filter(Year >= YEAR_MIN, Year <= YEAR_MAX) %>%
  group_by(Year, Region, Scenario) %>%
  summarise(
    Median = median(Prevalence, na.rm = TRUE),
    Lower  = qL(Prevalence),
    Upper  = qU(Prevalence),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# DHS POINTS
# -----------------------------------------------------------------------------
dhs_1314_file <- if (file.exists(paths$dhs_prev_1314_xlsx_user)) paths$dhs_prev_1314_xlsx_user else paths$dhs_prev_1314_xlsx_fallback
dhs_2324_file <- if (file.exists(paths$dhs_prev_2324_csv_user))  paths$dhs_prev_2324_csv_user  else paths$dhs_prev_2324_csv_fallback

if (!file.exists(dhs_1314_file)) stop("Missing DHS 2013–2014 prevalence file: ", dhs_1314_file)
if (!file.exists(dhs_2324_file)) stop("Missing DHS 2023–2024 prevalence file: ", dhs_2324_file)

dhs_prev_1314 <- load_dhs_prev(dhs_1314_file, point_year = 2014, source_label = "DHS 2013–2014", norm_key_fun = norm_key)
dhs_prev_2324 <- load_dhs_prev(dhs_2324_file, point_year = 2024, source_label = "DHS 2023–2024", norm_key_fun = norm_key)
dhs_prev_all  <- bind_rows(dhs_prev_1314, dhs_prev_2324)

# Align DHS province labels to plot labels (handle accents/hyphens)
label_lookup <- tibble(
  Region_label = unique(summary_data$Region),
  ProvKey      = norm_key(unique(summary_data$Region))
) %>% distinct(ProvKey, .keep_all = TRUE)

dhs_prev_all <- dhs_prev_all %>%
  left_join(label_lookup, by = "ProvKey") %>%
  mutate(Region = if_else(!is.na(Region_label), Region_label, Region)) %>%
  select(-Region_label)

# -----------------------------------------------------------------------------
# SUBSETS FOR FIGURES
# -----------------------------------------------------------------------------
region_levels <- sort(setdiff(unique(summary_data$Region), "DRC (mean)"))
top10_regions <- head(region_levels, 10)

four_targets      <- c("Kinshasa", "Haut-Katanga", "Kasaï-Central", "Sud-Kivu")
four_targets_keys <- norm_key(four_targets)

summary_top10 <- summary_data %>%
  filter(Region %in% top10_regions) %>%
  mutate(Region = factor(Region, levels = top10_regions))

summary_nat <- summary_data %>% filter(Region == "DRC (mean)")

summary_four <- summary_data %>%
  mutate(RegionKey = norm_key(Region)) %>%
  filter(RegionKey %in% four_targets_keys) %>%
  select(-RegionKey)

# Preserve order of four_targets but keep exact labels present in data
data_names <- unique(summary_four$Region)
map_tab  <- tibble(target_name = four_targets, key = norm_key(four_targets))
data_tab <- tibble(data_name   = data_names,   key = norm_key(data_names))
merge_tab <- inner_join(map_tab, data_tab, by = "key")

four_levels <- merge_tab$data_name[match(four_targets, merge_tab$target_name)]
four_levels <- four_levels[!is.na(four_levels)]

summary_four <- summary_four %>%
  mutate(Region = factor(Region, levels = four_levels))

# campaign lines for each figure
bednet_lines_pyrethroid_top10 <- bednet_lines_pyrethroid %>% filter(Region %in% top10_regions)
bednet_lines_pbo_top10        <- bednet_lines_pbo        %>% filter(Region %in% top10_regions)

bednet_lines_pyrethroid_four  <- bednet_lines_pyrethroid %>% filter(norm_key(Region) %in% norm_key(four_levels))
bednet_lines_pbo_four         <- bednet_lines_pbo        %>% filter(norm_key(Region) %in% norm_key(four_levels))

# DHS subsets
dhs_prev_top10 <- dhs_prev_all %>%
  filter(norm_key(Region) %in% norm_key(top10_regions)) %>%
  mutate(Region = factor(Region, levels = top10_regions))

dhs_prev_four <- dhs_prev_all %>%
  filter(norm_key(Region) %in% norm_key(four_levels)) %>%
  mutate(Region = factor(Region, levels = four_levels))

dhs_prev_nat <- dhs_prev_all %>% filter(Region == "DRC (mean)")

# -----------------------------------------------------------------------------
# PLOTTING
# -----------------------------------------------------------------------------
pal_scen <- c("With resistance" = "#D55E00", "Without resistance" = "#0072B2")

theme_white_all <- theme(
  panel.background      = element_rect(fill = "white", colour = NA),
  plot.background       = element_rect(fill = "white", colour = NA),
  legend.background     = element_rect(fill = "white", colour = NA),
  legend.box.background = element_rect(fill = "white", colour = NA),
  strip.background      = element_rect(fill = "skyblue", colour = "black"),
  strip.text.y.right    = element_text(angle = 90, face = "bold", colour = "black"),
  strip.text.x          = element_text(face = "bold", colour = "black"),
  panel.grid.major      = element_blank(),
  panel.grid.minor      = element_blank(),
  panel.border          = element_rect(fill = NA, colour = "black"),
  axis.text.x           = element_text(angle = 45, hjust = 1),
  legend.position       = "bottom",
  strip.placement       = "outside"
)

label_pyre <- tibble(x = 2012, y = 95, lab = "Pyrethroid-only")
label_pbo  <- tibble(x = 2021, y = 95, lab = "PBO")

plot_facets <- function(df_summary, df_lines_py, df_lines_pbo, df_dhs, ncol_facets, ylab) {
  ggplot(df_summary, aes(x = Year, y = Median, colour = Scenario, fill = Scenario)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.25, colour = NA) +
    geom_line(linewidth = 0.6) +
    geom_vline(data = df_lines_py, aes(xintercept = Year),
               colour = "black", linetype = "dashed", linewidth = 0.3) +
    geom_vline(data = df_lines_pbo, aes(xintercept = Year),
               colour = "darkgreen", linetype = "solid", linewidth = 0.3) +
    geom_text(data = label_pyre, aes(x = x, y = y, label = lab),
              inherit.aes = FALSE, colour = "black", fontface = "bold", size = 4) +
    geom_text(data = label_pbo, aes(x = x, y = y, label = lab),
              inherit.aes = FALSE, colour = "darkgreen", fontface = "bold", size = 4) +
    geom_errorbar(data = df_dhs, aes(x = Year, ymin = Lower, ymax = Upper),
                  inherit.aes = FALSE, width = 0.15, linewidth = 0.5, colour = "black") +
    geom_point(data = df_dhs, aes(x = Year, y = Prev),
               inherit.aes = FALSE, shape = 16, size = 2.4, colour = "black") +
    scale_color_manual(values = pal_scen) +
    scale_fill_manual(values  = pal_scen) +
    scale_x_continuous(limits = c(YEAR_MIN, YEAR_MAX), breaks = seq(YEAR_MIN, YEAR_MAX, by = 3)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
    labs(x = "Year", y = ylab) +
    theme_classic(base_size = 16) + theme_white_all +
    facet_wrap(~ Region, ncol = ncol_facets, strip.position = "top")
}

plot_national <- function(df_summary_nat, df_dhs_nat) {
  ggplot(df_summary_nat, aes(x = Year, y = Median, colour = Scenario, fill = Scenario)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.25, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_text(data = label_pyre, aes(x = x, y = y, label = lab),
              inherit.aes = FALSE, colour = "black", fontface = "bold", size = 4.5) +
    geom_text(data = label_pbo, aes(x = x, y = y, label = lab),
              inherit.aes = FALSE, colour = "darkgreen", fontface = "bold", size = 4.5) +
    geom_errorbar(data = df_dhs_nat, aes(x = Year, ymin = Lower, ymax = Upper),
                  inherit.aes = FALSE, width = 0.15, linewidth = 0.6, colour = "black") +
    geom_point(data = df_dhs_nat, aes(x = Year, y = Prev),
               inherit.aes = FALSE, shape = 16, size = 3.0, colour = "black") +
    scale_color_manual(values = pal_scen) +
    scale_fill_manual(values  = pal_scen) +
    scale_x_continuous(limits = c(YEAR_MIN, YEAR_MAX), breaks = seq(YEAR_MIN, YEAR_MAX, by = 3)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
    labs(x = "Year", y = "DRC prevalence (%)") +
    theme_classic(base_size = 16) + theme_white_all
}

p_prev_top10 <- plot_facets(
  df_summary   = summary_top10,
  df_lines_py  = bednet_lines_pyrethroid_top10,
  df_lines_pbo = bednet_lines_pbo_top10,
  df_dhs       = dhs_prev_top10,
  ncol_facets  = 5,
  ylab         = "Prevalence (%)"
)

p_prev_four <- plot_facets(
  df_summary   = summary_four,
  df_lines_py  = bednet_lines_pyrethroid_four,
  df_lines_pbo = bednet_lines_pbo_four,
  df_dhs       = dhs_prev_four,
  ncol_facets  = 2,
  ylab         = "Prevalence (%)"
)

p_prev_nat <- plot_national(summary_nat, dhs_prev_nat)

print(p_prev_top10)
print(p_prev_four)
print(p_prev_nat)

ggsave(file.path(out_dir, "prev_top10_with_DHS_points.png"), p_prev_top10, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "prev_four_with_DHS_points.png"),  p_prev_four,  width = 11, height = 7, dpi = 300)
ggsave(file.path(out_dir, "prev_nat_with_DHS_points.png"),   p_prev_nat,   width = 10, height = 6, dpi = 300)

message("Done. Figures saved to: ", normalizePath(out_dir))