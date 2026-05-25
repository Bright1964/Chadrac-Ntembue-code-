# =============================================================================
# DRC malaria ITN simulation: resistance, PBO nets and WHO-calibrated burden
# =============================================================================
#
# Purpose
#   This script runs province-level malaria simulations for the Democratic
#   Republic of the Congo (DRC), summarises national and subnational burden,
#   and estimates:
#     1. ITN impact compared with no nets,
#     2. pyrethroid-resistance-attributable cases and deaths,
#     3. additional cases and deaths averted by pyrethroid-PBO nets.
#
# Important modelling convention
#   The observed intervention scenario with resistance is treated as the
#   status-quo scenario. WHO case/death matching is applied only to
#   resistance-attributable burden outputs. ITN and PBO comparison plots/tables
#   use the raw population-weighted model outputs.
#
# Author: Chadrac Ntembue Tshiakuisha
# Last updated: 2026-05-25
# =============================================================================


# =============================================================================
# 0. Reproducible setup
# =============================================================================

required_packages <- c(
  "malariasimulation",
  "malariaEquilibrium",
  "ggplot2",
  "readxl",
  "readr",
  "dplyr",
  "tidyr",
  "future",
  "furrr",
  "janitor",
  "stringr",
  "purrr",
  "tibble",
  "stringi",
  "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following R packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  invisible(lapply(required_packages, library, character.only = TRUE))
})

set.seed(123)


# =============================================================================
# 1. User configuration
# =============================================================================

# Core simulation settings ------------------------------------------------------
years <- 25L
days_per_year <- 365L
sim_length <- years * days_per_year
human_population <- 10000L
ANALYSIS_YEARS <- 2009:2024
eps <- 1e-9

# Number of posterior draws and parallel workers can be changed without editing
# the code, for example:
#   Sys.setenv(N_DRAWS = "100", N_WORKERS = "8")
n_draws <- suppressWarnings(as.integer(Sys.getenv("N_DRAWS", unset = "10")))
if (is.na(n_draws) || n_draws < 1L) n_draws <- 10L

available_cores <- parallel::detectCores(logical = TRUE)
if (is.na(available_cores) || available_cores < 1L) available_cores <- 1L

workers <- suppressWarnings(as.integer(Sys.getenv("N_WORKERS", unset = "8")))
if (is.na(workers) || workers < 1L) workers <- 1L
workers <- min(workers, available_cores)

# Input and output paths --------------------------------------------------------
# GitHub/repository recommendation:
#   data/raw/      input files, not committed if confidential or very large
#   outputs/tables generated summary tables
#
# The function below uses the first existing candidate path. It lets the script
# run on your Windows machine, on another computer using project-relative paths,
# or inside a temporary analysis environment such as /mnt/data.

path_first_existing <- function(..., label = "file") {
  candidates <- unlist(list(...), use.names = FALSE)
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]

  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) {
    return(existing[[1]])
  }

  warning(
    "No existing path found for ", label, ". Using first candidate:\n  ",
    candidates[[1]],
    "\nThe script will stop later if this required file is still missing.",
    call. = FALSE
  )

  candidates[[1]]
}

data_raw_dir <- file.path("data", "raw")
output_dir <- file.path("outputs")
tables_output_dir <- file.path(output_dir, "tables")
dir.create(tables_output_dir, showWarnings = FALSE, recursive = TRUE)

efficacy_draw_file <- path_first_existing(
  Sys.getenv("EFFICACY_DRAWS_FILE", unset = ""),
  file.path(data_raw_dir, "itn_efficacy_posterior_draws_100samples.csv"),
  "D:/fichier csv/itn_efficacy_posterior_draws_100samples.csv",
  label = "pyrethroid-only ITN efficacy posterior draws"
)

pbo_file <- path_first_existing(
  Sys.getenv("PBO_EFFICACY_DRAWS_FILE", unset = ""),
  file.path(data_raw_dir, "itn_efficacy_pbo_posterior_draws_100samples.csv"),
  "D:/fichier csv/itn_efficacy_pbo_posterior_draws_100samples.csv",
  label = "PBO ITN efficacy posterior draws"
)

gamma_draw_file <- path_first_existing(
  Sys.getenv("GAMMA_DRAWS_FILE", unset = ""),
  file.path(data_raw_dir, "gamma_p_first100_samples.csv"),
  "D:/fichier csv/gamma_p_first100_samples.csv",
  label = "gamma posterior draws"
)

excel_file <- path_first_existing(
  Sys.getenv("DRC_DATA_XLSX", unset = ""),
  file.path(data_raw_dir, "DRC_Data.xlsx"),
  "/mnt/data/DRC_Data.xlsx",
  "D:/Exel file/DRC_Data.xlsx",
  label = "DRC model input workbook"
)

dist_file <- path_first_existing(
  Sys.getenv("BEDNET_DISTRIBUTION_XLSX", unset = ""),
  file.path(data_raw_dir, "bednet_distribution_by_province_2008_2022.xlsx"),
  "D:/Exel file/bednet_distribution_by_province_2008_2022.xlsx",
  label = "bed-net distribution calendar"
)

usage_file <- path_first_existing(
  Sys.getenv("ITN_USAGE_DRAWS_FILE", unset = ""),
  file.path(data_raw_dir, "Fully_Corrected_ITN_Usage_Coverage_Data.csv"),
  "D:/fichier csv/Fully_Corrected_ITN_Usage_Coverage_Data.csv",
  label = "ITN usage posterior draws"
)

pbo_mix_xlsx <- path_first_existing(
  Sys.getenv("PBO_MIX_XLSX", unset = ""),
  file.path(data_raw_dir, "type_itn_pbo.xlsx"),
  "/mnt/data/type_itn_pbo.xlsx",
  "D:/Exel file/type_itn_pbo.xlsx",
  label = "PBO/pyrethroid ITN mix workbook"
)

who_burden_csv <- path_first_existing(
  Sys.getenv("WHO_BURDEN_CSV", unset = ""),
  file.path(data_raw_dir, "DRC_WHO_cases_deaths_2009_2024.csv"),
  "/mnt/data/DRC_WHO_cases_deaths_2009_2024.csv",
  "D:/fichier csv/DRC_WHO_cases_deaths_2009_2024.csv",
  label = "WHO burden CSV"
)

who_burden_xlsx <- path_first_existing(
  Sys.getenv("WHO_BURDEN_XLSX", unset = ""),
  file.path(data_raw_dir, "who_new_data_2024.xlsx"),
  "/mnt/data/who_new_data_2024.xlsx",
  "D:/Exel file/who_new_data_2024.xlsx",
  label = "WHO burden workbook"
)

# Manual fallback table. Fill this only when no WHO CSV/XLSX is available.
# Values must be absolute national counts, not millions.
manual_who_obs <- tibble::tribble(
  ~Year, ~WHO_cases, ~WHO_deaths,
  2009L, NA_real_, NA_real_,
  2010L, NA_real_, NA_real_,
  2011L, NA_real_, NA_real_,
  2012L, NA_real_, NA_real_,
  2013L, NA_real_, NA_real_,
  2014L, NA_real_, NA_real_,
  2015L, NA_real_, NA_real_,
  2016L, NA_real_, NA_real_,
  2017L, NA_real_, NA_real_,
  2018L, NA_real_, NA_real_,
  2019L, NA_real_, NA_real_,
  2020L, NA_real_, NA_real_,
  2021L, NA_real_, NA_real_,
  2022L, NA_real_, NA_real_,
  2023L, NA_real_, NA_real_,
  2024L, NA_real_, NA_real_
)


# =============================================================================
# 2. Parallel backend
# =============================================================================

try(future::plan(future::sequential), silent = TRUE)
gc()

future::plan(future::multisession, workers = workers)
options(
  future.rng.onMisuse = "ignore",
  future.globals.maxSize = 2 * 1024^3
)


# =============================================================================
# 3. Helper functions
# =============================================================================
norm_key <- function(x){
  x <- ifelse(is.na(x), "", as.character(x))
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- tolower(x)
  gsub("[^a-z0-9]", "", x)
}
intersect_all <- function(lst) Reduce(intersect, lst)
qL  <- function(x) quantile(x, 0.025, na.rm = TRUE, names = FALSE)
qU  <- function(x) quantile(x, 0.975, na.rm = TRUE, names = FALSE)
to01 <- function(x) ifelse(x > 1, x/100, x)
norm_year <- function(x) as.integer(x)

safe_parse_number <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  readr::parse_number(as.character(x), na = c("", "NA", "NaN", "null", "NULL"))
}

standardise_who_burden <- function(df, source_name = "") {
  if (is.null(df) || nrow(df) == 0) return(tibble::tibble())

  df <- janitor::clean_names(df)

  year_col <- intersect(
    c("year", "years", "annee", "annees", "calendar_year"),
    names(df)
  )[1]

  if (is.na(year_col)) return(tibble::tibble())

  # Prefer explicit WHO columns when present. Otherwise look for columns
  # containing "case" and "death".
  case_candidates <- intersect(
    c("who_cases", "cases", "malaria_cases", "estimated_cases",
      "estimated_malaria_cases", "total_cases", "cases_reported",
      "reported_cases", "cas", "cas_malaria"),
    names(df)
  )

  death_candidates <- intersect(
    c("who_deaths", "deaths", "malaria_deaths", "estimated_deaths",
      "estimated_malaria_deaths", "total_deaths", "deaths_reported",
      "reported_deaths", "deces", "deaths_malaria"),
    names(df)
  )

  if (length(case_candidates) == 0) {
    case_candidates <- names(df)[
      grepl("case|cas", names(df), ignore.case = TRUE) &
        !grepl("rate|ratio|cfr|per|percent|pct", names(df), ignore.case = TRUE)
    ]
  }

  if (length(death_candidates) == 0) {
    death_candidates <- names(df)[
      grepl("death|deces|mort", names(df), ignore.case = TRUE) &
        !grepl("rate|ratio|cfr|per|percent|pct", names(df), ignore.case = TRUE)
    ]
  }

  if (length(case_candidates) == 0 || length(death_candidates) == 0) {
    return(tibble::tibble())
  }

  out <- tibble::tibble(
    Year       = as.integer(safe_parse_number(df[[year_col]])),
    WHO_cases  = safe_parse_number(df[[case_candidates[1]]]),
    WHO_deaths = safe_parse_number(df[[death_candidates[1]]])
  ) %>%
    dplyr::filter(
      !is.na(Year),
      !is.na(WHO_cases),
      !is.na(WHO_deaths),
      WHO_cases > 0,
      WHO_deaths >= 0
    ) %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(
      WHO_cases  = dplyr::first(WHO_cases),
      WHO_deaths = dplyr::first(WHO_deaths),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Year)

  if (nrow(out) > 0) {
    message("WHO burden data read from: ", source_name)
  }

  out
}

standardise_who_multilevel_burden <- function(raw_df, source_name = "") {
  # This parser is for WHO workbooks with a two-row header, e.g.
  # DRC | Year | Population at risk | Cases: Lower/Point/Upper | Deaths: Lower/Point/Upper
  # It uses POINT estimates:
  #   Year        = column 2
  #   WHO_cases   = column 5
  #   WHO_deaths  = column 8
  if (is.null(raw_df) || nrow(raw_df) == 0 || ncol(raw_df) < 8) {
    return(tibble::tibble())
  }

  out <- tibble::tibble(
    Year       = as.integer(safe_parse_number(raw_df[[2]])),
    WHO_cases  = safe_parse_number(raw_df[[5]]),
    WHO_deaths = safe_parse_number(raw_df[[8]])
  ) %>%
    dplyr::filter(
      !is.na(Year),
      Year >= 1900,
      Year <= 2100,
      !is.na(WHO_cases),
      !is.na(WHO_deaths),
      WHO_cases > 0,
      WHO_deaths >= 0
    ) %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(
      WHO_cases  = dplyr::first(WHO_cases),
      WHO_deaths = dplyr::first(WHO_deaths),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Year)

  if (nrow(out) > 0) {
    message("WHO burden data read from WHO two-row header file: ", source_name)
  }

  out
}

read_who_burden <- function(
    who_csv,
    who_xlsx,
    drc_excel_file,
    manual_table,
    analysis_years
) {
  # 1) Prefer the uploaded/new WHO workbook.
  #    We first try the two-row header parser because the file
  #    who_new_data_2024.xlsx stores Cases and Deaths as Lower/Point/Upper.
  if (!is.null(who_xlsx) && file.exists(who_xlsx)) {
    xlsx_sheets <- readxl::excel_sheets(who_xlsx)

    for (sh in xlsx_sheets) {
      raw_no_header <- suppressMessages(readxl::read_excel(
        who_xlsx,
        sheet = sh,
        col_names = FALSE
      ))

      out <- standardise_who_multilevel_burden(
        raw_no_header,
        source_name = paste0("XLSX: ", who_xlsx, " | sheet: ", sh)
      )

      out <- out %>%
        dplyr::filter(Year %in% analysis_years) %>%
        dplyr::mutate(CFR_WHO = WHO_deaths / WHO_cases)

      if (nrow(out) > 0) return(out)
    }

    # If it was not the WHO two-row format, try ordinary named columns.
    for (sh in xlsx_sheets) {
      regular_sheet <- suppressMessages(readxl::read_excel(who_xlsx, sheet = sh))
      out <- standardise_who_burden(
        regular_sheet,
        source_name = paste0("XLSX: ", who_xlsx, " | sheet: ", sh)
      )

      out <- out %>%
        dplyr::filter(Year %in% analysis_years) %>%
        dplyr::mutate(CFR_WHO = WHO_deaths / WHO_cases)

      if (nrow(out) > 0) return(out)
    }
  }

  # 2) CSV fallback with ordinary named columns.
  if (!is.null(who_csv) && file.exists(who_csv)) {
    out <- standardise_who_burden(
      readr::read_csv(who_csv, show_col_types = FALSE),
      source_name = paste0("CSV: ", who_csv)
    ) %>%
      dplyr::filter(Year %in% analysis_years) %>%
      dplyr::mutate(CFR_WHO = WHO_deaths / WHO_cases)

    if (nrow(out) > 0) return(out)
  }

  # 3) Optional fallback: try to find a burden sheet inside DRC_Data.xlsx.
  if (!is.null(drc_excel_file) && file.exists(drc_excel_file)) {
    shs <- readxl::excel_sheets(drc_excel_file)
    shs_try <- shs[grepl("case|death|burden|who|malaria|data", shs, ignore.case = TRUE)]
    if (length(shs_try) == 0) shs_try <- shs

    for (sh in shs_try) {
      regular_sheet <- suppressMessages(readxl::read_excel(drc_excel_file, sheet = sh))
      out <- standardise_who_burden(
        regular_sheet,
        source_name = paste0("DRC_Data.xlsx: ", drc_excel_file, " | sheet: ", sh)
      ) %>%
        dplyr::filter(Year %in% analysis_years) %>%
        dplyr::mutate(CFR_WHO = WHO_deaths / WHO_cases)

      if (nrow(out) > 0) return(out)
    }
  }

  # 4) Manual fallback.
  manual_out <- standardise_who_burden(manual_table, source_name = "manual_who_obs") %>%
    dplyr::filter(Year %in% analysis_years) %>%
    dplyr::mutate(CFR_WHO = WHO_deaths / WHO_cases)

  if (nrow(manual_out) > 0) return(manual_out)

  stop(
    "WHO burden data were not found. Please use the file ",
    "'D:/Exel file/who_new_data_2024.xlsx' or create a CSV/XLSX file with ",
    "columns Year, WHO_cases, WHO_deaths. The with-resistance scenario cannot ",
    "be calibrated without WHO national burden values."
  )
}

extend_to_2024 <- function(df, year_col, cols_to_repeat){
  years_av <- sort(unique(df[[year_col]]))
  if (any(years_av == 2021) && !all(c(2022, 2023, 2024) %in% years_av)) {
    row2021 <- df %>% dplyr::filter(.data[[year_col]] == 2021) %>% dplyr::slice(1)
    ext <- dplyr::bind_rows(
      dplyr::mutate(row2021, !!year_col := 2022),
      dplyr::mutate(row2021, !!year_col := 2023),
      dplyr::mutate(row2021, !!year_col := 2024)
    )
    for (cc in cols_to_repeat) {
      if (cc %in% names(ext)) ext[[cc]] <- row2021[[cc]]
    }
    df <- dplyr::bind_rows(df, ext)
  }
  df
}

has_any_na <- function(x) {
  any(is.na(unlist(x)), na.rm = TRUE)
}

safe_run_sim <- function(timesteps, parameters, label = "") {
  if (has_any_na(parameters)) {
    message("Skipping simulation ", label, " because parameters contain NA.")
    return(NULL)
  }
  out <- try(run_simulation(timesteps = timesteps, parameters = parameters),
             silent = TRUE)
  if (inherits(out, "try-error")) {
    message("run_simulation failed for ", label, ": ", as.character(out))
    return(NULL)
  }
  out
}

empty_region_df <- tibble::tibble(
  Year                    = integer(),
  Cases_NoNets            = double(),
  Cases_WithNets_Resist   = double(),
  Cases_WithNets_PyOnly   = double(),
  Cases_WithNets_NoResist = double(),
  ProvKey                 = character(),
  Draw                    = integer()
)

collapse_eff <- function(df){
  if (!"rn0" %in% names(df) && "rn" %in% names(df)) {
    df <- df %>% dplyr::rename(rn0 = rn)
  }
  ycol <- intersect(c("year","Year"), names(df))[1]
  dcol <- intersect(c("dn0","DN0"), names(df))[1]
  rcol <- intersect(c("rn0","RN0"), names(df))[1]
  if (is.na(ycol) || is.na(dcol) || is.na(rcol)) {
    stop("collapse_eff: missing year/dn0/rn0 columns.")
  }
  df %>%
    dplyr::transmute(Year = norm_year(.data[[ycol]]),
                     dn0  = .data[[dcol]],
                     rn0  = .data[[rcol]]) %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(
      dn0 = median(dn0, na.rm = TRUE),
      rn0 = median(rn0, na.rm = TRUE),
      .groups = "drop"
    )
}

read_pbo_mix <- function(path) {
  x <- readxl::read_excel(path)
  ycol <- intersect(c("Perc_PBO_net","PBO_weight","%PBO","PBO","pbo_share","pbo","%y","y"),
                    names(x))[1]
  if (is.na(ycol)) {
    stop("PBO % column not found in ", path,
         ". Columns: ", paste(names(x), collapse = ", "))
  }
  yval <- suppressWarnings(as.numeric(x[[ycol]]))
  yval <- to01(yval)
  ycol_year <- intersect(c("Year","year","annee"), names(x))[1]
  if (is.na(ycol_year)) stop("Year column not found in PBO mix file.")
  
  tibble::tibble(
    Year = norm_year(x[[ycol_year]]),
    y    = pmin(pmax(yval, 0), 1)
  ) %>%
    dplyr::mutate(y = dplyr::case_when(
      Year < 2019 ~ 0,   # <-- PBO starts in 2019, not 2018
      Year > 2024 ~ 0,
      TRUE ~ y
    ))
}

# ----------------------------
# LOAD DATA
# ----------------------------
if (!file.exists(excel_file))   stop("DRC_Data.xlsx not found: ", excel_file)
if (!file.exists(dist_file))    stop("Distribution calendar not found: ", dist_file)
if (!file.exists(pbo_mix_xlsx)) stop("type_itn_pbo.xlsx not found: ", pbo_mix_xlsx)

efficacy_draws <- readr::read_csv(efficacy_draw_file, show_col_types = FALSE) %>%
  janitor::clean_names()
pbo_draws      <- readr::read_csv(pbo_file,           show_col_types = FALSE) %>%
  janitor::clean_names()
gamma_draws    <- readr::read_csv(gamma_draw_file,    show_col_types = FALSE) %>%
  janitor::clean_names()

# Standardise draw column names across posterior files
if (!"draw_id" %in% names(efficacy_draws) && "draw" %in% names(efficacy_draws)) {
  efficacy_draws <- efficacy_draws %>% dplyr::rename(draw_id = draw)
}
if (!"draw_id" %in% names(pbo_draws) && "draw" %in% names(pbo_draws)) {
  pbo_draws <- pbo_draws %>% dplyr::rename(draw_id = draw)
}
if (!"draw_id" %in% names(gamma_draws) && "draw" %in% names(gamma_draws)) {
  gamma_draws <- gamma_draws %>% dplyr::rename(draw_id = draw)
}
if (!"draw_id" %in% names(efficacy_draws)) efficacy_draws$draw_id <- 1L
if (!"draw_id" %in% names(pbo_draws))      pbo_draws$draw_id      <- 1L
if (!"draw_id" %in% names(gamma_draws))    gamma_draws$draw_id    <- 1L

usage_draws <- readr::read_csv(usage_file, show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    province = as.character(province),
    year     = as.integer(year),
    provkey  = norm_key(province)
  )

cov_candidates <- c("p","usage","coverage","cov","itn_usage","itn_coverage","u")
cov_col <- intersect(cov_candidates, names(usage_draws))[1]
if (is.na(cov_col)) {
  stop("No coverage/usage column found in usage_draws. Available: ",
       paste(names(usage_draws), collapse = ", "),
       " | expected one of: ", paste(cov_candidates, collapse = ", "))
}
usage_draws <- usage_draws %>%
  dplyr::mutate(
    cov_val_raw = suppressWarnings(as.numeric(.data[[cov_col]])),
    cov_val = dplyr::case_when(
      is.na(cov_val_raw) ~ NA_real_,
      cov_val_raw > 1    ~ cov_val_raw / 100,
      TRUE               ~ cov_val_raw
    ),
    cov_val = pmin(pmax(cov_val, 0), 1)
  )
if (!"sample_id" %in% names(usage_draws)) usage_draws$sample_id <- 1L

# Excel sheets
pop_raw            <- readxl::read_excel(excel_file, sheet = "Population")
interventions_data <- readxl::read_excel(excel_file, sheet = "Interventions")
eir_data           <- readxl::read_excel(excel_file, sheet = "EIR")
dist_calendar_raw  <- readxl::read_excel(dist_file) %>% dplyr::rename(province = 1)

if ("urban_rural" %in% names(interventions_data)) {
  interventions_data <- interventions_data %>% dplyr::filter(urban_rural == "rural")
}
interventions_data <- interventions_data %>%
  dplyr::mutate(
    name_1  = as.character(name_1),
    provkey = norm_key(name_1),
    year    = as.integer(year)
  )

eir_data <- eir_data %>%
  dplyr::mutate(
    name_1  = as.character(name_1),
    provkey = norm_key(name_1)
  )

dist_calendar <- dist_calendar_raw %>%
  dplyr::mutate(
    province = as.character(province),
    provkey  = norm_key(province)
  )
names(dist_calendar) <- gsub("^x", "", names(dist_calendar))

# Names dictionary
pop_names <- pop_raw %>%
  dplyr::mutate(
    name_1  = as.character(name_1),
    provkey = norm_key(name_1)
  ) %>%
  dplyr::distinct(provkey, name_1) %>%
  dplyr::rename(region = name_1)

cal_names <- dist_calendar %>%
  dplyr::group_by(provkey) %>%
  dplyr::summarise(region = dplyr::first(province), .groups = "drop")

name_dict <- dplyr::bind_rows(pop_names, cal_names) %>%
  dplyr::group_by(provkey) %>%
  dplyr::summarise(region = dplyr::first(na.omit(region)), .groups = "drop")

# Population (analysis years)
pop_df <- {
  p0 <- if ("urban_rural" %in% names(pop_raw) && any(pop_raw$urban_rural == "rural", na.rm = TRUE)) {
    pop_raw %>% dplyr::filter(urban_rural == "rural")
  } else {
    pop_raw
  }
  p0 %>%
    dplyr::transmute(
      provkey = norm_key(as.character(name_1)),
      Year    = as.integer(year),
      Pop     = as.numeric(pop)
    ) %>%
    dplyr::filter(!is.na(Pop), Year %in% ANALYSIS_YEARS) %>%
    dplyr::group_by(provkey, Year) %>%
    dplyr::summarise(Pop = sum(Pop, na.rm = TRUE), .groups = "drop")
}

# --- helpers: mean populations (2009–2024) for per 10 000 ---

# Province-level mean pop (sur toutes les années d'analyse)
prov_pop_mean <- pop_df %>%
  dplyr::group_by(ProvKey = provkey) %>%
  dplyr::summarise(
    MeanPop_2009_2024        = mean(Pop, na.rm = TRUE),
    TotalPopYears_2009_2024  = sum(Pop,  na.rm = TRUE),
    .groups = "drop"
  )

# National mean pop (sur toutes les années d'analyse)
nat_pop_mean <- pop_df %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(PopYear = sum(Pop, na.rm = TRUE), .groups = "drop") %>%
  dplyr::summarise(
    MeanNatPop_2009_2024       = mean(PopYear, na.rm = TRUE),
    TotalNatPopYears_2009_2024 = sum(PopYear,  na.rm = TRUE),
    .groups = "drop"
  )

# Regions to simulate
regions_to_plot <- intersect_all(list(
  unique(dist_calendar$provkey),
  unique(interventions_data$provkey),
  unique(eir_data$provkey)
)) |> sort()

# Draw IDs
draw_ids <- unique(efficacy_draws$draw_id)
if (length(draw_ids) == 0 && "draw" %in% names(efficacy_draws)) {
  draw_ids <- unique(efficacy_draws$draw)
}
if (length(draw_ids) == 0) draw_ids <- 1L
draw_ids <- draw_ids[seq_len(min(n_draws, length(draw_ids)))]

sample_ids <- unique(usage_draws$sample_id)
if (length(sample_ids) == 0) sample_ids <- seq_along(draw_ids)

# PBO mix
pbo_mix <- read_pbo_mix(pbo_mix_xlsx)

# ----------------------------
# PER-DRAW WORKER
# ----------------------------
safe_run_one_draw <- function(i) {
  draw <- draw_ids[i]
  message("Running draw ", draw, " (index ", i, ")")
  
  # --- Filter data for this draw ---
  eff_py   <- efficacy_draws %>% dplyr::filter(is.na(draw_id) | draw_id == draw)
  eff_pbo  <- pbo_draws      %>% dplyr::filter(is.na(draw_id) | draw_id == draw)
  gam_df0  <- gamma_draws    %>% dplyr::filter(is.na(draw_id) | draw_id == draw)
  
  sample_id <- sample_ids[(i - 1) %% length(sample_ids) + 1]
  usage_sample <- usage_draws %>%
    dplyr::filter(sample_id == !!sample_id)
  if (nrow(usage_sample) == 0) usage_sample <- usage_draws
  
  # Collapse efficacy + extend
  eff_py_med  <- collapse_eff(eff_py)
  eff_pbo_med <- collapse_eff(eff_pbo)
  eff_py_med  <- extend_to_2024(eff_py_med,  "Year", c("dn0","rn0"))
  eff_pbo_med <- extend_to_2024(eff_pbo_med, "Year", c("dn0","rn0"))
  
  gam_df0 <- extend_to_2024(gam_df0, "year", c("gamma_p","gamman"))
  if (!"gamma_p" %in% names(gam_df0) && "gamman" %in% names(gam_df0)) {
    gam_df0 <- gam_df0 %>% dplyr::rename(gamma_p = gamman)
  }
  
  # Per-year blended efficacy (py + PBO) 2009–2024
  eff_years <- tibble::tibble(Year = ANALYSIS_YEARS) %>%
    dplyr::left_join(
      eff_py_med %>% dplyr::rename(dn0_py = dn0, rn0_py = rn0),
      by = "Year"
    ) %>%
    dplyr::left_join(
      eff_pbo_med %>% dplyr::rename(dn0_pbo = dn0, rn0_pbo = rn0),
      by = "Year"
    ) %>%
    dplyr::left_join(pbo_mix, by = "Year") %>%
    dplyr::mutate(
      y = dplyr::coalesce(y, 0),
      dn0_blend = dplyr::if_else(
        Year >= 2019 & Year <= 2024,        # PBO is active during 2019–2024
        y * dn0_pbo + (1 - y) * dn0_py,
        dn0_py
      ),
      rn0_blend = dplyr::if_else(
        Year >= 2019 & Year <= 2024,
        y * rn0_pbo + (1 - y) * rn0_py,
        rn0_py
      )
    )
  
  # --- Per-province loop ---
  reg_out <- lapply(regions_to_plot, function(pk) {
    # EIR
    reg_eir <- eir_data %>% dplyr::filter(provkey == pk)
    EIR0 <- NA_real_
    for (cand in c("eir_drc","EIR_drc","eir","EIR")) {
      if (cand %in% names(reg_eir)) {
        v <- reg_eir[[cand]]
        v <- v[!is.na(v)]
        if (length(v) > 0) {
          EIR0 <- mean(v)
          break
        }
      }
    }
    if (is.na(EIR0) || !is.finite(EIR0) || EIR0 <= 0) {
      message("Skipping prov ", pk, " (draw ", draw, "): invalid EIR0 = ", EIR0)
      return(empty_region_df)
    }
    
    # Interventions & calendar
    reg_int <- interventions_data %>% dplyr::filter(provkey == pk)
    if (nrow(reg_int) == 0) {
      message("No interventions for prov ", pk, ", draw ", draw)
      return(empty_region_df)
    }
    
    cal_row <- dist_calendar %>% dplyr::filter(provkey == pk)
    if (nrow(cal_row) == 0) {
      message("No calendar row for prov ", pk, ", draw ", draw)
      return(empty_region_df)
    }
    
    flag_cols <- setdiff(names(cal_row), c("province","provkey"))
    flags_num <- suppressWarnings(sapply(cal_row[1, flag_cols], as.numeric))
    campaign_years <- sort(as.integer(names(flags_num)[which(flags_num == 1)]))
    campaign_years <- campaign_years[!is.na(campaign_years)]
    if (length(campaign_years) == 0) {
      message("No campaign years for prov ", pk, ", draw ", draw)
      return(empty_region_df)
    }
    
    # Gamma by year
    gam_small <- gam_df0 %>%
      dplyr::rename(Year = year) %>%
      dplyr::group_by(Year) %>%
      dplyr::summarise(
        gamma_p = median(gamma_p, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Usage by year (median)
    region_usage <- usage_sample %>%
      dplyr::filter(provkey == pk, year %in% campaign_years) %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(coverages = median(cov_val, na.rm = TRUE), .groups = "drop")
    
    if (nrow(region_usage) == 0) {
      region_usage <- tibble::tibble(
        year      = campaign_years,
        coverages = rep(0.6, length(campaign_years))
      )
    }
    
    # Fill missing campaign years
    missing_years <- setdiff(campaign_years, region_usage$year)
    if (length(missing_years) > 0) {
      last_cov <- tail(region_usage$coverages[!is.na(region_usage$coverages)], 1)
      if (length(last_cov) == 0) last_cov <- 0.6
      region_usage <- dplyr::bind_rows(
        region_usage,
        tibble::tibble(
          year      = missing_years,
          coverages = rep(last_cov, length(missing_years))
        )
      )
    }
    
    # Replace remaining NA in coverage
    if (any(is.na(region_usage$coverages))) {
      cov_default <- ifelse(
        all(is.na(region_usage$coverages)),
        0.6,
        stats::median(region_usage$coverages, na.rm = TRUE)
      )
      region_usage$coverages[is.na(region_usage$coverages)] <- cov_default
    }
    
    # Efficacy for campaign years
    # dn0_blend/rn0_blend = observed pyrethroid + PBO mix.
    # dn0_py/rn0_py       = pyrethroid-only counterfactual.
    eff_small <- eff_years %>%
      dplyr::select(Year, dn0_py, rn0_py, dn0_blend, rn0_blend) %>%
      dplyr::group_by(Year) %>%
      dplyr::summarise(
        dn0_py    = median(dn0_py,    na.rm = TRUE),
        rn0_py    = median(rn0_py,    na.rm = TRUE),
        dn0_blend = median(dn0_blend, na.rm = TRUE),
        rn0_blend = median(rn0_blend, na.rm = TRUE),
        .groups   = "drop"
      )
    
    usage_small <- region_usage %>%
      dplyr::rename(Year = year, usage = coverages)
    
    # Valid years where all data present
    valid_years <- intersect(intersect(campaign_years, eff_small$Year), usage_small$Year)
    valid_years <- intersect(valid_years, gam_small$Year)
    
    if (length(valid_years) == 0) {
      message("No valid years for prov ", pk, ", draw ", draw)
      return(empty_region_df)
    }
    
    region_params <- tibble::tibble(Year = sort(valid_years)) %>%
      dplyr::left_join(eff_small,   by = "Year") %>%
      dplyr::left_join(usage_small, by = "Year") %>%
      dplyr::left_join(gam_small,   by = "Year") %>%
      dplyr::filter(
        !is.na(dn0_py),
        !is.na(rn0_py),
        !is.na(dn0_blend),
        !is.na(rn0_blend),
        !is.na(usage),
        !is.na(gamma_p)
      )
    
    if (nrow(region_params) == 0) {
      message("Region params empty after NA filter for prov ", pk, ", draw ", draw)
      return(empty_region_df)
    }
    
    # Timesteps
    start_year          <- min(region_params$Year)
    bednet_timesteps    <- (region_params$Year - start_year) * days_per_year
    region_usage_vec    <- region_params$usage

    # Observed pyrethroid + PBO mix
    region_dn0_vec      <- region_params$dn0_blend
    region_rn0_vec      <- region_params$rn0_blend

    # Pyrethroid-only counterfactual
    region_dn0_py_vec   <- region_params$dn0_py
    region_rn0_py_vec   <- region_params$rn0_py

    region_gamma_vec    <- region_params$gamma_p
    
    # Base parameters (no nets)
    base <- get_parameters(list(human_population = human_population))
    base <- set_equilibrium(parameters = base, init_EIR = EIR0)
    if (has_any_na(base)) {
      message("Base parameters contain NA for prov ", pk, ", draw ", draw,
              ". Skipping region.")
      return(empty_region_df)
    }
    
    # Scenario 1: NO NETS
    out_no_nets <- safe_run_sim(
      timesteps  = sim_length,
      parameters = base,
      label      = paste0("no_nets prov=", pk, " draw=", draw)
    )
    if (is.null(out_no_nets)) return(empty_region_df)
    
    # Scenario 2: NETS WITH RESISTANCE (py + PBO blend)
    params_resist <- set_bednets(
      parameters = base,
      timesteps  = bednet_timesteps,
      coverages  = matrix(region_usage_vec, ncol = 1),
      retention  = 19.8 / 12 * days_per_year,
      dn0        = matrix(region_dn0_vec,   ncol = 1),
      rn         = matrix(region_rn0_vec,   ncol = 1),
      rnm        = matrix(rep(0.24, length(bednet_timesteps)), ncol = 1),
      gamman     = matrix(region_gamma_vec * days_per_year,    ncol = 1)
    )
    out_resist <- safe_run_sim(
      timesteps  = sim_length,
      parameters = params_resist,
      label      = paste0("resist prov=", pk, " draw=", draw)
    )
    if (is.null(out_resist)) return(empty_region_df)

    # Scenario 2B: NETS WITH RESISTANCE, PYRETHROID-ONLY COUNTERFACTUAL
    # This keeps the same ITN coverage and resistance decay but removes PBO.
    # It is required to estimate the additional impact of switching from
    # pyrethroid-only nets to the observed pyrethroid-PBO mix.
    params_pyonly <- set_bednets(
      parameters = base,
      timesteps  = bednet_timesteps,
      coverages  = matrix(region_usage_vec, ncol = 1),
      retention  = 19.8 / 12 * days_per_year,
      dn0        = matrix(region_dn0_py_vec, ncol = 1),
      rn         = matrix(region_rn0_py_vec, ncol = 1),
      rnm        = matrix(rep(0.24, length(bednet_timesteps)), ncol = 1),
      gamman     = matrix(region_gamma_vec * days_per_year, ncol = 1)
    )

    out_pyonly <- safe_run_sim(
      timesteps  = sim_length,
      parameters = params_pyonly,
      label      = paste0("pyonly prov=", pk, " draw=", draw)
    )
    if (is.null(out_pyonly)) return(empty_region_df)
    
    # Scenario 3: NETS WITHOUT RESISTANCE (fixed efficacy)
    fixed_dn0 <- 0.533
    fixed_rn  <- 0.56
    fixed_gam <- 2.64 * days_per_year
    
    params_noresist <- set_bednets(
      parameters = base,
      timesteps  = bednet_timesteps,
      coverages  = matrix(region_usage_vec, ncol = 1),
      retention  = 19.8 / 12 * days_per_year,
      dn0        = matrix(rep(fixed_dn0, length(bednet_timesteps)), ncol = 1),
      rn         = matrix(rep(fixed_rn,  length(bednet_timesteps)), ncol = 1),
      rnm        = matrix(rep(0.24,     length(bednet_timesteps)), ncol = 1),
      gamman     = matrix(rep(fixed_gam,length(bednet_timesteps)), ncol = 1)
    )
    out_noresist <- safe_run_sim(
      timesteps  = sim_length,
      parameters = params_noresist,
      label      = paste0("no_resist prov=", pk, " draw=", draw)
    )
    if (is.null(out_noresist)) return(empty_region_df)
    
    # Daily → Annual
    idx_to_year <- function(i, start_year) start_year + floor((i - 1L) / days_per_year)
    
    yrs_no     <- idx_to_year(seq_len(nrow(out_no_nets)),  start_year)
    yrs_res    <- idx_to_year(seq_len(nrow(out_resist)),   start_year)
    yrs_pyonly <- idx_to_year(seq_len(nrow(out_pyonly)),   start_year)
    yrs_nres   <- idx_to_year(seq_len(nrow(out_noresist)), start_year)
    
    df_no <- tibble::tibble(
      Year  = yrs_no,
      Cases = out_no_nets$n_detect_730_3650
    ) %>%
      dplyr::filter(Year %in% ANALYSIS_YEARS) %>%
      dplyr::group_by(Year) %>%
      dplyr::summarise(
        Cases_NoNets = sum(Cases, na.rm = TRUE),
        .groups = "drop"
      )
    
    df_res <- tibble::tibble(
      Year  = yrs_res,
      Cases = out_resist$n_detect_730_3650
    ) %>%
      dplyr::filter(Year %in% ANALYSIS_YEARS) %>%
      dplyr::group_by(Year) %>%
      dplyr::summarise(
        Cases_WithNets_Resist = sum(Cases, na.rm = TRUE),
        .groups = "drop"
      )
    
    df_pyonly <- tibble::tibble(
      Year  = yrs_pyonly,
      Cases = out_pyonly$n_detect_730_3650
    ) %>%
      dplyr::filter(Year %in% ANALYSIS_YEARS) %>%
      dplyr::group_by(Year) %>%
      dplyr::summarise(
        Cases_WithNets_PyOnly = sum(Cases, na.rm = TRUE),
        .groups = "drop"
      )

    df_nres <- tibble::tibble(
      Year  = yrs_nres,
      Cases = out_noresist$n_detect_730_3650
    ) %>%
      dplyr::filter(Year %in% ANALYSIS_YEARS) %>%
      dplyr::group_by(Year) %>%
      dplyr::summarise(
        Cases_WithNets_NoResist = sum(Cases, na.rm = TRUE),
        .groups = "drop"
      )
    
    df_yearly <- df_no %>%
      dplyr::full_join(df_res,    by = "Year") %>%
      dplyr::full_join(df_pyonly, by = "Year") %>%
      dplyr::full_join(df_nres,   by = "Year") %>%
      dplyr::mutate(
        ProvKey = pk,
        Draw    = as.integer(draw)
      )
    
    df_yearly
  }) # end lapply(provkey)
  
  dplyr::bind_rows(reg_out)
}

# ----------------------------
# RUN ALL DRAWS (PARALLEL + FALLBACK)
# ----------------------------
results_list <- try({
  furrr::future_map(
    seq_along(draw_ids),
    ~ safe_run_one_draw(.x),
    .options = furrr::furrr_options(seed = TRUE, scheduling = 1)
  )
}, silent = TRUE)

if (inherits(results_list, "try-error")) {
  message("Parallel run failed; retrying sequentially...")
  future::plan(future::sequential)
  results_list <- purrr::map(seq_along(draw_ids), ~ safe_run_one_draw(.x))
}

core_yearly <- dplyr::bind_rows(results_list)

# Attach region names
core_yearly <- core_yearly %>%
  dplyr::left_join(name_dict, by = c("ProvKey" = "provkey")) %>%
  dplyr::mutate(
    Region = ifelse(is.na(region), ProvKey, region)
  ) %>%
  dplyr::select(
    Region, ProvKey, Draw, Year,
    Cases_NoNets,
    Cases_WithNets_Resist,      # observed pyrethroid + PBO mix
    Cases_WithNets_PyOnly,      # pyrethroid-only counterfactual
    Cases_WithNets_NoResist
  )

# ----------------------------
# POPULATION-WEIGHTED CASES + WHO CALIBRATION
# ----------------------------
# Important change:
#   1) ITN and PBO plots use the unscaled population-weighted model outputs.
#   2) WHO case/death calibration is used only for the resistance-attributable
#      outputs:
#        - Case increase due to resistance
#        - Deaths caused by resistance
#
# This keeps the public-health burden caused by resistance anchored to the
# latest WHO national burden, without forcing every ITN/PBO counterfactual plot
# to be WHO-scaled.

# Read WHO national burden data
who_obs <- read_who_burden(
  who_csv        = who_burden_csv,
  who_xlsx       = who_burden_xlsx,
  drc_excel_file = excel_file,
  manual_table   = manual_who_obs,
  analysis_years = ANALYSIS_YEARS
)

latest_who_year <- max(who_obs$Year, na.rm = TRUE)
fit_years <- ANALYSIS_YEARS[ANALYSIS_YEARS <= latest_who_year]
missing_fit_years <- setdiff(fit_years, who_obs$Year)

if (length(missing_fit_years) > 0) {
  stop(
    "WHO burden data are missing for these fitted years: ",
    paste(missing_fit_years, collapse = ", "),
    ". Add them to the WHO burden CSV/XLSX or manual_who_obs."
  )
}

if (any(ANALYSIS_YEARS > latest_who_year)) {
  message(
    "Years after ", latest_who_year,
    " are treated as projections. The last WHO scale and CFR are carried forward."
  )
}

# 1) Raw population-weighted model outputs.
#    These are used for ITN and PBO plots.
core_yearly_pw_raw <- core_yearly %>%
  dplyr::inner_join(pop_df, by = c("ProvKey" = "provkey", "Year" = "Year")) %>%
  dplyr::mutate(
    # Convert simulated cases per 10,000 people to absolute province-level cases
    Cases_NoNets_s_raw            = pmax((Cases_NoNets            / human_population) * Pop, 0),
    Cases_WithNets_Resist_s_raw   = pmax((Cases_WithNets_Resist   / human_population) * Pop, 0),
    Cases_WithNets_PyOnly_s_raw   = pmax((Cases_WithNets_PyOnly   / human_population) * Pop, 0),
    Cases_WithNets_NoResist_s_raw = pmax((Cases_WithNets_NoResist / human_population) * Pop, 0)
  )

# 2) National raw status-quo burden by draw and year.
#    The status-quo is observed ITN use with resistance.
model_nat_statusquo_raw <- core_yearly_pw_raw %>%
  dplyr::group_by(Year, Draw) %>%
  dplyr::summarise(
    Model_Cases_WithNets_Resist_raw =
      sum(Cases_WithNets_Resist_s_raw, na.rm = TRUE),
    .groups = "drop"
  )

# 3) Draw-specific annual WHO calibration factors.
#    These are used only for resistance-attributable calculations.
calibration_factors_fit <- model_nat_statusquo_raw %>%
  dplyr::filter(Year <= latest_who_year) %>%
  dplyr::left_join(who_obs, by = "Year") %>%
  dplyr::mutate(
    case_scale = WHO_cases / pmax(Model_Cases_WithNets_Resist_raw, eps),
    calibration_type = "matched_to_WHO_for_resistance_outputs_only"
  )

# 4) For projection years, carry forward last available WHO scale and CFR.
projection_factors <- tibble::tibble()

if (any(model_nat_statusquo_raw$Year > latest_who_year)) {
  last_factor_by_draw <- calibration_factors_fit %>%
    dplyr::filter(Year == latest_who_year) %>%
    dplyr::select(Draw, case_scale, CFR_WHO)

  projection_factors <- model_nat_statusquo_raw %>%
    dplyr::filter(Year > latest_who_year) %>%
    dplyr::left_join(last_factor_by_draw, by = "Draw") %>%
    dplyr::mutate(
      WHO_cases  = NA_real_,
      WHO_deaths = NA_real_,
      calibration_type = paste0("projection_using_", latest_who_year, "_scale")
    )
}

calibration_factors <- calibration_factors_fit %>%
  dplyr::bind_rows(projection_factors) %>%
  dplyr::select(
    Year, Draw,
    Model_Cases_WithNets_Resist_raw,
    WHO_cases, WHO_deaths, CFR_WHO,
    case_scale, calibration_type
  )

# 5) Keep two versions of burden outputs:
#    *_s   = unscaled population-weighted model outputs, used for ITN/PBO plots.
#    *_who = WHO-calibrated outputs, used only for resistance-attributable plots.
core_yearly_pw <- core_yearly_pw_raw %>%
  dplyr::left_join(
    calibration_factors %>%
      dplyr::select(Year, Draw, case_scale, CFR_WHO, calibration_type),
    by = c("Year", "Draw")
  ) %>%
  dplyr::mutate(
    # Unscaled model outputs for ITN/PBO comparisons
    Cases_NoNets_s            = pmax(Cases_NoNets_s_raw, 0),
    Cases_WithNets_Resist_s   = pmax(Cases_WithNets_Resist_s_raw, 0),
    Cases_WithNets_PyOnly_s   = pmax(Cases_WithNets_PyOnly_s_raw, 0),
    Cases_WithNets_NoResist_s = pmax(Cases_WithNets_NoResist_s_raw, 0),

    # Model-implied deaths without WHO case matching.
    # The annual CFR is used only as a conversion rate; the case totals remain raw.
    Deaths_NoNets_s            = pmax(Cases_NoNets_s            * CFR_WHO, 0),
    Deaths_WithNets_Resist_s   = pmax(Cases_WithNets_Resist_s   * CFR_WHO, 0),
    Deaths_WithNets_PyOnly_s   = pmax(Cases_WithNets_PyOnly_s   * CFR_WHO, 0),
    Deaths_WithNets_NoResist_s = pmax(Cases_WithNets_NoResist_s * CFR_WHO, 0),

    # WHO-calibrated outputs for resistance-attributable burden only
    Cases_NoNets_who            = pmax(Cases_NoNets_s_raw            * case_scale, 0),
    Cases_WithNets_Resist_who   = pmax(Cases_WithNets_Resist_s_raw   * case_scale, 0),
    Cases_WithNets_PyOnly_who   = pmax(Cases_WithNets_PyOnly_s_raw   * case_scale, 0),
    Cases_WithNets_NoResist_who = pmax(Cases_WithNets_NoResist_s_raw * case_scale, 0),

    Deaths_NoNets_who            = pmax(Cases_NoNets_who            * CFR_WHO, 0),
    Deaths_WithNets_Resist_who   = pmax(Cases_WithNets_Resist_who   * CFR_WHO, 0),
    Deaths_WithNets_PyOnly_who   = pmax(Cases_WithNets_PyOnly_who   * CFR_WHO, 0),
    Deaths_WithNets_NoResist_who = pmax(Cases_WithNets_NoResist_who * CFR_WHO, 0)
  )

# 6) Calibration checks: for fitted years, the WHO-calibrated status-quo
#    resistance scenario should match national WHO cases/deaths.
calibration_check <- core_yearly_pw %>%
  dplyr::group_by(Year, Draw) %>%
  dplyr::summarise(
    Model_cases_with_resistance_who =
      sum(Cases_WithNets_Resist_who, na.rm = TRUE),
    Model_deaths_with_resistance_who =
      sum(Deaths_WithNets_Resist_who, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(who_obs, by = "Year") %>%
  dplyr::mutate(
    cases_ratio_to_WHO = Model_cases_with_resistance_who / WHO_cases,
    deaths_ratio_to_WHO = Model_deaths_with_resistance_who / WHO_deaths
  )

calibration_check_summary <- calibration_check %>%
  dplyr::filter(!is.na(WHO_cases), !is.na(WHO_deaths)) %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(
    median_cases_ratio_to_WHO = median(cases_ratio_to_WHO, na.rm = TRUE),
    median_deaths_ratio_to_WHO = median(deaths_ratio_to_WHO, na.rm = TRUE),
    lower_cases_ratio_to_WHO = qL(cases_ratio_to_WHO),
    upper_cases_ratio_to_WHO = qU(cases_ratio_to_WHO),
    lower_deaths_ratio_to_WHO = qL(deaths_ratio_to_WHO),
    upper_deaths_ratio_to_WHO = qU(deaths_ratio_to_WHO),
    .groups = "drop"
  )

message("Calibration check: WHO-calibrated resistance outputs should have median case/death ratios close to 1 for fitted WHO years.")
print(calibration_check_summary)

# ----------------------------
# PROVINCE TOTALS (2009–2024, POP-WEIGHTED)
# ----------------------------
prov_totals_pw <- core_yearly_pw %>%
  dplyr::group_by(Region, ProvKey, Draw) %>%
  dplyr::summarise(
    # Raw cases for ITN/PBO outputs
    Sum_NoNets_s            = sum(Cases_NoNets_s,            na.rm = TRUE),
    Sum_WithNets_Resist_s   = sum(Cases_WithNets_Resist_s,   na.rm = TRUE),
    Sum_WithNets_PyOnly_s   = sum(Cases_WithNets_PyOnly_s,   na.rm = TRUE),
    Sum_WithNets_NoResist_s = sum(Cases_WithNets_NoResist_s, na.rm = TRUE),

    # Raw model-implied deaths for ITN/PBO outputs
    Sum_Deaths_NoNets_s            = sum(Deaths_NoNets_s,            na.rm = TRUE),
    Sum_Deaths_WithNets_Resist_s   = sum(Deaths_WithNets_Resist_s,   na.rm = TRUE),
    Sum_Deaths_WithNets_PyOnly_s   = sum(Deaths_WithNets_PyOnly_s,   na.rm = TRUE),
    Sum_Deaths_WithNets_NoResist_s = sum(Deaths_WithNets_NoResist_s, na.rm = TRUE),

    # WHO-calibrated burden for resistance-attributable outputs only
    Sum_WithNets_Resist_who   = sum(Cases_WithNets_Resist_who,   na.rm = TRUE),
    Sum_WithNets_NoResist_who = sum(Cases_WithNets_NoResist_who, na.rm = TRUE),
    Sum_Deaths_WithNets_Resist_who   = sum(Deaths_WithNets_Resist_who,   na.rm = TRUE),
    Sum_Deaths_WithNets_NoResist_who = sum(Deaths_WithNets_NoResist_who, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  dplyr::left_join(prov_pop_mean, by = "ProvKey") %>%
  dplyr::mutate(
    # Lives saved due to ITNs: raw/unscaled model comparison
    LivesSaved_ITN_abs =
      pmax(Sum_Deaths_NoNets_s - Sum_Deaths_WithNets_Resist_s, 0),

    # Cases/deaths caused by resistance: WHO-calibrated comparison only
    CasesCaused_Resist_abs =
      pmax(Sum_WithNets_Resist_who - Sum_WithNets_NoResist_who, 0),
    DeathsCaused_Resist_abs =
      pmax(Sum_Deaths_WithNets_Resist_who - Sum_Deaths_WithNets_NoResist_who, 0),

    LivesSaved_ITN_per10k =
      LivesSaved_ITN_abs / pmax(MeanPop_2009_2024 / 10000, eps),
    DeathsCaused_Resist_per10k =
      DeathsCaused_Resist_abs / pmax(MeanPop_2009_2024 / 10000, eps)
  )

# ----------------------------
# NATIONAL YEARLY (POP-WEIGHTED)
# ----------------------------
national_yearly_pw <- core_yearly_pw %>%
  dplyr::group_by(Year, Draw) %>%
  dplyr::summarise(
    # Raw cases for ITN/PBO outputs
    Sum_NoNets_s            = sum(Cases_NoNets_s,            na.rm = TRUE),
    Sum_WithNets_Resist_s   = sum(Cases_WithNets_Resist_s,   na.rm = TRUE),
    Sum_WithNets_PyOnly_s   = sum(Cases_WithNets_PyOnly_s,   na.rm = TRUE),
    Sum_WithNets_NoResist_s = sum(Cases_WithNets_NoResist_s, na.rm = TRUE),

    # Raw model-implied deaths for ITN/PBO outputs
    Sum_Deaths_NoNets_s            = sum(Deaths_NoNets_s,            na.rm = TRUE),
    Sum_Deaths_WithNets_Resist_s   = sum(Deaths_WithNets_Resist_s,   na.rm = TRUE),
    Sum_Deaths_WithNets_PyOnly_s   = sum(Deaths_WithNets_PyOnly_s,   na.rm = TRUE),
    Sum_Deaths_WithNets_NoResist_s = sum(Deaths_WithNets_NoResist_s, na.rm = TRUE),

    # WHO-calibrated burden for resistance-attributable outputs only
    Sum_WithNets_Resist_who   = sum(Cases_WithNets_Resist_who,   na.rm = TRUE),
    Sum_WithNets_NoResist_who = sum(Cases_WithNets_NoResist_who, na.rm = TRUE),
    Sum_Deaths_WithNets_Resist_who   = sum(Deaths_WithNets_Resist_who,   na.rm = TRUE),
    Sum_Deaths_WithNets_NoResist_who = sum(Deaths_WithNets_NoResist_who, na.rm = TRUE),

    .groups = "drop"
  )

# ----------------------------
# NATIONAL TOTAL (2009–2024, POP-WEIGHTED)
# ----------------------------
national_total_pw <- national_yearly_pw %>%
  dplyr::group_by(Draw) %>%
  dplyr::summarise(
    # Raw cases for ITN/PBO outputs
    Nat_Sum_NoNets_s            = sum(Sum_NoNets_s,            na.rm = TRUE),
    Nat_Sum_WithNets_Resist_s   = sum(Sum_WithNets_Resist_s,   na.rm = TRUE),
    Nat_Sum_WithNets_PyOnly_s   = sum(Sum_WithNets_PyOnly_s,   na.rm = TRUE),
    Nat_Sum_WithNets_NoResist_s = sum(Sum_WithNets_NoResist_s, na.rm = TRUE),

    # Raw model-implied deaths for ITN/PBO outputs
    Nat_Sum_Deaths_NoNets_s            = sum(Sum_Deaths_NoNets_s,            na.rm = TRUE),
    Nat_Sum_Deaths_WithNets_Resist_s   = sum(Sum_Deaths_WithNets_Resist_s,   na.rm = TRUE),
    Nat_Sum_Deaths_WithNets_PyOnly_s   = sum(Sum_Deaths_WithNets_PyOnly_s,   na.rm = TRUE),
    Nat_Sum_Deaths_WithNets_NoResist_s = sum(Sum_Deaths_WithNets_NoResist_s, na.rm = TRUE),

    # WHO-calibrated burden for resistance-attributable outputs only
    Nat_Sum_WithNets_Resist_who   = sum(Sum_WithNets_Resist_who,   na.rm = TRUE),
    Nat_Sum_WithNets_NoResist_who = sum(Sum_WithNets_NoResist_who, na.rm = TRUE),
    Nat_Sum_Deaths_WithNets_Resist_who   = sum(Sum_Deaths_WithNets_Resist_who,   na.rm = TRUE),
    Nat_Sum_Deaths_WithNets_NoResist_who = sum(Sum_Deaths_WithNets_NoResist_who, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  tidyr::crossing(nat_pop_mean) %>%
  dplyr::mutate(
    # Lives saved due to ITNs: raw/unscaled model comparison
    Nat_LivesSaved_ITN_abs =
      pmax(Nat_Sum_Deaths_NoNets_s - Nat_Sum_Deaths_WithNets_Resist_s, 0),

    # Cases/deaths caused by resistance: WHO-calibrated comparison only
    Nat_CasesCaused_Resist_abs =
      pmax(Nat_Sum_WithNets_Resist_who - Nat_Sum_WithNets_NoResist_who, 0),
    Nat_DeathsCaused_Resist_abs =
      pmax(Nat_Sum_Deaths_WithNets_Resist_who - Nat_Sum_Deaths_WithNets_NoResist_who, 0),

    Nat_LivesSaved_ITN_per10k =
      Nat_LivesSaved_ITN_abs / pmax(MeanNatPop_2009_2024 / 10000, eps),
    Nat_DeathsCaused_Resist_per10k =
      Nat_DeathsCaused_Resist_abs / pmax(MeanNatPop_2009_2024 / 10000, eps)
  )


# ----------------------------
# PBO period for 'cases averted' & 'lives saved due to PBO'
# ----------------------------

pbo_years_active <- pbo_mix %>%
  dplyr::filter(!is.na(y), y > 0, Year %in% ANALYSIS_YEARS) %>%
  dplyr::pull(Year) %>%
  unique()

if (length(pbo_years_active) == 0) {
  warning("No year with PBO share > 0 was found in pbo_mix; PBO-specific calculations are skipped.")
} else {
  PBO_START_YEAR <- min(pbo_years_active)
  PBO_END_YEAR   <- min(2024L, max(pbo_years_active))
  
  message(
    "PBO period for 'cases averted' & 'lives saved due to PBO': ",
    PBO_START_YEAR, "–", PBO_END_YEAR
  )
  
  # Subset model outputs to the PBO period
  core_pbo_period <- core_yearly_pw %>%
    dplyr::filter(Year >= PBO_START_YEAR, Year <= PBO_END_YEAR)
  
  # Mean populations over the PBO period only
  prov_pop_pbo_mean <- pop_df %>%
    dplyr::filter(Year >= PBO_START_YEAR, Year <= PBO_END_YEAR) %>%
    dplyr::group_by(ProvKey = provkey) %>%
    dplyr::summarise(
      MeanPop_PBO       = mean(Pop, na.rm = TRUE),
      TotalPopYears_PBO = sum(Pop,  na.rm = TRUE),
      .groups = "drop"
    )
  
  nat_pop_pbo_mean <- pop_df %>%
    dplyr::filter(Year >= PBO_START_YEAR, Year <= PBO_END_YEAR) %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(PopYear = sum(Pop, na.rm = TRUE), .groups = "drop") %>%
    dplyr::summarise(
      MeanNatPop_PBO       = mean(PopYear, na.rm = TRUE),
      TotalNatPopYears_PBO = sum(PopYear,  na.rm = TRUE),
      .groups = "drop"
    )
  
  # Province-level totals over the PBO period
  prov_pbo_totals_pw <- core_pbo_period %>%
    dplyr::group_by(Region, ProvKey, Draw) %>%
    dplyr::summarise(
      # Cases summed over PBO years, absolute counts
      Cases_NoNets_PBOperiod          = sum(Cases_NoNets_s,          na.rm = TRUE),
      Cases_WithNets_Resist_PBOperiod = sum(Cases_WithNets_Resist_s, na.rm = TRUE),
      # Deaths summed over PBO years, absolute counts
      Deaths_NoNets_PBOperiod          = sum(Deaths_NoNets_s,          na.rm = TRUE),
      Deaths_WithNets_Resist_PBOperiod = sum(Deaths_WithNets_Resist_s, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # Cases averted from the start of PBO use, absolute counts
      Cases_Averted_PBOperiod =
        Cases_NoNets_PBOperiod - Cases_WithNets_Resist_PBOperiod,
      # Lives saved due to PBO over the PBO period, absolute counts
      LivesSaved_PBOperiod =
        Deaths_NoNets_PBOperiod - Deaths_WithNets_Resist_PBOperiod
    ) %>%
    # Values per 10,000 population over the PBO period
    dplyr::left_join(prov_pop_pbo_mean, by = "ProvKey") %>%
    dplyr::mutate(
      LivesSaved_PBOperiod_per10k =
        LivesSaved_PBOperiod / (MeanPop_PBO / 10000)
    )
  
  # National yearly sur la période PBO
  nat_pbo_yearly_pw <- national_yearly_pw %>%
    dplyr::filter(Year >= PBO_START_YEAR, Year <= PBO_END_YEAR) %>%
    dplyr::mutate(
      # Cases averted vs No nets sur la période PBO
      Cases_Averted_PBO = Sum_NoNets_s - Sum_WithNets_Resist_s,
      # Lives saved vs No nets sur la période PBO
      LivesSaved_PBO    = Sum_Deaths_NoNets_s - Sum_Deaths_WithNets_Resist_s
    )
  
  # National total PBO (somme sur PBO_START_YEAR–2024)
  nat_pbo_total_pw <- nat_pbo_yearly_pw %>%
    dplyr::group_by(Draw) %>%
    dplyr::summarise(
      Nat_Cases_Averted_PBO = sum(Cases_Averted_PBO, na.rm = TRUE),
      Nat_LivesSaved_PBO    = sum(LivesSaved_PBO,    na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Add national mean population over the PBO period
    tidyr::crossing(nat_pop_pbo_mean) %>%
    dplyr::mutate(
      Nat_LivesSaved_PBO_per10k =
        Nat_LivesSaved_PBO / (MeanNatPop_PBO / 10000)
    )
  
  message("
PBO-specific objects created:
- prov_pbo_totals_pw : province/draw, cases and deaths over [", PBO_START_YEAR, "–", PBO_END_YEAR, "]
  * Cases_Averted_PBOperiod
  * LivesSaved_PBOperiod (absolu) + LivesSaved_PBOperiod_per10k
- nat_pbo_yearly_pw  : national, par année sur cette période
  * Cases_Averted_PBO, LivesSaved_PBO
- nat_pbo_total_pw   : national, totals over the PBO period
  * Nat_Cases_Averted_PBO, Nat_LivesSaved_PBO, Nat_LivesSaved_PBO_per10k
")
}


# ============================================================
# FINAL PLOT SECTION
# ============================================================

# ============================================================
# PLOTS + SUMMARIES: PROVINCES, FOUR PROVINCES + DRC, NATIONAL
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(scales)
})

# ------------------------------------------------------------
# THEME + FORMATTING + COLOURS
# ------------------------------------------------------------

base_theme_box <- theme_bw(base_size = 20) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )

fmt_pct1 <- function(x) sprintf("%.1f", x)
fmt_int  <- function(x) format(round(x), big.mark = " ", scientific = FALSE)

# four_labels is assumed to exist already, but we re-declare to be safe
four_labels <- c("Haut-Katanga", "Kasaï-Central", "Kinshasa", "Sud-Kivu")

# Palette for provinces (including others as in your earlier code)
palette_table <- tibble::tribble(
  ~Province,        ~FillColor, ~PointColor, ~n_sites,
  "Equateur",       "#7CE3BB",  "#1B9E77",   3,
  "Haut-Katanga",   "#FFBCAB",  "#D95F02",   16,
  "Ituri",          "#C4C0FF",  "#7570B3",   3,
  "Kasaï-Central",  "#FFB4CE",  "#E7298A",   73,
  "Kasaï-Oriental", "#A5E575",  "#66A61E",   3,
  "Kinshasa",       "#FFDBAA",  "#E6AB02",   19,
  "Kongo-Central",  "#F5C17E",  "#A6761D",   6,
  "Kwilu",          "#BFBFBF",  "#666666",   4,
  "Maniema",        "#FFFFE0",  "#FFFFB3",   5,
  "Maï-Ndombe",     "#B0F5E9",  "#8DD3C7",   3,
  "Nord-Kivu",      "#E5E1FF",  "#BEBADA",   3,
  "Nord-Ubangi",    "#FFCECB",  "#FB8072",   9,
  "Sankuru",        "#BBE2FF",  "#80B1D3",   9,
  "Sud-Kivu",       "#FFE1CB",  "#FDB462",   18,
  "Tanganyika",     "#CFFB89",  "#B3DE69",   6,
  "Tshopo",         "#FFEBF5",  "#FCCDE5",   13
)

four_palette <- palette_table %>% filter(Province %in% four_labels)

col_four_regions_fill  <- four_palette$FillColor
names(col_four_regions_fill) <- four_palette$Province
# DRC (mean) in blue
col_four_regions_fill  <- c(col_four_regions_fill, "DRC (mean)" = "#0072B2")

col_four_regions_point <- four_palette$PointColor
names(col_four_regions_point) <- four_palette$Province
col_four_regions_point <- c(col_four_regions_point, "DRC (mean)" = "#0072B2")

# National overtime bars in blue
col_nat_overtime <- "#0072B2"

# small epsilon
eps <- 1e-9
eps_bar <- 1e-6

qL <- function(x) quantile(x, 0.025, na.rm = TRUE, names = FALSE)
qU <- function(x) quantile(x, 0.975, na.rm = TRUE, names = FALSE)

# Repair uncertainty intervals before plotting.
# This avoids negative uncertainty bars and prevents lower bounds from sitting
# exactly on zero when a quantity must be non-negative.
fix_ci <- function(data, median_col, lower_col, upper_col, floor = 0, cap = Inf) {
  med <- as.numeric(data[[median_col]])
  lo  <- as.numeric(data[[lower_col]])
  up  <- as.numeric(data[[upper_col]])

  med <- pmin(pmax(med, floor), cap)
  lo  <- pmin(pmax(lo,  floor), cap)
  up  <- pmin(pmax(up,  floor), cap)

  lo <- pmin(lo, med)
  up <- pmax(up, med)

  data[[median_col]] <- med
  data[[lower_col]]  <- lo
  data[[upper_col]]  <- up
  data
}

# Floors used only for display/summary repair. They are intentionally tiny
# so they remove impossible negative/zero intervals without changing the result.
pct_floor   <- 0.01
count_floor <- 1e-6
per10k_floor <- 1e-6

theme_x_oblique <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1)
)

# ============================================================
# 1) PROVINCE + NATIONAL SUMMARIES (ITNs / RESISTANCE)
# ============================================================

# prov_totals_pw is assumed to exist (from your core script)
# with columns:
#  Region, ProvKey, Draw,
#  Sum_NoNets_s, Sum_WithNets_Resist_s, Sum_WithNets_NoResist_s,
#  Sum_Deaths_NoNets_s, Sum_Deaths_WithNets_Resist_s, Sum_Deaths_WithNets_NoResist_s,
#  LivesSaved_ITN_abs, DeathsCaused_Resist_abs,
#  LivesSaved_ITN_per10k, DeathsCaused_Resist_per10k

prov_summary <- prov_totals_pw %>%
  group_by(Region, ProvKey) %>%
  summarise(
    # ITN effects use the raw/unscaled model outputs
    PercAverted_ITNs = 100 *
      pmin(
        pmax(Sum_NoNets_s - Sum_WithNets_Resist_s, 0) /
          pmax(Sum_NoNets_s, eps),
        1
      ),

    LivesSaved_ITNs_abs    = pmax(LivesSaved_ITN_abs, 0),
    LivesSaved_ITNs_per10k = pmax(LivesSaved_ITN_per10k, 0),

    # Resistance-attributable effects use WHO-calibrated outputs only
    PercIncrease_Resist = 100 *
      pmax(Sum_WithNets_Resist_who - Sum_WithNets_NoResist_who, 0) /
      pmax(Sum_WithNets_NoResist_who, eps),

    DeathsCaused_Resist_abs    = pmax(DeathsCaused_Resist_abs, 0),
    DeathsCaused_Resist_per10k = pmax(DeathsCaused_Resist_per10k, 0),

    .groups = "drop_last"
  ) %>%
  summarise(
    # % cases averted ITNs
    PercAverted_ITNs_median = median(PercAverted_ITNs, na.rm = TRUE),
    PercAverted_ITNs_lower  = qL(PercAverted_ITNs),
    PercAverted_ITNs_upper  = qU(PercAverted_ITNs),

    # % case increase due to resistance
    PercIncrease_Resist_median = median(PercIncrease_Resist, na.rm = TRUE),
    PercIncrease_Resist_lower  = qL(PercIncrease_Resist),
    PercIncrease_Resist_upper  = qU(PercIncrease_Resist),

    # Lives saved due to ITNs (absolute)
    LivesSaved_ITNs_median = median(LivesSaved_ITNs_abs, na.rm = TRUE),
    LivesSaved_ITNs_lower  = qL(LivesSaved_ITNs_abs),
    LivesSaved_ITNs_upper  = qU(LivesSaved_ITNs_abs),

    # Lives saved per 10 000
    LivesSaved_ITNs_per10k_median = median(LivesSaved_ITNs_per10k, na.rm = TRUE),
    LivesSaved_ITNs_per10k_lower  = qL(LivesSaved_ITNs_per10k),
    LivesSaved_ITNs_per10k_upper  = qU(LivesSaved_ITNs_per10k),

    # Deaths caused by resistance (absolute)
    DeathsCaused_Resist_median = median(DeathsCaused_Resist_abs, na.rm = TRUE),
    DeathsCaused_Resist_lower  = qL(DeathsCaused_Resist_abs),
    DeathsCaused_Resist_upper  = qU(DeathsCaused_Resist_abs),

    # Deaths per 10 000
    DeathsCaused_Resist_per10k_median =
      median(DeathsCaused_Resist_per10k, na.rm = TRUE),
    DeathsCaused_Resist_per10k_lower  = qL(DeathsCaused_Resist_per10k),
    DeathsCaused_Resist_per10k_upper  = qU(DeathsCaused_Resist_per10k),

    .groups = "drop"
  ) %>%
  fix_ci("PercAverted_ITNs_median", "PercAverted_ITNs_lower", "PercAverted_ITNs_upper",
         floor = pct_floor, cap = 100) %>%
  fix_ci("PercIncrease_Resist_median", "PercIncrease_Resist_lower", "PercIncrease_Resist_upper",
         floor = pct_floor) %>%
  fix_ci("LivesSaved_ITNs_median", "LivesSaved_ITNs_lower", "LivesSaved_ITNs_upper",
         floor = count_floor) %>%
  fix_ci("LivesSaved_ITNs_per10k_median", "LivesSaved_ITNs_per10k_lower", "LivesSaved_ITNs_per10k_upper",
         floor = per10k_floor) %>%
  fix_ci("DeathsCaused_Resist_median", "DeathsCaused_Resist_lower", "DeathsCaused_Resist_upper",
         floor = count_floor) %>%
  fix_ci("DeathsCaused_Resist_per10k_median", "DeathsCaused_Resist_per10k_lower", "DeathsCaused_Resist_per10k_upper",
         floor = per10k_floor)

# national_total_pw assumed to exist with Nat_* columns as before
nat_total_summary <- national_total_pw %>%
  mutate(
    # ITN effects use raw/unscaled model outputs
    Nat_PercAverted_ITNs = 100 *
      pmin(
        pmax(Nat_Sum_NoNets_s - Nat_Sum_WithNets_Resist_s, 0) /
          pmax(Nat_Sum_NoNets_s, eps),
        1
      ),

    # Resistance-attributable effects use WHO-calibrated outputs only
    Nat_PercIncrease_Resist = 100 *
      pmax(Nat_Sum_WithNets_Resist_who - Nat_Sum_WithNets_NoResist_who, 0) /
      pmax(Nat_Sum_WithNets_NoResist_who, eps)
  ) %>%
  summarise(
    # % cases averted ITNs
    PercAverted_ITNs_median = median(Nat_PercAverted_ITNs, na.rm = TRUE),
    PercAverted_ITNs_lower  = qL(Nat_PercAverted_ITNs),
    PercAverted_ITNs_upper  = qU(Nat_PercAverted_ITNs),

    # % case increase due to resistance
    PercIncrease_Resist_median = median(Nat_PercIncrease_Resist, na.rm = TRUE),
    PercIncrease_Resist_lower  = qL(Nat_PercIncrease_Resist),
    PercIncrease_Resist_upper  = qU(Nat_PercIncrease_Resist),

    # Lives saved ITNs (absolute national; raw/unscaled)
    Nat_LivesSaved_ITNs_median = median(Nat_LivesSaved_ITN_abs, na.rm = TRUE),
    Nat_LivesSaved_ITNs_lower  = qL(Nat_LivesSaved_ITN_abs),
    Nat_LivesSaved_ITNs_upper  = qU(Nat_LivesSaved_ITN_abs),

    # Lives saved ITNs per 10k (national; raw/unscaled)
    Nat_LivesSaved_ITNs_per10k_median =
      median(Nat_LivesSaved_ITN_per10k, na.rm = TRUE),
    Nat_LivesSaved_ITNs_per10k_lower  = qL(Nat_LivesSaved_ITN_per10k),
    Nat_LivesSaved_ITNs_per10k_upper  = qU(Nat_LivesSaved_ITN_per10k),

    # Deaths caused by resistance (absolute national; WHO-calibrated)
    Nat_DeathsCaused_Resist_median =
      median(Nat_DeathsCaused_Resist_abs, na.rm = TRUE),
    Nat_DeathsCaused_Resist_lower  = qL(Nat_DeathsCaused_Resist_abs),
    Nat_DeathsCaused_Resist_upper  = qU(Nat_DeathsCaused_Resist_abs),

    # Deaths per 10k (national; WHO-calibrated)
    Nat_DeathsCaused_Resist_per10k_median =
      median(Nat_DeathsCaused_Resist_per10k, na.rm = TRUE),
    Nat_DeathsCaused_Resist_per10k_lower  = qL(Nat_DeathsCaused_Resist_per10k),
    Nat_DeathsCaused_Resist_per10k_upper  = qU(Nat_DeathsCaused_Resist_per10k),

    .groups = "drop"
  ) %>%
  fix_ci("PercAverted_ITNs_median", "PercAverted_ITNs_lower", "PercAverted_ITNs_upper",
         floor = pct_floor, cap = 100) %>%
  fix_ci("PercIncrease_Resist_median", "PercIncrease_Resist_lower", "PercIncrease_Resist_upper",
         floor = pct_floor) %>%
  fix_ci("Nat_LivesSaved_ITNs_median", "Nat_LivesSaved_ITNs_lower", "Nat_LivesSaved_ITNs_upper",
         floor = count_floor) %>%
  fix_ci("Nat_LivesSaved_ITNs_per10k_median", "Nat_LivesSaved_ITNs_per10k_lower", "Nat_LivesSaved_ITNs_per10k_upper",
         floor = per10k_floor) %>%
  fix_ci("Nat_DeathsCaused_Resist_median", "Nat_DeathsCaused_Resist_lower", "Nat_DeathsCaused_Resist_upper",
         floor = count_floor) %>%
  fix_ci("Nat_DeathsCaused_Resist_per10k_median", "Nat_DeathsCaused_Resist_per10k_lower", "Nat_DeathsCaused_Resist_per10k_upper",
         floor = per10k_floor)

# ============================================================
# 2) PROVINCE + NATIONAL SUMMARIES (PBO PERIOD)
# ============================================================

# prov_pbo_totals_pw assumed to exist:
#  Region, ProvKey, Draw,
#  Cases_NoNets_PBOperiod, Cases_WithNets_Resist_PBOperiod,
#  Deaths_NoNets_PBOperiod, Deaths_WithNets_Resist_PBOperiod,
#  Cases_Averted_PBOperiod, LivesSaved_PBOperiod,
#  LivesSaved_PBOperiod_per10k

prov_pbo_summary <- prov_pbo_totals_pw %>%
  group_by(Region, ProvKey) %>%
  summarise(
    # PBO summaries use the raw/unscaled model outputs
    PercAverted_PBO = 100 *
      pmin(
        pmax(Cases_Averted_PBOperiod, 0) /
          pmax(Cases_NoNets_PBOperiod, eps),
        1
      ),

    LivesSaved_PBOperiod_abs    = pmax(LivesSaved_PBOperiod, 0),
    LivesSaved_PBOperiod_per10k = pmax(LivesSaved_PBOperiod_per10k, 0),
    .groups = "drop_last"
  ) %>%
  summarise(
    PercAverted_PBO_median = median(PercAverted_PBO, na.rm = TRUE),
    PercAverted_PBO_lower  = qL(PercAverted_PBO),
    PercAverted_PBO_upper  = qU(PercAverted_PBO),

    LivesSaved_PBOperiod_median =
      median(LivesSaved_PBOperiod_abs, na.rm = TRUE),
    LivesSaved_PBOperiod_lower  = qL(LivesSaved_PBOperiod_abs),
    LivesSaved_PBOperiod_upper  = qU(LivesSaved_PBOperiod_abs),

    LivesSaved_PBOperiod_per10k_median =
      median(LivesSaved_PBOperiod_per10k, na.rm = TRUE),
    LivesSaved_PBOperiod_per10k_lower  = qL(LivesSaved_PBOperiod_per10k),
    LivesSaved_PBOperiod_per10k_upper  = qU(LivesSaved_PBOperiod_per10k),

    .groups = "drop"
  ) %>%
  fix_ci("PercAverted_PBO_median", "PercAverted_PBO_lower", "PercAverted_PBO_upper",
         floor = pct_floor, cap = 100) %>%
  fix_ci("LivesSaved_PBOperiod_median", "LivesSaved_PBOperiod_lower", "LivesSaved_PBOperiod_upper",
         floor = count_floor) %>%
  fix_ci("LivesSaved_PBOperiod_per10k_median", "LivesSaved_PBOperiod_per10k_lower", "LivesSaved_PBOperiod_per10k_upper",
         floor = per10k_floor)

# nat_pbo_yearly_pw assumed to exist & nat_pop_pbo_mean from your PBO block
# Build national totals over PBO period
nat_pbo_total_summary <- nat_pbo_yearly_pw %>%
  group_by(Draw) %>%
  summarise(
    Nat_NoNets_PBO     = sum(Sum_NoNets_s,      na.rm = TRUE),
    Nat_Cases_Avt_PBO  = sum(pmax(Cases_Averted_PBO, 0), na.rm = TRUE),
    Nat_LivesSaved_PBO = sum(pmax(LivesSaved_PBO,    0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::crossing(nat_pop_pbo_mean) %>%
  mutate(
    Nat_PercAverted_PBO = 100 *
      pmin(Nat_Cases_Avt_PBO / pmax(Nat_NoNets_PBO, eps), 1),
    Nat_LivesSaved_PBO_per10k =
      Nat_LivesSaved_PBO / pmax(MeanNatPop_PBO / 10000, eps)
  ) %>%
  summarise(
    Nat_PercAverted_PBO_median = median(Nat_PercAverted_PBO, na.rm = TRUE),
    Nat_PercAverted_PBO_lower  = qL(Nat_PercAverted_PBO),
    Nat_PercAverted_PBO_upper  = qU(Nat_PercAverted_PBO),

    Nat_LivesSaved_PBO_median = median(Nat_LivesSaved_PBO, na.rm = TRUE),
    Nat_LivesSaved_PBO_lower  = qL(Nat_LivesSaved_PBO),
    Nat_LivesSaved_PBO_upper  = qU(Nat_LivesSaved_PBO),

    Nat_LivesSaved_PBO_per10k_median =
      median(Nat_LivesSaved_PBO_per10k, na.rm = TRUE),
    Nat_LivesSaved_PBO_per10k_lower  = qL(Nat_LivesSaved_PBO_per10k),
    Nat_LivesSaved_PBO_per10k_upper  = qU(Nat_LivesSaved_PBO_per10k),

    .groups = "drop"
  ) %>%
  fix_ci("Nat_PercAverted_PBO_median", "Nat_PercAverted_PBO_lower", "Nat_PercAverted_PBO_upper",
         floor = pct_floor, cap = 100) %>%
  fix_ci("Nat_LivesSaved_PBO_median", "Nat_LivesSaved_PBO_lower", "Nat_LivesSaved_PBO_upper",
         floor = count_floor) %>%
  fix_ci("Nat_LivesSaved_PBO_per10k_median", "Nat_LivesSaved_PBO_per10k_lower", "Nat_LivesSaved_PBO_per10k_upper",
         floor = per10k_floor)

# ------------------------------------------------------------
# Province ordering for the 26-province plots
# ------------------------------------------------------------
# ggplot draws the first factor level at the bottom of a horizontal plot.
# Therefore, using A-Z levels gives A-Z order from bottom to top.
province_levels_az <- stringi::stri_sort(
  unique(c(as.character(prov_summary$Region), as.character(prov_pbo_summary$Region))),
  locale = "fr"
)

prov_summary <- prov_summary %>%
  mutate(Region = factor(as.character(Region), levels = province_levels_az)) %>%
  arrange(Region)

prov_pbo_summary <- prov_pbo_summary %>%
  mutate(Region = factor(as.character(Region), levels = province_levels_az)) %>%
  arrange(Region)

# ============================================================
# 3) FOUR PROVINCES + DRC (MEAN) TABLES
# ============================================================

# ITNs / Resistance
four_plus_drc <- prov_summary %>%
  filter(Region %in% four_labels) %>%
  bind_rows(
    tibble(
      Region = "DRC (mean)",
      
      PercAverted_ITNs_median = nat_total_summary$PercAverted_ITNs_median,
      PercAverted_ITNs_lower  = nat_total_summary$PercAverted_ITNs_lower,
      PercAverted_ITNs_upper  = nat_total_summary$PercAverted_ITNs_upper,
      
      LivesSaved_ITNs_median  = nat_total_summary$Nat_LivesSaved_ITNs_median,
      LivesSaved_ITNs_lower   = nat_total_summary$Nat_LivesSaved_ITNs_lower,
      LivesSaved_ITNs_upper   = nat_total_summary$Nat_LivesSaved_ITNs_upper,
      
      LivesSaved_ITNs_per10k_median =
        nat_total_summary$Nat_LivesSaved_ITNs_per10k_median,
      LivesSaved_ITNs_per10k_lower  =
        nat_total_summary$Nat_LivesSaved_ITNs_per10k_lower,
      LivesSaved_ITNs_per10k_upper  =
        nat_total_summary$Nat_LivesSaved_ITNs_per10k_upper,
      
      PercIncrease_Resist_median =
        nat_total_summary$PercIncrease_Resist_median,
      PercIncrease_Resist_lower  =
        nat_total_summary$PercIncrease_Resist_lower,
      PercIncrease_Resist_upper  =
        nat_total_summary$PercIncrease_Resist_upper,
      
      DeathsCaused_Resist_median =
        nat_total_summary$Nat_DeathsCaused_Resist_median,
      DeathsCaused_Resist_lower  =
        nat_total_summary$Nat_DeathsCaused_Resist_lower,
      DeathsCaused_Resist_upper  =
        nat_total_summary$Nat_DeathsCaused_Resist_upper,
      
      DeathsCaused_Resist_per10k_median =
        nat_total_summary$Nat_DeathsCaused_Resist_per10k_median,
      DeathsCaused_Resist_per10k_lower  =
        nat_total_summary$Nat_DeathsCaused_Resist_per10k_lower,
      DeathsCaused_Resist_per10k_upper  =
        nat_total_summary$Nat_DeathsCaused_Resist_per10k_upper
    )
  ) %>%
  mutate(
    Region = factor(Region, levels = c(four_labels, "DRC (mean)"))
  )

# PBO
four_plus_drc_pbo <- prov_pbo_summary %>%
  filter(Region %in% four_labels) %>%
  bind_rows(
    tibble(
      Region = "DRC (mean)",
      
      PercAverted_PBO_median =
        nat_pbo_total_summary$Nat_PercAverted_PBO_median,
      PercAverted_PBO_lower  =
        nat_pbo_total_summary$Nat_PercAverted_PBO_lower,
      PercAverted_PBO_upper  =
        nat_pbo_total_summary$Nat_PercAverted_PBO_upper,
      
      LivesSaved_PBOperiod_median =
        nat_pbo_total_summary$Nat_LivesSaved_PBO_median,
      LivesSaved_PBOperiod_lower  =
        nat_pbo_total_summary$Nat_LivesSaved_PBO_lower,
      LivesSaved_PBOperiod_upper  =
        nat_pbo_total_summary$Nat_LivesSaved_PBO_upper,
      
      LivesSaved_PBOperiod_per10k_median =
        nat_pbo_total_summary$Nat_LivesSaved_PBO_per10k_median,
      LivesSaved_PBOperiod_per10k_lower  =
        nat_pbo_total_summary$Nat_LivesSaved_PBO_per10k_lower,
      LivesSaved_PBOperiod_per10k_upper  =
        nat_pbo_total_summary$Nat_LivesSaved_PBO_per10k_upper
    )
  ) %>%
  mutate(
    Region = factor(Region, levels = c(four_labels, "DRC (mean)"))
  )

# ============================================================
# 4) PROVINCES (26) – HORIZONTAL DOT + CI PLOTS
# ============================================================

# 4.A Cases averted due to ITNs (%)
p_prov_cases_averted <- ggplot(
  prov_summary,
  aes(x = PercAverted_ITNs_median,
      y = Region)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(
    aes(xmin = PercAverted_ITNs_lower,
        xmax = PercAverted_ITNs_upper),
    height = 0.25,
    alpha  = 0.8
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(
      x     = PercAverted_ITNs_median,
      label = fmt_pct1(PercAverted_ITNs_median)
    ),
    hjust = -0.2, size = 3, colour = "black"
  ) +
  scale_x_continuous("Cases averted due to ITNs (%)") +
  ylab("Province") +
  base_theme_box

# 4.B Lives saved due to ITNs – per 10 000, label = total
p_prov_lives_saved <- ggplot(
  prov_summary,
  aes(x = LivesSaved_ITNs_per10k_median,
      y = Region)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(
    aes(xmin = LivesSaved_ITNs_per10k_lower,
        xmax = LivesSaved_ITNs_per10k_upper),
    height = 0.25,
    alpha  = 0.8
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(
      x     = LivesSaved_ITNs_per10k_median,
      label = fmt_int(LivesSaved_ITNs_median)
    ),
    hjust = -0.2, size = 3, colour = "black"
  ) +
  scale_x_continuous(
    "Lives saved due to ITNs (per 10 000 population)",
    labels = label_number(accuracy = 0.1)
  ) +
  ylab("Province") +
  base_theme_box

# 4.C Case increase due to resistance (%)
p_prov_cases_increase_res <- ggplot(
  prov_summary,
  aes(x = PercIncrease_Resist_median,
      y = Region)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(
    aes(xmin = PercIncrease_Resist_lower,
        xmax = PercIncrease_Resist_upper),
    height = 0.25,
    alpha  = 0.8
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(
      x     = PercIncrease_Resist_median,
      label = fmt_pct1(PercIncrease_Resist_median)
    ),
    hjust = -0.2, size = 3, colour = "black"
  ) +
  scale_x_continuous("Case increase due to resistance (%)") +
  ylab("Province") +
  base_theme_box

# 4.D Deaths caused by resistance – per 10 000, label = total
p_prov_deaths_res <- ggplot(
  prov_summary,
  aes(x = DeathsCaused_Resist_per10k_median,
      y = Region)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(
    aes(
      xmin = pmax(DeathsCaused_Resist_per10k_lower, 0),
      xmax = pmax(DeathsCaused_Resist_per10k_upper, 0)
    ),
    height = 0.25,
    alpha  = 0.8
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(
      x     = pmax(DeathsCaused_Resist_per10k_median, 0),
      label = fmt_int(DeathsCaused_Resist_median)
    ),
    hjust = -0.2, size = 3, colour = "black"
  ) +
  scale_x_continuous(
    "Deaths caused by resistance (per 10 000 population)",
    labels = label_number(
      accuracy  = 0.1,
      scale_cut = cut_short_scale()
    )
  ) +
  ylab("Province") +
  base_theme_box

# 4.E Cases averted due to PBO (%)
p_prov_cases_averted_pbo <- ggplot(
  prov_pbo_summary,
  aes(x = PercAverted_PBO_median,
      y = Region)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(
    aes(xmin = PercAverted_PBO_lower,
        xmax = PercAverted_PBO_upper),
    height = 0.25,
    alpha  = 0.8
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(
      x     = PercAverted_PBO_median,
      label = fmt_pct1(PercAverted_PBO_median)
    ),
    hjust = -0.2, size = 3, colour = "black"
  ) +
  scale_x_continuous("Cases averted due to PBO (%) ") +
  ylab("Province") +
  base_theme_box

# 4.F Lives saved due to PBO – per 10 000, label = total
p_prov_lives_saved_pbo <- ggplot(
  prov_pbo_summary,
  aes(x = LivesSaved_PBOperiod_per10k_median,
      y = Region)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(
    aes(xmin = LivesSaved_PBOperiod_per10k_lower,
        xmax = LivesSaved_PBOperiod_per10k_upper),
    height = 0.25,
    alpha  = 0.8
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(
      x     = LivesSaved_PBOperiod_per10k_median,
      label = fmt_int(LivesSaved_PBOperiod_median)
    ),
    hjust = -0.2, size = 3, colour = "black"
  ) +
  scale_x_continuous(
    "Lives saved due to PBO (per 10 000 population)",
    labels = label_number(
      accuracy  = 0.1,
      scale_cut = cut_short_scale()
    )
  ) +
  ylab("Province") +
  base_theme_box

# ============================================================
# 5) FOUR REGIONS + DRC (MEAN)
# ============================================================

# 5.A Cases averted due to ITNs (%)
p_four_cases_averted <- ggplot(
  four_plus_drc,
  aes(x = PercAverted_ITNs_median, y = Region, fill = Region)
) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_errorbarh(
    aes(xmin = PercAverted_ITNs_lower,
        xmax = PercAverted_ITNs_upper),
    height = 0.2
  ) +
  geom_text(
    aes(
      x     = PercAverted_ITNs_median,
      label = fmt_pct1(PercAverted_ITNs_median)
    ),
    hjust = -0.2, size = 3, colour = "black"
  ) +
  scale_x_continuous("Cases averted due to ITNs (%)") +
  scale_fill_manual(values = col_four_regions_fill) +
  ylab("Province") +
  base_theme_box +
  theme(legend.position = "none")

# 5.B Lives saved due to ITNs – per 10 000, label = total
p_four_lives_saved <- ggplot(
  four_plus_drc,
  aes(x = LivesSaved_ITNs_per10k_median, y = Region, colour = Region)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(
    aes(xmin = LivesSaved_ITNs_per10k_lower,
        xmax = LivesSaved_ITNs_per10k_upper),
    height = 0.2
  ) +
  geom_point(size = 3) +
  geom_text(
    aes(
      x     = LivesSaved_ITNs_per10k_median,
      label = fmt_int(LivesSaved_ITNs_median)
    ),
    hjust = -0.2, size = 3, colour = "black"
  ) +
  scale_x_continuous(
    "Lives saved due to ITNs (per 10 000 population)",
    labels = label_number(
      accuracy  = 0.1,
      scale_cut = cut_short_scale()
    )
  ) +
  scale_colour_manual(values = col_four_regions_point) +
  ylab("Province") +
  base_theme_box +
  theme(legend.position = "none")

# 5.C Case increase due to resistance (%)
p_four_cases_increase_res <- ggplot(
  four_plus_drc,
  aes(x = PercIncrease_Resist_median, y = Region, fill = Region)
) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_errorbarh(
    aes(xmin = PercIncrease_Resist_lower,
        xmax = PercIncrease_Resist_upper),
    height = 0.2
  ) +
  geom_text(
    aes(
      x     = PercIncrease_Resist_median,
      label = fmt_pct1(PercIncrease_Resist_median)
    ),
    hjust = -0.2, size = 3, colour = "black"
  ) +
  scale_x_continuous("Case increase due to resistance (%)") +
  scale_fill_manual(values = col_four_regions_fill) +
  ylab("Province") +
  base_theme_box +
  theme(legend.position = "none")

# 5.D Deaths caused by resistance – per 10 000, label = total
p_four_deaths_res <- ggplot(
  four_plus_drc,
  aes(x = DeathsCaused_Resist_per10k_median, y = Region, colour = Region)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(
    aes(
      xmin = pmax(DeathsCaused_Resist_per10k_lower, 0),
      xmax = pmax(DeathsCaused_Resist_per10k_upper, 0)
    ),
    height = 0.2
  ) +
  geom_point(size = 3) +
  geom_text(
    aes(
      x     = pmax(DeathsCaused_Resist_per10k_median, 0),
      label = fmt_int(DeathsCaused_Resist_median)
    ),
    hjust = -0.2, size = 3, colour = "black"
  ) +
  scale_x_continuous(
    "Deaths caused by resistance (per 10 000 population)",
    labels = label_number(
      accuracy  = 0.1,
      scale_cut = cut_short_scale()
    )
  ) +
  scale_colour_manual(values = col_four_regions_point) +
  ylab("Province") +
  base_theme_box +
  theme(legend.position = "none")

# 5.E Cases averted due to PBO (%)
p_four_cases_averted_pbo <- ggplot(
  four_plus_drc_pbo,
  aes(x = PercAverted_PBO_median, y = Region, fill = Region)
) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_errorbarh(
    aes(xmin = PercAverted_PBO_lower,
        xmax = PercAverted_PBO_upper),
    height = 0.2
  ) +
  geom_text(
    aes(
      x     = PercAverted_PBO_median,
      label = fmt_pct1(PercAverted_PBO_median)
    ),
    hjust = -0.2, size = 3, colour = "black"
  ) +
  scale_x_continuous("Cases averted due to PBO (%) ") +
  scale_fill_manual(values = col_four_regions_fill) +
  ylab("Province") +
  base_theme_box +
  theme(legend.position = "none")

# 5.F Lives saved due to PBO – per 10 000, label = total
p_four_lives_saved_pbo <- ggplot(
  four_plus_drc_pbo,
  aes(x = LivesSaved_PBOperiod_per10k_median, y = Region, colour = Region)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(
    aes(xmin = LivesSaved_PBOperiod_per10k_lower,
        xmax = LivesSaved_PBOperiod_per10k_upper),
    height = 0.2
  ) +
  geom_point(size = 3) +
  geom_text(
    aes(
      x     = LivesSaved_PBOperiod_per10k_median,
      label = fmt_int(LivesSaved_PBOperiod_median)
    ),
    hjust = -0.2, size = 3, colour = "black"
  ) +
  scale_x_continuous(
    "Lives saved due to PBO (per 10 000 population)",
    labels = label_number(
      accuracy  = 0.1,
      scale_cut = cut_short_scale()
    )
  ) +
  scale_colour_manual(values = col_four_regions_point) +
  ylab("Province") +
  base_theme_box +
  theme(legend.position = "none")

# ============================================================
# 6) NATIONAL OVER TIME (2009–2024) – BLUE BARS
# ============================================================

# national_yearly_pw assumed to exist with Sum_* columns

nat_yearly_summary_raw <- national_yearly_pw %>%
  filter(Year >= 2009, Year <= 2024) %>%
  mutate(
    # ITN yearly outputs use raw/unscaled model values
    Cases_Averted_ITNs =
      pmax(Sum_NoNets_s - Sum_WithNets_Resist_s, 0),
    LivesSaved_ITNs =
      pmax(Sum_Deaths_NoNets_s - Sum_Deaths_WithNets_Resist_s, 0),
    PercAverted_ITNs =
      100 * pmin(Cases_Averted_ITNs / pmax(Sum_NoNets_s, eps), 1),

    # Resistance-attributable yearly outputs use WHO-calibrated values only
    DeathsCaused_Resist =
      pmax(Sum_Deaths_WithNets_Resist_who - Sum_Deaths_WithNets_NoResist_who, 0),
    PercIncrease_Resist =
      100 * pmax(Sum_WithNets_Resist_who - Sum_WithNets_NoResist_who, 0) /
      pmax(Sum_WithNets_NoResist_who, eps)
  ) %>%
  group_by(Year) %>%
  summarise(
    PercAverted_ITNs_median = median(PercAverted_ITNs, na.rm = TRUE),
    PercAverted_ITNs_lower  = qL(PercAverted_ITNs),
    PercAverted_ITNs_upper  = qU(PercAverted_ITNs),

    LivesSaved_ITNs_median  = median(LivesSaved_ITNs, na.rm = TRUE),
    LivesSaved_ITNs_lower   = qL(LivesSaved_ITNs),
    LivesSaved_ITNs_upper   = qU(LivesSaved_ITNs),

    PercIncrease_Resist_median = median(PercIncrease_Resist, na.rm = TRUE),
    PercIncrease_Resist_lower  = qL(PercIncrease_Resist),
    PercIncrease_Resist_upper  = qU(PercIncrease_Resist),

    DeathsCaused_Resist_median = median(DeathsCaused_Resist, na.rm = TRUE),
    DeathsCaused_Resist_lower  = qL(DeathsCaused_Resist),
    DeathsCaused_Resist_upper  = qU(DeathsCaused_Resist),

    .groups = "drop"
  ) %>%
  fix_ci("PercAverted_ITNs_median", "PercAverted_ITNs_lower", "PercAverted_ITNs_upper",
         floor = pct_floor, cap = 100) %>%
  fix_ci("LivesSaved_ITNs_median", "LivesSaved_ITNs_lower", "LivesSaved_ITNs_upper",
         floor = count_floor) %>%
  fix_ci("PercIncrease_Resist_median", "PercIncrease_Resist_lower", "PercIncrease_Resist_upper",
         floor = pct_floor) %>%
  fix_ci("DeathsCaused_Resist_median", "DeathsCaused_Resist_lower", "DeathsCaused_Resist_upper",
         floor = count_floor)

full_years_nat <- tibble(Year = 2009:2024)

nat_yearly_summary <- full_years_nat %>%
  left_join(nat_yearly_summary_raw, by = "Year") %>%
  mutate(
    across(where(is.numeric), ~replace_na(.x, 0))
  ) %>%
  fix_ci("PercAverted_ITNs_median", "PercAverted_ITNs_lower", "PercAverted_ITNs_upper",
         floor = pct_floor, cap = 100) %>%
  fix_ci("LivesSaved_ITNs_median", "LivesSaved_ITNs_lower", "LivesSaved_ITNs_upper",
         floor = count_floor) %>%
  fix_ci("PercIncrease_Resist_median", "PercIncrease_Resist_lower", "PercIncrease_Resist_upper",
         floor = pct_floor) %>%
  fix_ci("DeathsCaused_Resist_median", "DeathsCaused_Resist_lower", "DeathsCaused_Resist_upper",
         floor = count_floor)

first_nonzero_year <- nat_yearly_summary %>%
  filter(
    LivesSaved_ITNs_median > 0 |
      PercAverted_ITNs_median > 0 |
      DeathsCaused_Resist_median > 0 |
      PercIncrease_Resist_median > 0
  ) %>%
  summarise(min_year = if_else(n() > 0, min(Year), 2009L)) %>%
  pull(min_year)

if (length(first_nonzero_year) == 0 || is.na(first_nonzero_year)) {
  first_nonzero_year <- 2009L
}

nat_yearly_summary <- nat_yearly_summary %>%
  mutate(
    PercAverted_ITNs_plot = if_else(
      Year < first_nonzero_year & PercAverted_ITNs_median == 0,
      eps_bar,
      PercAverted_ITNs_median
    ),
    LivesSaved_ITNs_plot = if_else(
      Year < first_nonzero_year & LivesSaved_ITNs_median == 0,
      eps_bar,
      LivesSaved_ITNs_median
    ),
    PercIncrease_Resist_plot = if_else(
      Year < first_nonzero_year & PercIncrease_Resist_median == 0,
      eps_bar,
      PercIncrease_Resist_median
    ),
    DeathsCaused_Resist_plot = if_else(
      Year < first_nonzero_year & DeathsCaused_Resist_median == 0,
      eps_bar,
      DeathsCaused_Resist_median
    )
  )

# nat_pbo_yearly_summary for PBO (2019–2024)
nat_pbo_yearly_summary_raw <- nat_pbo_yearly_pw %>%
  filter(Year >= 2019, Year <= 2024) %>%
  mutate(
    Cases_Averted_PBO = pmax(Cases_Averted_PBO, 0),
    LivesSaved_PBO    = pmax(LivesSaved_PBO, 0)
  ) %>%
  group_by(Year) %>%
  summarise(
    Nat_Cases_Averted_PBO_median = median(Cases_Averted_PBO, na.rm = TRUE),
    Nat_Cases_Averted_PBO_lower  = qL(Cases_Averted_PBO),
    Nat_Cases_Averted_PBO_upper  = qU(Cases_Averted_PBO),

    Nat_LivesSaved_PBO_median    = median(LivesSaved_PBO, na.rm = TRUE),
    Nat_LivesSaved_PBO_lower     = qL(LivesSaved_PBO),
    Nat_LivesSaved_PBO_upper     = qU(LivesSaved_PBO),

    .groups = "drop"
  ) %>%
  fix_ci("Nat_Cases_Averted_PBO_median", "Nat_Cases_Averted_PBO_lower", "Nat_Cases_Averted_PBO_upper",
         floor = count_floor) %>%
  fix_ci("Nat_LivesSaved_PBO_median", "Nat_LivesSaved_PBO_lower", "Nat_LivesSaved_PBO_upper",
         floor = count_floor)

full_years_pbo <- tibble(Year = 2019:2024)

nat_pbo_yearly_summary <- full_years_pbo %>%
  left_join(nat_pbo_yearly_summary_raw, by = "Year") %>%
  mutate(
    across(where(is.numeric), ~replace_na(.x, 0))
  ) %>%
  fix_ci("Nat_Cases_Averted_PBO_median", "Nat_Cases_Averted_PBO_lower", "Nat_Cases_Averted_PBO_upper",
         floor = count_floor) %>%
  fix_ci("Nat_LivesSaved_PBO_median", "Nat_LivesSaved_PBO_lower", "Nat_LivesSaved_PBO_upper",
         floor = count_floor)

# 6.A
p_ITN_nat_time <- ggplot(
  nat_yearly_summary,
  aes(x = Year, y = PercAverted_ITNs_plot)
) +
  geom_col(fill = col_nat_overtime, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = PercAverted_ITNs_lower,
        ymax = PercAverted_ITNs_upper),
    width = 0.25
  ) +
  ylab("DRC Cases averted due to ITNs (%)") +
  xlab("Year") +
  base_theme_box +
  scale_x_continuous(breaks = 2009:2024) +
  theme_x_oblique

# 6.B
p_LS_nat_time <- ggplot(
  nat_yearly_summary,
  aes(x = Year, y = LivesSaved_ITNs_plot)
) +
  geom_col(fill = col_nat_overtime, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = LivesSaved_ITNs_lower,
        ymax = LivesSaved_ITNs_upper),
    width = 0.25
  ) +
  ylab("DRC Lives saved due to ITNs") +
  xlab("Year") +
  base_theme_box +
  scale_x_continuous(breaks = 2009:2024) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_x_oblique

# 6.C
p_RES_nat_time <- ggplot(
  nat_yearly_summary,
  aes(x = Year, y = PercIncrease_Resist_plot)
) +
  geom_col(fill = col_nat_overtime, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = PercIncrease_Resist_lower,
        ymax = PercIncrease_Resist_upper),
    width = 0.25
  ) +
  ylab("DRC Case increase due to resistance (%)") +
  xlab("Year") +
  base_theme_box +
  scale_x_continuous(breaks = 2009:2024) +
  theme_x_oblique

# 6.D
p_DEAD_nat_time <- ggplot(
  nat_yearly_summary,
  aes(x = Year, y = DeathsCaused_Resist_plot)
) +
  geom_col(fill = col_nat_overtime, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = DeathsCaused_Resist_lower,
        ymax = DeathsCaused_Resist_upper),
    width = 0.25
  ) +
  ylab("DRC Deaths caused by resistance") +
  xlab("Year") +
  base_theme_box +
  scale_x_continuous(breaks = 2009:2024) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_x_oblique

# 6.E
p_PBO_nat_time <- ggplot(
  nat_pbo_yearly_summary,
  aes(x = Year, y = Nat_Cases_Averted_PBO_median)
) +
  geom_col(fill = col_nat_overtime, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = Nat_Cases_Averted_PBO_lower,
        ymax = Nat_Cases_Averted_PBO_upper),
    width = 0.25
  ) +
  ylab("DRC Cases averted due to PBO") +
  xlab("Year") +
  base_theme_box +
  scale_x_continuous(breaks = 2019:2024) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_x_oblique

# 6.F
p_PBO_LS_nat_time <- ggplot(
  nat_pbo_yearly_summary,
  aes(x = Year, y = Nat_LivesSaved_PBO_median)
) +
  geom_col(fill = col_nat_overtime, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = Nat_LivesSaved_PBO_lower,
        ymax = Nat_LivesSaved_PBO_upper),
    width = 0.25
  ) +
  ylab("DRC Lives saved due to PBO") +
  xlab("Year") +
  base_theme_box +
  scale_x_continuous(breaks = 2019:2024) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_x_oblique

# ============================================================
# 7) ALIASES + PRINT
# ============================================================

p_ITN_prov      <- p_prov_cases_averted
p_ITN_five      <- p_four_cases_averted

p_LS_prov       <- p_prov_lives_saved
p_LS_five       <- p_four_lives_saved

p_RES_prov      <- p_prov_cases_increase_res
p_RES_five      <- p_four_cases_increase_res

p_DEAD_prov     <- p_prov_deaths_res
p_DEAD_five     <- p_four_deaths_res

p_PBO_prov      <- p_prov_cases_averted_pbo
p_PBO_five      <- p_four_cases_averted_pbo

p_PBO_LS_prov   <- p_prov_lives_saved_pbo
p_PBO_LS_five   <- p_four_lives_saved_pbo

p_ITN_nat       <- p_ITN_nat_time
p_LS_nat        <- p_LS_nat_time
p_RES_nat       <- p_RES_nat_time
p_DEAD_nat      <- p_DEAD_nat_time
p_PBO_nat       <- p_PBO_nat_time
p_PBO_LS_nat    <- p_PBO_LS_nat_time


# ============================================================
# 8) PRINTABLE TABLES FOR EACH PLOT
# ============================================================
# These tables give the numerical values behind the plots.
# For each output, the table reports the posterior mean estimate and the
# 95% credible interval (2.5% and 97.5% posterior quantiles).
# ITN and PBO tables use raw population-weighted model outputs.
# Resistance-attributable case/death tables use WHO-calibrated outputs only.

fmt_num_table <- function(x, digits = 0) {
  ifelse(
    is.na(x),
    NA_character_,
    formatC(round(x, digits), format = "f", digits = digits, big.mark = " ")
  )
}

fmt_ci_table <- function(est, lo, up, digits = 0, suffix = "") {
  ifelse(
    is.na(est) | is.na(lo) | is.na(up),
    NA_character_,
    paste0(
      fmt_num_table(est, digits), suffix,
      " (", fmt_num_table(lo, digits), "–", fmt_num_table(up, digits), suffix, ")"
    )
  )
}

repair_interval_values <- function(est, lo, up, floor = 0, cap = Inf) {
  est <- pmin(pmax(as.numeric(est), floor), cap)
  lo  <- pmin(pmax(as.numeric(lo),  floor), cap)
  up  <- pmin(pmax(as.numeric(up),  floor), cap)
  lo  <- pmin(lo, est)
  up  <- pmax(up, est)
  list(est = est, lo = lo, up = up)
}

make_abs_pct_table <- function(data, group_var, group_label,
                               abs_col, pct_col,
                               abs_name, pct_name,
                               abs_digits = 0,
                               pct_digits = 1,
                               pct_cap = Inf) {
  out <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_var))) %>%
    dplyr::summarise(
      abs_mean  = mean(.data[[abs_col]], na.rm = TRUE),
      abs_lower = qL(.data[[abs_col]]),
      abs_upper = qU(.data[[abs_col]]),
      pct_mean  = mean(.data[[pct_col]], na.rm = TRUE),
      pct_lower = qL(.data[[pct_col]]),
      pct_upper = qU(.data[[pct_col]]),
      .groups = "drop"
    )

  names(out)[names(out) == group_var] <- group_label

  abs_fixed <- repair_interval_values(out$abs_mean, out$abs_lower, out$abs_upper,
                                      floor = count_floor, cap = Inf)
  pct_fixed <- repair_interval_values(out$pct_mean, out$pct_lower, out$pct_upper,
                                      floor = pct_floor, cap = pct_cap)

  out[[paste0(abs_name, "_mean")]]      <- round(abs_fixed$est, abs_digits)
  out[[paste0(abs_name, "_lower95")]]   <- round(abs_fixed$lo,  abs_digits)
  out[[paste0(abs_name, "_upper95")]]   <- round(abs_fixed$up,  abs_digits)
  out[[paste0(abs_name, "_95CrI")]]     <- fmt_ci_table(abs_fixed$est, abs_fixed$lo, abs_fixed$up,
                                                        digits = abs_digits)

  out[[paste0(pct_name, "_mean")]]      <- round(pct_fixed$est, pct_digits)
  out[[paste0(pct_name, "_lower95")]]   <- round(pct_fixed$lo,  pct_digits)
  out[[paste0(pct_name, "_upper95")]]   <- round(pct_fixed$up,  pct_digits)
  out[[paste0(pct_name, "_95CrI")]]     <- fmt_ci_table(pct_fixed$est, pct_fixed$lo, pct_fixed$up,
                                                        digits = pct_digits, suffix = "%")

  out %>%
    dplyr::select(-abs_mean, -abs_lower, -abs_upper,
                  -pct_mean, -pct_lower, -pct_upper)
}

make_abs_per10k_table <- function(data, group_var, group_label,
                                  abs_col, per10k_col,
                                  abs_name, per10k_name,
                                  abs_digits = 0,
                                  per10k_digits = 1) {
  out <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_var))) %>%
    dplyr::summarise(
      abs_mean     = mean(.data[[abs_col]], na.rm = TRUE),
      abs_lower    = qL(.data[[abs_col]]),
      abs_upper    = qU(.data[[abs_col]]),
      per10k_mean  = mean(.data[[per10k_col]], na.rm = TRUE),
      per10k_lower = qL(.data[[per10k_col]]),
      per10k_upper = qU(.data[[per10k_col]]),
      .groups = "drop"
    )

  names(out)[names(out) == group_var] <- group_label

  abs_fixed <- repair_interval_values(out$abs_mean, out$abs_lower, out$abs_upper,
                                      floor = count_floor, cap = Inf)
  per10k_fixed <- repair_interval_values(out$per10k_mean, out$per10k_lower, out$per10k_upper,
                                         floor = per10k_floor, cap = Inf)

  out[[paste0(abs_name, "_mean")]]      <- round(abs_fixed$est, abs_digits)
  out[[paste0(abs_name, "_lower95")]]   <- round(abs_fixed$lo,  abs_digits)
  out[[paste0(abs_name, "_upper95")]]   <- round(abs_fixed$up,  abs_digits)
  out[[paste0(abs_name, "_95CrI")]]     <- fmt_ci_table(abs_fixed$est, abs_fixed$lo, abs_fixed$up,
                                                        digits = abs_digits)

  out[[paste0(per10k_name, "_mean")]]      <- round(per10k_fixed$est, per10k_digits)
  out[[paste0(per10k_name, "_lower95")]]   <- round(per10k_fixed$lo,  per10k_digits)
  out[[paste0(per10k_name, "_upper95")]]   <- round(per10k_fixed$up,  per10k_digits)
  out[[paste0(per10k_name, "_95CrI")]]     <- fmt_ci_table(per10k_fixed$est, per10k_fixed$lo, per10k_fixed$up,
                                                           digits = per10k_digits)

  out %>%
    dplyr::select(-abs_mean, -abs_lower, -abs_upper,
                  -per10k_mean, -per10k_lower, -per10k_upper)
}

order_province_table <- function(tab, label_col = "Province") {
  tab %>%
    dplyr::mutate(
      !!label_col := as.character(.data[[label_col]]),
      .province_order = factor(.data[[label_col]], levels = c(as.character(province_levels_az), "DRC (mean)"))
    ) %>%
    dplyr::arrange(.province_order) %>%
    dplyr::select(-.province_order)
}

order_five_table <- function(tab, label_col = "Province") {
  tab %>%
    dplyr::filter(as.character(.data[[label_col]]) %in% c(four_labels, "DRC (mean)")) %>%
    dplyr::mutate(
      !!label_col := as.character(.data[[label_col]]),
      .province_order = factor(.data[[label_col]], levels = c(four_labels, "DRC (mean)"))
    ) %>%
    dplyr::arrange(.province_order) %>%
    dplyr::select(-.province_order)
}

# Province-level draw data plus one DRC aggregate row.
prov_metric_draws <- prov_totals_pw %>%
  dplyr::mutate(
    Region = as.character(Region),
    Cases_Averted_ITNs_abs = pmax(Sum_NoNets_s - Sum_WithNets_Resist_s, 0),
    PercAverted_ITNs = 100 * pmin(Cases_Averted_ITNs_abs / pmax(Sum_NoNets_s, eps), 1),

    LivesSaved_ITNs_abs = pmax(LivesSaved_ITN_abs, 0),
    LivesSaved_ITNs_per10k = pmax(LivesSaved_ITN_per10k, 0),

    CasesCaused_Resist_abs = pmax(CasesCaused_Resist_abs, 0),
    PercIncrease_Resist = 100 *
      pmax(Sum_WithNets_Resist_who - Sum_WithNets_NoResist_who, 0) /
      pmax(Sum_WithNets_NoResist_who, eps),

    DeathsCaused_Resist_abs = pmax(DeathsCaused_Resist_abs, 0),
    DeathsCaused_Resist_per10k = pmax(DeathsCaused_Resist_per10k, 0)
  )

drc_metric_draws <- national_total_pw %>%
  dplyr::transmute(
    Region = "DRC (mean)",
    ProvKey = "DRC_mean",
    Draw,

    Cases_Averted_ITNs_abs = pmax(Nat_Sum_NoNets_s - Nat_Sum_WithNets_Resist_s, 0),
    PercAverted_ITNs = 100 *
      pmin(Cases_Averted_ITNs_abs / pmax(Nat_Sum_NoNets_s, eps), 1),

    LivesSaved_ITNs_abs = pmax(Nat_LivesSaved_ITN_abs, 0),
    LivesSaved_ITNs_per10k = pmax(Nat_LivesSaved_ITN_per10k, 0),

    CasesCaused_Resist_abs = pmax(Nat_CasesCaused_Resist_abs, 0),
    PercIncrease_Resist = 100 *
      pmax(Nat_Sum_WithNets_Resist_who - Nat_Sum_WithNets_NoResist_who, 0) /
      pmax(Nat_Sum_WithNets_NoResist_who, eps),

    DeathsCaused_Resist_abs = pmax(Nat_DeathsCaused_Resist_abs, 0),
    DeathsCaused_Resist_per10k = pmax(Nat_DeathsCaused_Resist_per10k, 0)
  )

prov_drc_metric_draws <- dplyr::bind_rows(
  prov_metric_draws %>%
    dplyr::select(Region, ProvKey, Draw,
                  Cases_Averted_ITNs_abs, PercAverted_ITNs,
                  LivesSaved_ITNs_abs, LivesSaved_ITNs_per10k,
                  CasesCaused_Resist_abs, PercIncrease_Resist,
                  DeathsCaused_Resist_abs, DeathsCaused_Resist_per10k),
  drc_metric_draws
)

# PBO draw data plus one DRC aggregate row.
prov_pbo_metric_draws <- prov_pbo_totals_pw %>%
  dplyr::mutate(
    Region = as.character(Region),
    Cases_Averted_PBO_abs = pmax(Cases_Averted_PBOperiod, 0),
    PercAverted_PBO = 100 *
      pmin(Cases_Averted_PBO_abs / pmax(Cases_NoNets_PBOperiod, eps), 1),
    LivesSaved_PBO_abs = pmax(LivesSaved_PBOperiod, 0),
    LivesSaved_PBO_per10k = pmax(LivesSaved_PBOperiod_per10k, 0)
  )

drc_pbo_metric_draws <- nat_pbo_yearly_pw %>%
  dplyr::group_by(Draw) %>%
  dplyr::summarise(
    Nat_NoNets_PBO = sum(Sum_NoNets_s, na.rm = TRUE),
    Cases_Averted_PBO_abs = sum(pmax(Cases_Averted_PBO, 0), na.rm = TRUE),
    LivesSaved_PBO_abs = sum(pmax(LivesSaved_PBO, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::crossing(nat_pop_pbo_mean) %>%
  dplyr::transmute(
    Region = "DRC (mean)",
    ProvKey = "DRC_mean",
    Draw,
    Cases_Averted_PBO_abs,
    PercAverted_PBO = 100 *
      pmin(Cases_Averted_PBO_abs / pmax(Nat_NoNets_PBO, eps), 1),
    LivesSaved_PBO_abs,
    LivesSaved_PBO_per10k = LivesSaved_PBO_abs / pmax(MeanNatPop_PBO / 10000, eps)
  )

prov_drc_pbo_metric_draws <- dplyr::bind_rows(
  prov_pbo_metric_draws %>%
    dplyr::select(Region, ProvKey, Draw,
                  Cases_Averted_PBO_abs, PercAverted_PBO,
                  LivesSaved_PBO_abs, LivesSaved_PBO_per10k),
  drc_pbo_metric_draws
)

# Province + DRC tables: one table for each province plot.
table_ITN_prov <- make_abs_pct_table(
  prov_drc_metric_draws, "Region", "Province",
  "Cases_Averted_ITNs_abs", "PercAverted_ITNs",
  "Cases_averted_ITNs", "Percent_cases_averted_ITNs",
  pct_cap = 100
) %>% order_province_table("Province")

table_LS_prov <- make_abs_per10k_table(
  prov_drc_metric_draws, "Region", "Province",
  "LivesSaved_ITNs_abs", "LivesSaved_ITNs_per10k",
  "Lives_saved_ITNs", "Lives_saved_ITNs_per10000"
) %>% order_province_table("Province")

table_RES_prov <- make_abs_pct_table(
  prov_drc_metric_draws, "Region", "Province",
  "CasesCaused_Resist_abs", "PercIncrease_Resist",
  "Additional_cases_due_resistance", "Percent_case_increase_resistance",
  pct_cap = Inf
) %>% order_province_table("Province")

table_DEAD_prov <- make_abs_per10k_table(
  prov_drc_metric_draws, "Region", "Province",
  "DeathsCaused_Resist_abs", "DeathsCaused_Resist_per10k",
  "Deaths_due_resistance", "Deaths_due_resistance_per10000"
) %>% order_province_table("Province")

table_PBO_prov <- make_abs_pct_table(
  prov_drc_pbo_metric_draws, "Region", "Province",
  "Cases_Averted_PBO_abs", "PercAverted_PBO",
  "Cases_averted_PBO", "Percent_cases_averted_PBO",
  pct_cap = 100
) %>% order_province_table("Province")

table_PBO_LS_prov <- make_abs_per10k_table(
  prov_drc_pbo_metric_draws, "Region", "Province",
  "LivesSaved_PBO_abs", "LivesSaved_PBO_per10k",
  "Lives_saved_PBO", "Lives_saved_PBO_per10000"
) %>% order_province_table("Province")

# Four selected provinces + DRC mean tables.
table_ITN_five    <- order_five_table(table_ITN_prov, "Province")
table_LS_five     <- order_five_table(table_LS_prov, "Province")
table_RES_five    <- order_five_table(table_RES_prov, "Province")
table_DEAD_five   <- order_five_table(table_DEAD_prov, "Province")
table_PBO_five    <- order_five_table(table_PBO_prov, "Province")
table_PBO_LS_five <- order_five_table(table_PBO_LS_prov, "Province")

# National yearly tables: one table for each national-over-time plot.
nat_yearly_metric_draws <- national_yearly_pw %>%
  dplyr::filter(Year >= 2009, Year <= 2024) %>%
  dplyr::mutate(
    Cases_Averted_ITNs = pmax(Sum_NoNets_s - Sum_WithNets_Resist_s, 0),
    PercAverted_ITNs = 100 * pmin(Cases_Averted_ITNs / pmax(Sum_NoNets_s, eps), 1),

    LivesSaved_ITNs = pmax(Sum_Deaths_NoNets_s - Sum_Deaths_WithNets_Resist_s, 0),

    CasesCaused_Resist = pmax(Sum_WithNets_Resist_who - Sum_WithNets_NoResist_who, 0),
    PercIncrease_Resist = 100 * CasesCaused_Resist / pmax(Sum_WithNets_NoResist_who, eps),

    DeathsCaused_Resist = pmax(Sum_Deaths_WithNets_Resist_who - Sum_Deaths_WithNets_NoResist_who, 0)
  )

table_ITN_nat <- make_abs_pct_table(
  nat_yearly_metric_draws, "Year", "Year",
  "Cases_Averted_ITNs", "PercAverted_ITNs",
  "Cases_averted_ITNs", "Percent_cases_averted_ITNs",
  pct_cap = 100
) %>% dplyr::arrange(Year)

table_LS_nat <- nat_yearly_metric_draws %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(
    Lives_saved_ITNs_mean = round(mean(LivesSaved_ITNs, na.rm = TRUE), 0),
    Lives_saved_ITNs_lower95 = round(qL(LivesSaved_ITNs), 0),
    Lives_saved_ITNs_upper95 = round(qU(LivesSaved_ITNs), 0),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Lives_saved_ITNs_95CrI = fmt_ci_table(
      Lives_saved_ITNs_mean,
      Lives_saved_ITNs_lower95,
      Lives_saved_ITNs_upper95,
      digits = 0
    )
  ) %>%
  dplyr::arrange(Year)

table_RES_nat <- make_abs_pct_table(
  nat_yearly_metric_draws, "Year", "Year",
  "CasesCaused_Resist", "PercIncrease_Resist",
  "Additional_cases_due_resistance", "Percent_case_increase_resistance",
  pct_cap = Inf
) %>% dplyr::arrange(Year)

table_DEAD_nat <- nat_yearly_metric_draws %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(
    Deaths_due_resistance_mean = round(mean(DeathsCaused_Resist, na.rm = TRUE), 0),
    Deaths_due_resistance_lower95 = round(qL(DeathsCaused_Resist), 0),
    Deaths_due_resistance_upper95 = round(qU(DeathsCaused_Resist), 0),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Deaths_due_resistance_95CrI = fmt_ci_table(
      Deaths_due_resistance_mean,
      Deaths_due_resistance_lower95,
      Deaths_due_resistance_upper95,
      digits = 0
    )
  ) %>%
  dplyr::arrange(Year)

table_PBO_nat <- nat_pbo_yearly_pw %>%
  dplyr::filter(Year >= 2019, Year <= 2024) %>%
  dplyr::mutate(
    Cases_Averted_PBO = pmax(Cases_Averted_PBO, 0),
    PercAverted_PBO = 100 * pmin(Cases_Averted_PBO / pmax(Sum_NoNets_s, eps), 1)
  ) %>%
  make_abs_pct_table(
    "Year", "Year",
    "Cases_Averted_PBO", "PercAverted_PBO",
    "Cases_averted_PBO", "Percent_cases_averted_PBO",
    pct_cap = 100
  ) %>%
  dplyr::arrange(Year)

table_PBO_LS_nat <- nat_pbo_yearly_pw %>%
  dplyr::filter(Year >= 2019, Year <= 2024) %>%
  dplyr::mutate(LivesSaved_PBO = pmax(LivesSaved_PBO, 0)) %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(
    Lives_saved_PBO_mean = round(mean(LivesSaved_PBO, na.rm = TRUE), 0),
    Lives_saved_PBO_lower95 = round(qL(LivesSaved_PBO), 0),
    Lives_saved_PBO_upper95 = round(qU(LivesSaved_PBO), 0),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Lives_saved_PBO_95CrI = fmt_ci_table(
      Lives_saved_PBO_mean,
      Lives_saved_PBO_lower95,
      Lives_saved_PBO_upper95,
      digits = 0
    )
  ) %>%
  dplyr::arrange(Year)

print_table_full <- function(title, tab) {
  cat("\n\n", title, "\n", sep = "")
  print(tab, n = Inf, width = Inf)
}

# Optional: save all printed tables as CSV files in your working directory.
dir.create(tables_output_dir, showWarnings = FALSE, recursive = TRUE)

readr::write_csv(table_ITN_prov,      file.path(tables_output_dir, "table_ITN_provinces_plus_DRC.csv"))
readr::write_csv(table_LS_prov,       file.path(tables_output_dir, "table_lives_saved_ITN_provinces_plus_DRC.csv"))
readr::write_csv(table_RES_prov,      file.path(tables_output_dir, "table_resistance_cases_provinces_plus_DRC.csv"))
readr::write_csv(table_DEAD_prov,     file.path(tables_output_dir, "table_resistance_deaths_provinces_plus_DRC.csv"))
readr::write_csv(table_PBO_prov,      file.path(tables_output_dir, "table_PBO_cases_provinces_plus_DRC.csv"))
readr::write_csv(table_PBO_LS_prov,   file.path(tables_output_dir, "table_PBO_lives_saved_provinces_plus_DRC.csv"))
readr::write_csv(table_ITN_five,      file.path(tables_output_dir, "table_ITN_four_provinces_plus_DRC.csv"))
readr::write_csv(table_LS_five,       file.path(tables_output_dir, "table_lives_saved_ITN_four_provinces_plus_DRC.csv"))
readr::write_csv(table_RES_five,      file.path(tables_output_dir, "table_resistance_cases_four_provinces_plus_DRC.csv"))
readr::write_csv(table_DEAD_five,     file.path(tables_output_dir, "table_resistance_deaths_four_provinces_plus_DRC.csv"))
readr::write_csv(table_PBO_five,      file.path(tables_output_dir, "table_PBO_cases_four_provinces_plus_DRC.csv"))
readr::write_csv(table_PBO_LS_five,   file.path(tables_output_dir, "table_PBO_lives_saved_four_provinces_plus_DRC.csv"))
readr::write_csv(table_ITN_nat,       file.path(tables_output_dir, "table_ITN_national_yearly.csv"))
readr::write_csv(table_LS_nat,        file.path(tables_output_dir, "table_lives_saved_ITN_national_yearly.csv"))
readr::write_csv(table_RES_nat,       file.path(tables_output_dir, "table_resistance_cases_national_yearly.csv"))
readr::write_csv(table_DEAD_nat,      file.path(tables_output_dir, "table_resistance_deaths_national_yearly.csv"))
readr::write_csv(table_PBO_nat,       file.path(tables_output_dir, "table_PBO_cases_national_yearly.csv"))
readr::write_csv(table_PBO_LS_nat,    file.path(tables_output_dir, "table_PBO_lives_saved_national_yearly.csv"))

# --- PRINT EVERYTHING ---

# Provinces (26)
print(p_ITN_prov)
print(p_LS_prov)
print(p_RES_prov)
print(p_DEAD_prov)
print(p_PBO_prov)
print(p_PBO_LS_prov)

# 4 provinces + DRC (mean)
print(p_ITN_five)
print(p_LS_five)
print(p_RES_five)
print(p_DEAD_five)
print(p_PBO_five)
print(p_PBO_LS_five)

# National over time
print(p_ITN_nat)
print(p_LS_nat)
print(p_RES_nat)
print(p_DEAD_nat)
print(p_PBO_nat)
print(p_PBO_LS_nat)

# --- PRINT TABLES BEHIND EACH PLOT ---

# Province tables: 26 provinces + DRC mean
print_table_full("=== TABLE 1A. Provinces + DRC mean: cases averted due to ITNs ===", table_ITN_prov)
print_table_full("=== TABLE 1B. Provinces + DRC mean: lives saved due to ITNs ===", table_LS_prov)
print_table_full("=== TABLE 1C. Provinces + DRC mean: case increase due to resistance ===", table_RES_prov)
print_table_full("=== TABLE 1D. Provinces + DRC mean: deaths caused by resistance ===", table_DEAD_prov)
print_table_full("=== TABLE 1E. Provinces + DRC mean: cases averted due to PBO ===", table_PBO_prov)
print_table_full("=== TABLE 1F. Provinces + DRC mean: lives saved due to PBO ===", table_PBO_LS_prov)

# Four selected provinces + DRC mean tables
print_table_full("=== TABLE 2A. Four provinces + DRC mean: cases averted due to ITNs ===", table_ITN_five)
print_table_full("=== TABLE 2B. Four provinces + DRC mean: lives saved due to ITNs ===", table_LS_five)
print_table_full("=== TABLE 2C. Four provinces + DRC mean: case increase due to resistance ===", table_RES_five)
print_table_full("=== TABLE 2D. Four provinces + DRC mean: deaths caused by resistance ===", table_DEAD_five)
print_table_full("=== TABLE 2E. Four provinces + DRC mean: cases averted due to PBO ===", table_PBO_five)
print_table_full("=== TABLE 2F. Four provinces + DRC mean: lives saved due to PBO ===", table_PBO_LS_five)

# National yearly tables
print_table_full("=== TABLE 3A. DRC national yearly: cases averted due to ITNs ===", table_ITN_nat)
print_table_full("=== TABLE 3B. DRC national yearly: lives saved due to ITNs ===", table_LS_nat)
print_table_full("=== TABLE 3C. DRC national yearly: case increase due to resistance ===", table_RES_nat)
print_table_full("=== TABLE 3D. DRC national yearly: deaths caused by resistance ===", table_DEAD_nat)
print_table_full("=== TABLE 3E. DRC national yearly: cases averted due to PBO ===", table_PBO_nat)
print_table_full("=== TABLE 3F. DRC national yearly: lives saved due to PBO ===", table_PBO_LS_nat)

message("All summary tables were also saved as CSV files in: ", tables_output_dir)

# ============================================================
# COMPLETE PBO IMPACT SECTION
# Observed pyrethroid-PBO mix vs pyrethroid-only counterfactual
# ============================================================
# This section prints additional cases and deaths averted by PBO nets.
#
# Correct comparison:
#   Pyrethroid-only counterfactual - observed pyrethroid/PBO mix
#
# Therefore:
#   Additional cases averted by PBO = Sum_WithNets_PyOnly_s - Sum_WithNets_Resist_s
#   Additional deaths averted by PBO = Sum_Deaths_WithNets_PyOnly_s - Sum_Deaths_WithNets_Resist_s
#
# Important wording:
#   Use "additional cases/deaths averted by PBO", not "case increase due to PBO".
#   PBO is protective, so the additional burden is avoided, not increased.
# ============================================================

# ----------------------------
# Output folder and print helper
# ----------------------------
if (!exists("tables_output_dir")) {
  }
dir.create(tables_output_dir, showWarnings = FALSE, recursive = TRUE)

if (!exists("print_table_full")) {
  print_table_full <- function(title, tab) {
    cat("\n\n", title, "\n", sep = "")
    print(tab, n = Inf, width = Inf)
  }
}

if (!exists("eps")) eps <- 1e-9

# ----------------------------
# Formatting helpers
# ----------------------------
qL <- function(x) quantile(x, 0.025, na.rm = TRUE, names = FALSE)
qU <- function(x) quantile(x, 0.975, na.rm = TRUE, names = FALSE)

fmt_big <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    format(round(x), big.mark = ",", scientific = FALSE)
  )
}

fmt_pct <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.1f", x))
}

fmt_ci_count <- function(mean, lower, upper) {
  paste0(
    fmt_big(mean),
    " (95% CrI: ",
    fmt_big(lower),
    "–",
    fmt_big(upper),
    ")"
  )
}

fmt_ci_pct <- function(mean, lower, upper) {
  paste0(
    fmt_pct(mean),
    "% (95% CrI: ",
    fmt_pct(lower),
    "–",
    fmt_pct(upper),
    "%)"
  )
}

# ----------------------------
# Safety checks
# ----------------------------
required_cols_nat_pbo <- c(
  "Year",
  "Draw",
  "Sum_WithNets_PyOnly_s",
  "Sum_WithNets_Resist_s",
  "Sum_Deaths_WithNets_PyOnly_s",
  "Sum_Deaths_WithNets_Resist_s"
)

missing_cols_nat_pbo <- setdiff(required_cols_nat_pbo, names(national_yearly_pw))

if (length(missing_cols_nat_pbo) > 0) {
  stop(
    "These columns are missing from national_yearly_pw: ",
    paste(missing_cols_nat_pbo, collapse = ", "),
    "\nRun the full updated simulation code with the pyrethroid-only counterfactual first."
  )
}

required_cols_prov_pbo <- c(
  "Region",
  "ProvKey",
  "Draw",
  "Year",
  "Cases_WithNets_PyOnly_s",
  "Cases_WithNets_Resist_s",
  "Deaths_WithNets_PyOnly_s",
  "Deaths_WithNets_Resist_s"
)

missing_cols_prov_pbo <- setdiff(required_cols_prov_pbo, names(core_yearly_pw))

if (length(missing_cols_prov_pbo) > 0) {
  stop(
    "These columns are missing from core_yearly_pw: ",
    paste(missing_cols_prov_pbo, collapse = ", "),
    "\nRun the full updated simulation code with the pyrethroid-only counterfactual first."
  )
}

# ----------------------------
# PBO years and PBO percentage by year
# ----------------------------
if (exists("pbo_years_active") && length(pbo_years_active) > 0) {
  pbo_years_for_summary <- sort(unique(pbo_years_active))
} else {
  pbo_years_for_summary <- 2019:2024
}

pbo_years_for_summary <- pbo_years_for_summary[
  pbo_years_for_summary %in% unique(national_yearly_pw$Year)
]

if (length(pbo_years_for_summary) == 0) {
  stop("No valid PBO years were found in national_yearly_pw.")
}

pbo_percent_yearly <- pbo_mix %>%
  dplyr::filter(Year %in% pbo_years_for_summary) %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(
    PBO_percent = 100 * mean(y, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 1) NATIONAL YEARLY ADDITIONAL CASES AND DEATHS AVERTED BY PBO
# ============================================================

pbo_additional_nat_yearly_draws <- national_yearly_pw %>%
  dplyr::filter(Year %in% pbo_years_for_summary) %>%
  dplyr::mutate(
    Additional_cases_averted_by_PBO =
      pmax(Sum_WithNets_PyOnly_s - Sum_WithNets_Resist_s, 0),

    Additional_deaths_averted_by_PBO =
      pmax(Sum_Deaths_WithNets_PyOnly_s - Sum_Deaths_WithNets_Resist_s, 0),

    Percent_cases_averted_by_PBO =
      100 * Additional_cases_averted_by_PBO /
      pmax(Sum_WithNets_PyOnly_s, eps),

    Percent_deaths_averted_by_PBO =
      100 * Additional_deaths_averted_by_PBO /
      pmax(Sum_Deaths_WithNets_PyOnly_s, eps)
  )

table_PBO_additional_nat_yearly <- pbo_additional_nat_yearly_draws %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(
    Additional_cases_averted_by_PBO_mean =
      round(mean(Additional_cases_averted_by_PBO, na.rm = TRUE), 0),
    Additional_cases_averted_by_PBO_lower95 =
      round(qL(Additional_cases_averted_by_PBO), 0),
    Additional_cases_averted_by_PBO_upper95 =
      round(qU(Additional_cases_averted_by_PBO), 0),

    Percent_cases_averted_by_PBO_mean =
      round(mean(Percent_cases_averted_by_PBO, na.rm = TRUE), 1),
    Percent_cases_averted_by_PBO_lower95 =
      round(qL(Percent_cases_averted_by_PBO), 1),
    Percent_cases_averted_by_PBO_upper95 =
      round(qU(Percent_cases_averted_by_PBO), 1),

    Additional_deaths_averted_by_PBO_mean =
      round(mean(Additional_deaths_averted_by_PBO, na.rm = TRUE), 0),
    Additional_deaths_averted_by_PBO_lower95 =
      round(qL(Additional_deaths_averted_by_PBO), 0),
    Additional_deaths_averted_by_PBO_upper95 =
      round(qU(Additional_deaths_averted_by_PBO), 0),

    Percent_deaths_averted_by_PBO_mean =
      round(mean(Percent_deaths_averted_by_PBO, na.rm = TRUE), 1),
    Percent_deaths_averted_by_PBO_lower95 =
      round(qL(Percent_deaths_averted_by_PBO), 1),
    Percent_deaths_averted_by_PBO_upper95 =
      round(qU(Percent_deaths_averted_by_PBO), 1),

    .groups = "drop"
  ) %>%
  dplyr::left_join(pbo_percent_yearly, by = "Year") %>%
  dplyr::mutate(
    PBO_percent = round(PBO_percent, 1),

    Additional_cases_averted_by_PBO_95CrI =
      fmt_ci_count(
        Additional_cases_averted_by_PBO_mean,
        Additional_cases_averted_by_PBO_lower95,
        Additional_cases_averted_by_PBO_upper95
      ),

    Percent_cases_averted_by_PBO_95CrI =
      fmt_ci_pct(
        Percent_cases_averted_by_PBO_mean,
        Percent_cases_averted_by_PBO_lower95,
        Percent_cases_averted_by_PBO_upper95
      ),

    Additional_deaths_averted_by_PBO_95CrI =
      fmt_ci_count(
        Additional_deaths_averted_by_PBO_mean,
        Additional_deaths_averted_by_PBO_lower95,
        Additional_deaths_averted_by_PBO_upper95
      ),

    Percent_deaths_averted_by_PBO_95CrI =
      fmt_ci_pct(
        Percent_deaths_averted_by_PBO_mean,
        Percent_deaths_averted_by_PBO_lower95,
        Percent_deaths_averted_by_PBO_upper95
      )
  ) %>%
  dplyr::arrange(Year)

print_table_full(
  "=== NATIONAL YEARLY ADDITIONAL CASES AND DEATHS AVERTED BY PBO ===",
  table_PBO_additional_nat_yearly
)

# ============================================================
# 2) NATIONAL TOTAL ADDITIONAL CASES AND DEATHS AVERTED BY PBO
#    OVER THE PBO PERIOD
# ============================================================

table_PBO_additional_nat_total <- pbo_additional_nat_yearly_draws %>%
  dplyr::group_by(Draw) %>%
  dplyr::summarise(
    Total_additional_cases_averted_by_PBO =
      sum(Additional_cases_averted_by_PBO, na.rm = TRUE),
    Total_additional_deaths_averted_by_PBO =
      sum(Additional_deaths_averted_by_PBO, na.rm = TRUE),

    Total_pyrethroid_only_cases =
      sum(Sum_WithNets_PyOnly_s, na.rm = TRUE),
    Total_pyrethroid_only_deaths =
      sum(Sum_Deaths_WithNets_PyOnly_s, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Percent_cases_averted_by_PBO =
      100 * Total_additional_cases_averted_by_PBO /
      pmax(Total_pyrethroid_only_cases, eps),

    Percent_deaths_averted_by_PBO =
      100 * Total_additional_deaths_averted_by_PBO /
      pmax(Total_pyrethroid_only_deaths, eps)
  ) %>%
  dplyr::summarise(
    Total_additional_cases_averted_by_PBO_mean =
      round(mean(Total_additional_cases_averted_by_PBO, na.rm = TRUE), 0),
    Total_additional_cases_averted_by_PBO_lower95 =
      round(qL(Total_additional_cases_averted_by_PBO), 0),
    Total_additional_cases_averted_by_PBO_upper95 =
      round(qU(Total_additional_cases_averted_by_PBO), 0),

    Percent_cases_averted_by_PBO_mean =
      round(mean(Percent_cases_averted_by_PBO, na.rm = TRUE), 1),
    Percent_cases_averted_by_PBO_lower95 =
      round(qL(Percent_cases_averted_by_PBO), 1),
    Percent_cases_averted_by_PBO_upper95 =
      round(qU(Percent_cases_averted_by_PBO), 1),

    Total_additional_deaths_averted_by_PBO_mean =
      round(mean(Total_additional_deaths_averted_by_PBO, na.rm = TRUE), 0),
    Total_additional_deaths_averted_by_PBO_lower95 =
      round(qL(Total_additional_deaths_averted_by_PBO), 0),
    Total_additional_deaths_averted_by_PBO_upper95 =
      round(qU(Total_additional_deaths_averted_by_PBO), 0),

    Percent_deaths_averted_by_PBO_mean =
      round(mean(Percent_deaths_averted_by_PBO, na.rm = TRUE), 1),
    Percent_deaths_averted_by_PBO_lower95 =
      round(qL(Percent_deaths_averted_by_PBO), 1),
    Percent_deaths_averted_by_PBO_upper95 =
      round(qU(Percent_deaths_averted_by_PBO), 1),

    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Years_included = paste(range(pbo_years_for_summary), collapse = "–"),

    Total_additional_cases_averted_by_PBO_95CrI =
      fmt_ci_count(
        Total_additional_cases_averted_by_PBO_mean,
        Total_additional_cases_averted_by_PBO_lower95,
        Total_additional_cases_averted_by_PBO_upper95
      ),

    Percent_cases_averted_by_PBO_95CrI =
      fmt_ci_pct(
        Percent_cases_averted_by_PBO_mean,
        Percent_cases_averted_by_PBO_lower95,
        Percent_cases_averted_by_PBO_upper95
      ),

    Total_additional_deaths_averted_by_PBO_95CrI =
      fmt_ci_count(
        Total_additional_deaths_averted_by_PBO_mean,
        Total_additional_deaths_averted_by_PBO_lower95,
        Total_additional_deaths_averted_by_PBO_upper95
      ),

    Percent_deaths_averted_by_PBO_95CrI =
      fmt_ci_pct(
        Percent_deaths_averted_by_PBO_mean,
        Percent_deaths_averted_by_PBO_lower95,
        Percent_deaths_averted_by_PBO_upper95
      )
  ) %>%
  dplyr::relocate(Years_included)

print_table_full(
  "=== NATIONAL TOTAL ADDITIONAL CASES AND DEATHS AVERTED BY PBO OVER THE PBO PERIOD ===",
  table_PBO_additional_nat_total
)

# ============================================================
# 3) NATIONAL 2024 ADDITIONAL CASES AND DEATHS AVERTED BY PBO
# ============================================================

table_PBO_additional_nat_2024 <- table_PBO_additional_nat_yearly %>%
  dplyr::filter(Year == 2024)

print_table_full(
  "=== NATIONAL 2024 ADDITIONAL CASES AND DEATHS AVERTED BY PBO ===",
  table_PBO_additional_nat_2024
)

# ============================================================
# 4) PROVINCE-YEAR ADDITIONAL CASES AND DEATHS AVERTED BY PBO
# ============================================================

pbo_additional_prov_yearly_draws <- core_yearly_pw %>%
  dplyr::filter(Year %in% pbo_years_for_summary) %>%
  dplyr::mutate(
    Additional_cases_averted_by_PBO =
      pmax(Cases_WithNets_PyOnly_s - Cases_WithNets_Resist_s, 0),

    Additional_deaths_averted_by_PBO =
      pmax(Deaths_WithNets_PyOnly_s - Deaths_WithNets_Resist_s, 0),

    Percent_cases_averted_by_PBO =
      100 * Additional_cases_averted_by_PBO /
      pmax(Cases_WithNets_PyOnly_s, eps),

    Percent_deaths_averted_by_PBO =
      100 * Additional_deaths_averted_by_PBO /
      pmax(Deaths_WithNets_PyOnly_s, eps)
  )

table_PBO_additional_prov_yearly <- pbo_additional_prov_yearly_draws %>%
  dplyr::group_by(Region, ProvKey, Year) %>%
  dplyr::summarise(
    Additional_cases_averted_by_PBO_mean =
      round(mean(Additional_cases_averted_by_PBO, na.rm = TRUE), 0),
    Additional_cases_averted_by_PBO_lower95 =
      round(qL(Additional_cases_averted_by_PBO), 0),
    Additional_cases_averted_by_PBO_upper95 =
      round(qU(Additional_cases_averted_by_PBO), 0),

    Percent_cases_averted_by_PBO_mean =
      round(mean(Percent_cases_averted_by_PBO, na.rm = TRUE), 1),
    Percent_cases_averted_by_PBO_lower95 =
      round(qL(Percent_cases_averted_by_PBO), 1),
    Percent_cases_averted_by_PBO_upper95 =
      round(qU(Percent_cases_averted_by_PBO), 1),

    Additional_deaths_averted_by_PBO_mean =
      round(mean(Additional_deaths_averted_by_PBO, na.rm = TRUE), 0),
    Additional_deaths_averted_by_PBO_lower95 =
      round(qL(Additional_deaths_averted_by_PBO), 0),
    Additional_deaths_averted_by_PBO_upper95 =
      round(qU(Additional_deaths_averted_by_PBO), 0),

    Percent_deaths_averted_by_PBO_mean =
      round(mean(Percent_deaths_averted_by_PBO, na.rm = TRUE), 1),
    Percent_deaths_averted_by_PBO_lower95 =
      round(qL(Percent_deaths_averted_by_PBO), 1),
    Percent_deaths_averted_by_PBO_upper95 =
      round(qU(Percent_deaths_averted_by_PBO), 1),

    .groups = "drop"
  ) %>%
  dplyr::left_join(pbo_percent_yearly, by = "Year") %>%
  dplyr::mutate(
    PBO_percent = round(PBO_percent, 1),

    Additional_cases_averted_by_PBO_95CrI =
      fmt_ci_count(
        Additional_cases_averted_by_PBO_mean,
        Additional_cases_averted_by_PBO_lower95,
        Additional_cases_averted_by_PBO_upper95
      ),

    Percent_cases_averted_by_PBO_95CrI =
      fmt_ci_pct(
        Percent_cases_averted_by_PBO_mean,
        Percent_cases_averted_by_PBO_lower95,
        Percent_cases_averted_by_PBO_upper95
      ),

    Additional_deaths_averted_by_PBO_95CrI =
      fmt_ci_count(
        Additional_deaths_averted_by_PBO_mean,
        Additional_deaths_averted_by_PBO_lower95,
        Additional_deaths_averted_by_PBO_upper95
      ),

    Percent_deaths_averted_by_PBO_95CrI =
      fmt_ci_pct(
        Percent_deaths_averted_by_PBO_mean,
        Percent_deaths_averted_by_PBO_lower95,
        Percent_deaths_averted_by_PBO_upper95
      )
  ) %>%
  dplyr::arrange(Region, Year)

print_table_full(
  "=== PROVINCE-YEAR ADDITIONAL CASES AND DEATHS AVERTED BY PBO ===",
  table_PBO_additional_prov_yearly
)

# ============================================================
# 5) PROVINCE-LEVEL 2024 ADDITIONAL CASES AND DEATHS AVERTED BY PBO
# ============================================================

table_PBO_additional_prov_2024 <- table_PBO_additional_prov_yearly %>%
  dplyr::filter(Year == 2024) %>%
  dplyr::arrange(Region)

print_table_full(
  "=== PROVINCE-LEVEL 2024 ADDITIONAL CASES AND DEATHS AVERTED BY PBO ===",
  table_PBO_additional_prov_2024
)

# ============================================================
# 6) FOUR SELECTED PROVINCES + DRC MEAN FOR 2024
# ============================================================

if (!exists("four_labels")) {
  four_labels <- c("Haut-Katanga", "Kasaï-Central", "Kinshasa", "Sud-Kivu")
}

drc_2024_row <- table_PBO_additional_nat_2024 %>%
  dplyr::mutate(
    Region = "DRC (mean)",
    ProvKey = "DRC_mean"
  ) %>%
  dplyr::select(
    Region, ProvKey, Year, PBO_percent,
    Additional_cases_averted_by_PBO_mean,
    Additional_cases_averted_by_PBO_lower95,
    Additional_cases_averted_by_PBO_upper95,
    Percent_cases_averted_by_PBO_mean,
    Percent_cases_averted_by_PBO_lower95,
    Percent_cases_averted_by_PBO_upper95,
    Additional_deaths_averted_by_PBO_mean,
    Additional_deaths_averted_by_PBO_lower95,
    Additional_deaths_averted_by_PBO_upper95,
    Percent_deaths_averted_by_PBO_mean,
    Percent_deaths_averted_by_PBO_lower95,
    Percent_deaths_averted_by_PBO_upper95,
    Additional_cases_averted_by_PBO_95CrI,
    Percent_cases_averted_by_PBO_95CrI,
    Additional_deaths_averted_by_PBO_95CrI,
    Percent_deaths_averted_by_PBO_95CrI
  )

table_PBO_additional_five_2024 <- table_PBO_additional_prov_2024 %>%
  dplyr::filter(Region %in% four_labels) %>%
  dplyr::bind_rows(drc_2024_row) %>%
  dplyr::mutate(
    Region = factor(Region, levels = c(four_labels, "DRC (mean)"))
  ) %>%
  dplyr::arrange(Region) %>%
  dplyr::mutate(Region = as.character(Region))

print_table_full(
  "=== FOUR PROVINCES + DRC MEAN, 2024: ADDITIONAL CASES AND DEATHS AVERTED BY PBO ===",
  table_PBO_additional_five_2024
)

# ============================================================
# 7) MANUSCRIPT SENTENCE FOR 2024
# ============================================================

pbo_percent_2024 <- pbo_percent_yearly %>%
  dplyr::filter(Year == 2024) %>%
  dplyr::pull(PBO_percent)

if (length(pbo_percent_2024) == 0 || is.na(pbo_percent_2024)) {
  warning("No valid PBO percentage found for 2024 in pbo_mix.")
  pbo_percent_2024 <- NA_real_
}

sentence_PBO_2024 <- paste0(
  "In 2024, ",
  fmt_pct(pbo_percent_2024),
  "% of the nets in circulation were predicted to be pyrethroid-PBO, ",
  "and these nets were predicted to have averted an additional ",
  table_PBO_additional_nat_2024$Additional_cases_averted_by_PBO_95CrI,
  " malaria cases and ",
  table_PBO_additional_nat_2024$Additional_deaths_averted_by_PBO_95CrI,
  " deaths compared with a counterfactual scenario in which pyrethroid-only nets had been used."
)

cat("\n\n", sentence_PBO_2024, "\n\n")

# ============================================================
# 8) SAVE PBO ADDITIONAL IMPACT TABLES
# ============================================================

readr::write_csv(
  table_PBO_additional_nat_yearly,
  file.path(tables_output_dir, "table_PBO_additional_cases_deaths_national_yearly.csv")
)

readr::write_csv(
  table_PBO_additional_nat_total,
  file.path(tables_output_dir, "table_PBO_additional_cases_deaths_national_total.csv")
)

readr::write_csv(
  table_PBO_additional_nat_2024,
  file.path(tables_output_dir, "table_PBO_additional_cases_deaths_national_2024.csv")
)

readr::write_csv(
  table_PBO_additional_prov_yearly,
  file.path(tables_output_dir, "table_PBO_additional_cases_deaths_province_yearly.csv")
)

readr::write_csv(
  table_PBO_additional_prov_2024,
  file.path(tables_output_dir, "table_PBO_additional_cases_deaths_province_2024.csv")
)

readr::write_csv(
  table_PBO_additional_five_2024,
  file.path(tables_output_dir, "table_PBO_additional_cases_deaths_four_provinces_plus_DRC_2024.csv")
)

writeLines(
  sentence_PBO_2024,
  con = file.path(tables_output_dir, "sentence_PBO_impact_2024.txt")
)

message("PBO additional cases/deaths tables and sentence were saved in: ", tables_output_dir)


########### -----------------------------------------------------------------------------
# ============================================================
# CORRECTED 2024 PBO IMPACT USING WHO-CALIBRATED OUTPUTS
# Compared with pyrethroid-only nets
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(scales)
})

if (!exists("eps")) eps <- 1e-9

qL <- function(x) quantile(x, 0.025, na.rm = TRUE, names = FALSE)
qU <- function(x) quantile(x, 0.975, na.rm = TRUE, names = FALSE)

fmt_big <- function(x) {
  format(round(x), big.mark = ",", scientific = FALSE)
}

fmt_pct <- function(x) {
  sprintf("%.1f", x)
}

# ------------------------------------------------------------
# Safety checks
# ------------------------------------------------------------

required_cols_core_pbo_who <- c(
  "Year",
  "Draw",
  "Cases_WithNets_PyOnly_who",
  "Cases_WithNets_Resist_who",
  "Deaths_WithNets_PyOnly_who",
  "Deaths_WithNets_Resist_who"
)

missing_cols_core_pbo_who <- setdiff(required_cols_core_pbo_who, names(core_yearly_pw))

if (length(missing_cols_core_pbo_who) > 0) {
  stop(
    "These WHO-calibrated columns are missing from core_yearly_pw: ",
    paste(missing_cols_core_pbo_who, collapse = ", "),
    "\nRun the full code until the section where core_yearly_pw is created with *_who variables."
  )
}

if (!exists("pbo_mix")) {
  stop("The object pbo_mix was not found. Run the section that reads type_itn_pbo.xlsx first.")
}

# ------------------------------------------------------------
# PBO percentage in 2024
# ------------------------------------------------------------

pbo_percent_2024 <- pbo_mix %>%
  dplyr::filter(Year == 2024) %>%
  dplyr::summarise(
    PBO_percent = 100 * mean(y, na.rm = TRUE)
  ) %>%
  dplyr::pull(PBO_percent)

# ------------------------------------------------------------
# National WHO-calibrated PBO impact in 2024
# ------------------------------------------------------------

pbo_impact_2024_who_draws <- core_yearly_pw %>%
  dplyr::filter(Year == 2024) %>%
  dplyr::group_by(Year, Draw) %>%
  dplyr::summarise(
    Cases_WithNets_PyOnly_who =
      sum(Cases_WithNets_PyOnly_who, na.rm = TRUE),
    
    Cases_WithNets_Resist_who =
      sum(Cases_WithNets_Resist_who, na.rm = TRUE),
    
    Deaths_WithNets_PyOnly_who =
      sum(Deaths_WithNets_PyOnly_who, na.rm = TRUE),
    
    Deaths_WithNets_Resist_who =
      sum(Deaths_WithNets_Resist_who, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Additional_cases_averted_by_PBO =
      pmax(Cases_WithNets_PyOnly_who - Cases_WithNets_Resist_who, 0),
    
    Additional_deaths_averted_by_PBO =
      pmax(Deaths_WithNets_PyOnly_who - Deaths_WithNets_Resist_who, 0),
    
    Percent_cases_averted_by_PBO =
      100 * Additional_cases_averted_by_PBO /
      pmax(Cases_WithNets_PyOnly_who, eps),
    
    Percent_deaths_averted_by_PBO =
      100 * Additional_deaths_averted_by_PBO /
      pmax(Deaths_WithNets_PyOnly_who, eps)
  )

# ------------------------------------------------------------
# Summary across posterior draws
# ------------------------------------------------------------

pbo_impact_2024_who_summary <- pbo_impact_2024_who_draws %>%
  dplyr::summarise(
    Year = 2024,
    PBO_percent = pbo_percent_2024,
    
    cases_mean = mean(Additional_cases_averted_by_PBO, na.rm = TRUE),
    cases_lower95 = qL(Additional_cases_averted_by_PBO),
    cases_upper95 = qU(Additional_cases_averted_by_PBO),
    
    deaths_mean = mean(Additional_deaths_averted_by_PBO, na.rm = TRUE),
    deaths_lower95 = qL(Additional_deaths_averted_by_PBO),
    deaths_upper95 = qU(Additional_deaths_averted_by_PBO),
    
    pct_cases_mean = mean(Percent_cases_averted_by_PBO, na.rm = TRUE),
    pct_cases_lower95 = qL(Percent_cases_averted_by_PBO),
    pct_cases_upper95 = qU(Percent_cases_averted_by_PBO),
    
    pct_deaths_mean = mean(Percent_deaths_averted_by_PBO, na.rm = TRUE),
    pct_deaths_lower95 = qL(Percent_deaths_averted_by_PBO),
    pct_deaths_upper95 = qU(Percent_deaths_averted_by_PBO)
  ) %>%
  dplyr::mutate(
    cases_95CrI = paste0(
      fmt_big(cases_mean),
      " (95% CrI: ",
      fmt_big(cases_lower95),
      "–",
      fmt_big(cases_upper95),
      ")"
    ),
    
    deaths_95CrI = paste0(
      fmt_big(deaths_mean),
      " (95% CrI: ",
      fmt_big(deaths_lower95),
      "–",
      fmt_big(deaths_upper95),
      ")"
    ),
    
    pct_cases_95CrI = paste0(
      fmt_pct(pct_cases_mean),
      "% (95% CrI: ",
      fmt_pct(pct_cases_lower95),
      "–",
      fmt_pct(pct_cases_upper95),
      "%)"
    )
  )

# ------------------------------------------------------------
# Print values for manuscript
# ------------------------------------------------------------

cat("\n=== CORRECTED VALUES FOR MANUSCRIPT SENTENCE ===\n")
cat("X   = ", fmt_pct(pbo_impact_2024_who_summary$PBO_percent), "%\n", sep = "")
cat("XXX = ", pbo_impact_2024_who_summary$cases_95CrI, " cases\n", sep = "")
cat("YYY = ", pbo_impact_2024_who_summary$deaths_95CrI, " deaths\n", sep = "")
cat("Percent additional cases averted = ",
    pbo_impact_2024_who_summary$pct_cases_95CrI, "\n", sep = "")

sentence_PBO_2024_corrected <- paste0(
  "In 2024, ",
  fmt_pct(pbo_impact_2024_who_summary$PBO_percent),
  "\\% of the nets in circulation were predicted to be pyrethroid-PBO and ",
  "were predicted to have averted an additional ",
  pbo_impact_2024_who_summary$cases_95CrI,
  " cases and ",
  pbo_impact_2024_who_summary$deaths_95CrI,
  " deaths compared with a counterfactual scenario in which pyrethroid-only nets had remained in use."
)

cat("\n=== CORRECTED MANUSCRIPT SENTENCE ===\n")
cat(sentence_PBO_2024_corrected, "\n")

########## -----------------------------------------------------------------------------

# ============================================================
# CORRECTED SENTENCE:
# Switch from pyrethroid-only to pyrethroid-PBO
# Additional cases and deaths averted, WHO-calibrated
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

if (!exists("eps")) eps <- 1e-9

qL <- function(x) quantile(x, 0.025, na.rm = TRUE, names = FALSE)
qU <- function(x) quantile(x, 0.975, na.rm = TRUE, names = FALSE)

fmt_big <- function(x) {
  format(round(x), big.mark = ",", scientific = FALSE)
}

fmt_pct <- function(x) {
  sprintf("%.1f", x)
}

fmt_ci_count <- function(mean, lower, upper) {
  paste0(
    fmt_big(mean),
    " (95% CrI: ",
    fmt_big(lower),
    "–",
    fmt_big(upper),
    ")"
  )
}

fmt_ci_pct <- function(mean, lower, upper) {
  paste0(
    fmt_pct(mean),
    "% (95% CrI: ",
    fmt_pct(lower),
    "–",
    fmt_pct(upper),
    "%)"
  )
}

# ------------------------------------------------------------
# Safety checks
# ------------------------------------------------------------

required_cols <- c(
  "Region",
  "ProvKey",
  "Year",
  "Draw",
  "Cases_NoNets_who",
  "Cases_WithNets_PyOnly_who",
  "Cases_WithNets_Resist_who",
  "Deaths_WithNets_PyOnly_who",
  "Deaths_WithNets_Resist_who"
)

missing_cols <- setdiff(required_cols, names(core_yearly_pw))

if (length(missing_cols) > 0) {
  stop(
    "These columns are missing from core_yearly_pw: ",
    paste(missing_cols, collapse = ", "),
    "\nRun the full simulation code until the *_who variables are created."
  )
}

# ------------------------------------------------------------
# Define PBO period
# ------------------------------------------------------------

if (exists("pbo_years_active") && length(pbo_years_active) > 0) {
  pbo_years_for_sentence <- sort(unique(pbo_years_active))
} else {
  pbo_years_for_sentence <- pbo_mix %>%
    dplyr::filter(!is.na(y), y > 0) %>%
    dplyr::pull(Year) %>%
    unique() %>%
    sort()
}

pbo_years_for_sentence <- pbo_years_for_sentence[
  pbo_years_for_sentence %in% unique(core_yearly_pw$Year)
]

if (length(pbo_years_for_sentence) == 0) {
  stop("No PBO years were found. Check pbo_mix or pbo_years_active.")
}

cat("\nPBO years used for this sentence: ",
    paste(pbo_years_for_sentence, collapse = ", "), "\n", sep = "")

# ------------------------------------------------------------
# National WHO-calibrated PBO impact over the PBO period
# ------------------------------------------------------------

switch_pbo_nat_draws <- core_yearly_pw %>%
  dplyr::filter(Year %in% pbo_years_for_sentence) %>%
  dplyr::group_by(Draw) %>%
  dplyr::summarise(
    Cases_NoNets_who =
      sum(Cases_NoNets_who, na.rm = TRUE),
    
    Cases_WithNets_PyOnly_who =
      sum(Cases_WithNets_PyOnly_who, na.rm = TRUE),
    
    Cases_WithNets_Resist_who =
      sum(Cases_WithNets_Resist_who, na.rm = TRUE),
    
    Deaths_WithNets_PyOnly_who =
      sum(Deaths_WithNets_PyOnly_who, na.rm = TRUE),
    
    Deaths_WithNets_Resist_who =
      sum(Deaths_WithNets_Resist_who, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # ITN cases averted under pyrethroid-only nets
    Cases_Averted_ITN_PyOnly =
      pmax(Cases_NoNets_who - Cases_WithNets_PyOnly_who, 0),
    
    # ITN cases averted under observed pyrethroid-PBO mix
    Cases_Averted_ITN_PBO =
      pmax(Cases_NoNets_who - Cases_WithNets_Resist_who, 0),
    
    # Percentage of cases averted by ITNs
    Percent_Cases_Averted_ITN_PyOnly =
      100 * Cases_Averted_ITN_PyOnly / pmax(Cases_NoNets_who, eps),
    
    Percent_Cases_Averted_ITN_PBO =
      100 * Cases_Averted_ITN_PBO / pmax(Cases_NoNets_who, eps),
    
    # Increase in percentage points
    National_Increase_Percentage_Points =
      Percent_Cases_Averted_ITN_PBO - Percent_Cases_Averted_ITN_PyOnly,
    
    # Additional cases and deaths averted by switching to PBO
    Additional_cases_averted_by_PBO =
      pmax(Cases_WithNets_PyOnly_who - Cases_WithNets_Resist_who, 0),
    
    Additional_deaths_averted_by_PBO =
      pmax(Deaths_WithNets_PyOnly_who - Deaths_WithNets_Resist_who, 0)
  )

switch_pbo_nat_summary <- switch_pbo_nat_draws %>%
  dplyr::summarise(
    PyOnly_percent_mean =
      mean(Percent_Cases_Averted_ITN_PyOnly, na.rm = TRUE),
    PyOnly_percent_lower95 =
      qL(Percent_Cases_Averted_ITN_PyOnly),
    PyOnly_percent_upper95 =
      qU(Percent_Cases_Averted_ITN_PyOnly),
    
    PBO_percent_mean =
      mean(Percent_Cases_Averted_ITN_PBO, na.rm = TRUE),
    PBO_percent_lower95 =
      qL(Percent_Cases_Averted_ITN_PBO),
    PBO_percent_upper95 =
      qU(Percent_Cases_Averted_ITN_PBO),
    
    National_increase_mean =
      mean(National_Increase_Percentage_Points, na.rm = TRUE),
    National_increase_lower95 =
      qL(National_Increase_Percentage_Points),
    National_increase_upper95 =
      qU(National_Increase_Percentage_Points),
    
    Additional_cases_mean =
      mean(Additional_cases_averted_by_PBO, na.rm = TRUE),
    Additional_cases_lower95 =
      qL(Additional_cases_averted_by_PBO),
    Additional_cases_upper95 =
      qU(Additional_cases_averted_by_PBO),
    
    Additional_deaths_mean =
      mean(Additional_deaths_averted_by_PBO, na.rm = TRUE),
    Additional_deaths_lower95 =
      qL(Additional_deaths_averted_by_PBO),
    Additional_deaths_upper95 =
      qU(Additional_deaths_averted_by_PBO)
  ) %>%
  dplyr::mutate(
    PyOnly_percent_95CrI =
      fmt_ci_pct(PyOnly_percent_mean, PyOnly_percent_lower95, PyOnly_percent_upper95),
    
    PBO_percent_95CrI =
      fmt_ci_pct(PBO_percent_mean, PBO_percent_lower95, PBO_percent_upper95),
    
    National_increase_95CrI =
      fmt_ci_pct(National_increase_mean, National_increase_lower95, National_increase_upper95),
    
    Additional_cases_95CrI =
      fmt_ci_count(Additional_cases_mean, Additional_cases_lower95, Additional_cases_upper95),
    
    Additional_deaths_95CrI =
      fmt_ci_count(Additional_deaths_mean, Additional_deaths_lower95, Additional_deaths_upper95)
  )

# ------------------------------------------------------------
# Print values to replace XXX and YYY
# ------------------------------------------------------------

cat("\n=== CORRECTED VALUES FOR THIS SENTENCE ===\n")
cat("Pyrethroid-only ITN cases averted (%) = ",
    switch_pbo_nat_summary$PyOnly_percent_95CrI, "\n", sep = "")

cat("Pyrethroid-PBO ITN cases averted (%) = ",
    switch_pbo_nat_summary$PBO_percent_95CrI, "\n", sep = "")

cat("National increase = ",
    switch_pbo_nat_summary$National_increase_95CrI,
    " percentage points\n", sep = "")

cat("XXX = ",
    switch_pbo_nat_summary$Additional_cases_95CrI,
    " cases\n", sep = "")

cat("YYY = ",
    switch_pbo_nat_summary$Additional_deaths_95CrI,
    " deaths\n", sep = "")

# ------------------------------------------------------------
# Manuscript sentence
# ------------------------------------------------------------

sentence_switch_pbo <- paste0(
  "The switch from pyrethroid-only to pyrethroid-PBO nets was predicted to have ",
  "substantially increased malaria prevention. During the PBO distribution period, ",
  "the percentage of cases averted by ITNs increased from ",
  fmt_pct(switch_pbo_nat_summary$PyOnly_percent_mean),
  "\\% under a pyrethroid-only counterfactual to ",
  fmt_pct(switch_pbo_nat_summary$PBO_percent_mean),
  "\\% under the observed pyrethroid-PBO mix. Nationally, this represented an increase of ",
  fmt_pct(switch_pbo_nat_summary$National_increase_mean),
  " percentage points, averting an additional ",
  switch_pbo_nat_summary$Additional_cases_95CrI,
  " cases and ",
  switch_pbo_nat_summary$Additional_deaths_95CrI,
  " deaths."
)

cat("\n=== CORRECTED MANUSCRIPT SENTENCE ===\n")
cat(sentence_switch_pbo, "\n")

###### -----------------------------------------------------------------------------

# ============================================================
# PRINT X FOR:
# In 2024 X% of the nets in circulation were predicted to be pyrethroid-PBO
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

fmt_pct <- function(x) sprintf("%.1f", x)

# If pbo_percent_yearly already exists, use it.
# If not, create it from pbo_mix.
if (!exists("pbo_percent_yearly")) {
  pbo_percent_yearly <- pbo_mix %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(
      PBO_percent = 100 * mean(y, na.rm = TRUE),
      .groups = "drop"
    )
}

# Extract 2024 value
pbo_percent_2024 <- pbo_percent_yearly %>%
  dplyr::filter(Year == 2024) %>%
  dplyr::pull(PBO_percent)

if (length(pbo_percent_2024) == 0 || is.na(pbo_percent_2024)) {
  stop("No valid PBO percentage found for 2024. Check pbo_mix or pbo_percent_yearly.")
}

# Print X
cat("\n=== VALUE TO REPLACE X ===\n")
cat("X = ", fmt_pct(pbo_percent_2024), "%\n", sep = "")

# Print manuscript sentence
sentence_pbo_percent_2024 <- paste0(
  "In 2024, ",
  fmt_pct(pbo_percent_2024),
  "\\% of the nets in circulation were predicted to be pyrethroid-PBO."
)

cat("\n=== MANUSCRIPT SENTENCE ===\n")
cat(sentence_pbo_percent_2024, "\n")

############### -----------------------------------------------------------------------------

# ============================================================
# RESISTANCE-ATTRIBUTABLE CASES AND DEATHS
# For manuscript sentence:
# Predicted malaria cases and deaths attributable to pyrethroid resistance
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringi)
})

# ----------------------------
# Helper functions
# ----------------------------

if (!exists("eps")) eps <- 1e-9

qL <- function(x) quantile(x, 0.025, na.rm = TRUE, names = FALSE)
qU <- function(x) quantile(x, 0.975, na.rm = TRUE, names = FALSE)

fmt_big <- function(x) {
  format(round(x), big.mark = ",", scientific = FALSE)
}

fmt_k <- function(x) {
  sprintf("%.1f", x / 1000)
}

fmt_pct <- function(x) {
  sprintf("%.1f", x)
}

fmt_num1 <- function(x) {
  sprintf("%.1f", x)
}

fmt_ci_count <- function(mean, lower, upper) {
  paste0(
    fmt_big(mean),
    " (95% CrI: ",
    fmt_big(lower),
    "–",
    fmt_big(upper),
    ")"
  )
}

fmt_ci_k <- function(mean, lower, upper) {
  paste0(
    fmt_k(mean),
    "k (95% CrI: ",
    fmt_k(lower),
    "–",
    fmt_k(upper),
    "k)"
  )
}

fmt_ci_pct <- function(mean, lower, upper) {
  paste0(
    fmt_pct(mean),
    "% (95% CrI: ",
    fmt_pct(lower),
    "–",
    fmt_pct(upper),
    "%)"
  )
}

fmt_ci_per10k <- function(mean, lower, upper) {
  paste0(
    fmt_num1(mean),
    " (95% CrI: ",
    fmt_num1(lower),
    "–",
    fmt_num1(upper),
    ")"
  )
}

norm_key2 <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- tolower(x)
  gsub("[^a-z0-9]", "", x)
}

# ----------------------------
# Safety checks
# ----------------------------

required_cols_res <- c(
  "Region",
  "ProvKey",
  "Year",
  "Draw",
  "Cases_WithNets_Resist_who",
  "Cases_WithNets_NoResist_who",
  "Deaths_WithNets_Resist_who",
  "Deaths_WithNets_NoResist_who"
)

missing_cols_res <- setdiff(required_cols_res, names(core_yearly_pw))

if (length(missing_cols_res) > 0) {
  stop(
    "These columns are missing from core_yearly_pw: ",
    paste(missing_cols_res, collapse = ", "),
    "\nRun the full code until the WHO-calibrated *_who variables are created."
  )
}

if (!exists("pop_df")) {
  stop("The object pop_df was not found. Run the population section first.")
}

# ----------------------------
# Study period
# ----------------------------

if (exists("ANALYSIS_YEARS")) {
  resistance_years <- ANALYSIS_YEARS
} else {
  resistance_years <- 2009:2024
}

resistance_years <- resistance_years[
  resistance_years %in% unique(core_yearly_pw$Year)
]

cat("\nResistance-attributable period used: ",
    min(resistance_years), "–", max(resistance_years), "\n", sep = "")

# ============================================================
# 1) NATIONAL RESISTANCE-ATTRIBUTABLE BURDEN
# ============================================================

national_resistance_draws <- core_yearly_pw %>%
  dplyr::filter(Year %in% resistance_years) %>%
  dplyr::group_by(Draw) %>%
  dplyr::summarise(
    Cases_WithNets_Resist_who =
      sum(Cases_WithNets_Resist_who, na.rm = TRUE),
    
    Cases_WithNets_NoResist_who =
      sum(Cases_WithNets_NoResist_who, na.rm = TRUE),
    
    Deaths_WithNets_Resist_who =
      sum(Deaths_WithNets_Resist_who, na.rm = TRUE),
    
    Deaths_WithNets_NoResist_who =
      sum(Deaths_WithNets_NoResist_who, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Additional_cases_due_resistance =
      pmax(Cases_WithNets_Resist_who - Cases_WithNets_NoResist_who, 0),
    
    Percent_case_increase_resistance =
      100 * Additional_cases_due_resistance /
      pmax(Cases_WithNets_NoResist_who, eps),
    
    Deaths_due_resistance =
      pmax(Deaths_WithNets_Resist_who - Deaths_WithNets_NoResist_who, 0)
  )

# National mean population over the study period
national_mean_population <- pop_df %>%
  dplyr::filter(Year %in% resistance_years) %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(
    PopYear = sum(Pop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::summarise(
    MeanNationalPopulation = mean(PopYear, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::pull(MeanNationalPopulation)

national_resistance_draws <- national_resistance_draws %>%
  dplyr::mutate(
    Deaths_due_resistance_per10k =
      Deaths_due_resistance / pmax(national_mean_population / 10000, eps)
  )

national_resistance_summary <- national_resistance_draws %>%
  dplyr::summarise(
    cases_mean =
      mean(Additional_cases_due_resistance, na.rm = TRUE),
    cases_lower95 =
      qL(Additional_cases_due_resistance),
    cases_upper95 =
      qU(Additional_cases_due_resistance),
    
    pct_cases_mean =
      mean(Percent_case_increase_resistance, na.rm = TRUE),
    pct_cases_lower95 =
      qL(Percent_case_increase_resistance),
    pct_cases_upper95 =
      qU(Percent_case_increase_resistance),
    
    deaths_mean =
      mean(Deaths_due_resistance, na.rm = TRUE),
    deaths_lower95 =
      qL(Deaths_due_resistance),
    deaths_upper95 =
      qU(Deaths_due_resistance),
    
    deaths_per10k_mean =
      mean(Deaths_due_resistance_per10k, na.rm = TRUE),
    deaths_per10k_lower95 =
      qL(Deaths_due_resistance_per10k),
    deaths_per10k_upper95 =
      qU(Deaths_due_resistance_per10k)
  ) %>%
  dplyr::mutate(
    cases_95CrI =
      fmt_ci_count(cases_mean, cases_lower95, cases_upper95),
    
    pct_cases_95CrI =
      fmt_ci_pct(pct_cases_mean, pct_cases_lower95, pct_cases_upper95),
    
    deaths_95CrI =
      fmt_ci_count(deaths_mean, deaths_lower95, deaths_upper95),
    
    deaths_k_95CrI =
      fmt_ci_k(deaths_mean, deaths_lower95, deaths_upper95),
    
    deaths_per10k_95CrI =
      fmt_ci_per10k(
        deaths_per10k_mean,
        deaths_per10k_lower95,
        deaths_per10k_upper95
      )
  )

# ============================================================
# 2) PROVINCE RESISTANCE-ATTRIBUTABLE BURDEN
# ============================================================

province_mean_population <- pop_df %>%
  dplyr::filter(Year %in% resistance_years) %>%
  dplyr::group_by(ProvKey = provkey) %>%
  dplyr::summarise(
    MeanProvincePopulation = mean(Pop, na.rm = TRUE),
    .groups = "drop"
  )

province_resistance_draws <- core_yearly_pw %>%
  dplyr::filter(Year %in% resistance_years) %>%
  dplyr::group_by(Region, ProvKey, Draw) %>%
  dplyr::summarise(
    Cases_WithNets_Resist_who =
      sum(Cases_WithNets_Resist_who, na.rm = TRUE),
    
    Cases_WithNets_NoResist_who =
      sum(Cases_WithNets_NoResist_who, na.rm = TRUE),
    
    Deaths_WithNets_Resist_who =
      sum(Deaths_WithNets_Resist_who, na.rm = TRUE),
    
    Deaths_WithNets_NoResist_who =
      sum(Deaths_WithNets_NoResist_who, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  dplyr::left_join(province_mean_population, by = "ProvKey") %>%
  dplyr::mutate(
    Additional_cases_due_resistance =
      pmax(Cases_WithNets_Resist_who - Cases_WithNets_NoResist_who, 0),
    
    Percent_case_increase_resistance =
      100 * Additional_cases_due_resistance /
      pmax(Cases_WithNets_NoResist_who, eps),
    
    Deaths_due_resistance =
      pmax(Deaths_WithNets_Resist_who - Deaths_WithNets_NoResist_who, 0),
    
    Deaths_due_resistance_per10k =
      Deaths_due_resistance /
      pmax(MeanProvincePopulation / 10000, eps)
  )

province_resistance_summary <- province_resistance_draws %>%
  dplyr::group_by(Region, ProvKey) %>%
  dplyr::summarise(
    cases_mean =
      mean(Additional_cases_due_resistance, na.rm = TRUE),
    cases_lower95 =
      qL(Additional_cases_due_resistance),
    cases_upper95 =
      qU(Additional_cases_due_resistance),
    
    pct_cases_mean =
      mean(Percent_case_increase_resistance, na.rm = TRUE),
    pct_cases_lower95 =
      qL(Percent_case_increase_resistance),
    pct_cases_upper95 =
      qU(Percent_case_increase_resistance),
    
    deaths_mean =
      mean(Deaths_due_resistance, na.rm = TRUE),
    deaths_lower95 =
      qL(Deaths_due_resistance),
    deaths_upper95 =
      qU(Deaths_due_resistance),
    
    deaths_per10k_mean =
      mean(Deaths_due_resistance_per10k, na.rm = TRUE),
    deaths_per10k_lower95 =
      qL(Deaths_due_resistance_per10k),
    deaths_per10k_upper95 =
      qU(Deaths_due_resistance_per10k),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Region_key = norm_key2(Region),
    
    cases_95CrI =
      fmt_ci_count(cases_mean, cases_lower95, cases_upper95),
    
    pct_cases_95CrI =
      fmt_ci_pct(pct_cases_mean, pct_cases_lower95, pct_cases_upper95),
    
    deaths_95CrI =
      fmt_ci_count(deaths_mean, deaths_lower95, deaths_upper95),
    
    deaths_k_95CrI =
      fmt_ci_k(deaths_mean, deaths_lower95, deaths_upper95),
    
    deaths_per10k_95CrI =
      fmt_ci_per10k(
        deaths_per10k_mean,
        deaths_per10k_lower95,
        deaths_per10k_upper95
      )
  )

# ============================================================
# 3) EXAMPLE FOUR PROVINCES
# ============================================================

wanted_four_keys <- c(
  "sudkivu",
  "kasaicentral",
  "kinshasa",
  "hautkatanga"
)

four_resistance_summary <- province_resistance_summary %>%
  dplyr::filter(Region_key %in% wanted_four_keys) %>%
  dplyr::arrange(match(Region_key, wanted_four_keys))

# Province with highest percentage case increase among the four provinces
highest_four_pct <- four_resistance_summary %>%
  dplyr::slice_max(
    order_by = pct_cases_mean,
    n = 1,
    with_ties = FALSE
  )

# Extract named province rows
sud_kivu_row <- four_resistance_summary %>%
  dplyr::filter(Region_key == "sudkivu")

kasai_central_row <- four_resistance_summary %>%
  dplyr::filter(Region_key == "kasaicentral")

kinshasa_row <- four_resistance_summary %>%
  dplyr::filter(Region_key == "kinshasa")

haut_katanga_row <- four_resistance_summary %>%
  dplyr::filter(Region_key == "hautkatanga")

# ============================================================
# 4) PRINT VALUES TO REPLACE XXX AND YYY
# ============================================================

cat("\n============================================================\n")
cat("NATIONAL VALUES: PYRETHROID RESISTANCE ATTRIBUTABLE BURDEN\n")
cat("============================================================\n")

cat("National percentage case increase = ",
    national_resistance_summary$pct_cases_95CrI, "\n", sep = "")

cat("XXX national cases due to resistance = ",
    national_resistance_summary$cases_95CrI, " cases\n", sep = "")

cat("National deaths per 10,000 population = ",
    national_resistance_summary$deaths_per10k_95CrI,
    " deaths per 10,000 population\n", sep = "")

cat("YYY national deaths due to resistance = ",
    national_resistance_summary$deaths_95CrI, " deaths\n", sep = "")

cat("YYY national deaths due to resistance in thousands = ",
    national_resistance_summary$deaths_k_95CrI, " deaths\n", sep = "")

cat("\n============================================================\n")
cat("FOUR PROVINCES\n")
cat("============================================================\n")

print(
  four_resistance_summary %>%
    dplyr::select(
      Region,
      pct_cases_95CrI,
      cases_95CrI,
      deaths_95CrI,
      deaths_k_95CrI,
      deaths_per10k_95CrI
    )
)

cat("\nHighest percentage case increase among the four provinces = ",
    highest_four_pct$Region,
    " (",
    highest_four_pct$pct_cases_95CrI,
    ")\n", sep = "")

cat("\nKasai-Central deaths attributed to resistance = ",
    kasai_central_row$deaths_95CrI,
    " deaths\n", sep = "")

cat("Kasai-Central deaths attributed to resistance in thousands = ",
    kasai_central_row$deaths_k_95CrI,
    " deaths\n", sep = "")

cat("\nKinshasa deaths attributed to resistance = ",
    kinshasa_row$deaths_k_95CrI,
    " deaths\n", sep = "")

cat("Haut-Katanga deaths attributed to resistance = ",
    haut_katanga_row$deaths_k_95CrI,
    " deaths\n", sep = "")

# ============================================================
# 5) MANUSCRIPT SENTENCE
# ============================================================

sentence_resistance <- paste0(
  "Predicted malaria cases and deaths attributable to pyrethroid resistance varied ",
  "across national estimates and the four example provinces. Nationally, resistance ",
  "was associated with a mean increase of ",
  fmt_pct(national_resistance_summary$pct_cases_mean),
  "\\% in malaria cases, equivalent to ",
  national_resistance_summary$cases_95CrI,
  " additional cases, and ",
  national_resistance_summary$deaths_per10k_95CrI,
  " deaths per 10,000 population, corresponding to ",
  national_resistance_summary$deaths_95CrI,
  " deaths over the study period. Among the four example provinces, ",
  highest_four_pct$Region,
  " had the highest modelled percentage increase in cases due to resistance, ",
  "whereas Kasaï-Central had the largest absolute mortality burden, with approximately ",
  kasai_central_row$deaths_k_95CrI,
  " deaths attributed to pyrethroid resistance. Kinshasa and Haut-Katanga had lower ",
  "resistance-attributable mortality burdens, estimated at approximately ",
  kinshasa_row$deaths_k_95CrI,
  " and ",
  haut_katanga_row$deaths_k_95CrI,
  " deaths, respectively."
)

cat("\n============================================================\n")
cat("MANUSCRIPT SENTENCE\n")
cat("============================================================\n")
cat(sentence_resistance, "\n")

# Optional: save table
write.csv(
  four_resistance_summary,
  "four_provinces_resistance_attributable_burden.csv",
  row.names = FALSE
)

write.csv(
  national_resistance_summary,
  "national_resistance_attributable_burden.csv",
  row.names = FALSE
)
# -----------------------------------------------------------------------------
# =========================================================
# PRINT ADDITIONAL CASES AND DEATHS DUE TO RESISTANCE
# =========================================================

library(dplyr)
library(scales)

# ---------------------------------------------------------
# 1. Choose period
# ---------------------------------------------------------
# For the full study period:
analysis_years <- 2009:2024

# If you want only one year, for example 2024, use:
# analysis_years <- 2024

# ---------------------------------------------------------
# 2. Choose WHO-calibrated or raw simulation outputs
# ---------------------------------------------------------
# TRUE  = use WHO-calibrated cases and deaths
# FALSE = use raw simulation outputs

use_who_calibrated <- TRUE

# ---------------------------------------------------------
# 3. Select correct columns
# ---------------------------------------------------------

if (use_who_calibrated) {
  
  cases_res_col    <- "Sum_WithNets_Resist_who"
  cases_nores_col  <- "Sum_WithNets_NoResist_who"
  deaths_res_col   <- "Sum_Deaths_WithNets_Resist_who"
  deaths_nores_col <- "Sum_Deaths_WithNets_NoResist_who"
  
} else {
  
  cases_res_col    <- "Sum_WithNets_Resist_s"
  cases_nores_col  <- "Sum_WithNets_NoResist_s"
  deaths_res_col   <- "Sum_Deaths_WithNets_Resist_s"
  deaths_nores_col <- "Sum_Deaths_WithNets_NoResist_s"
}

required_cols <- c(
  "Year", "Draw",
  cases_res_col, cases_nores_col,
  deaths_res_col, deaths_nores_col
)

missing_cols <- setdiff(required_cols, names(national_yearly_pw))

if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing from national_yearly_pw:\n",
    paste(missing_cols, collapse = ", "),
    "\n\nAvailable columns are:\n",
    paste(names(national_yearly_pw), collapse = ", ")
  )
}

# ---------------------------------------------------------
# 4. Calculate additional cases and deaths due to resistance
# ---------------------------------------------------------

resistance_burden_yearly <- national_yearly_pw %>%
  filter(Year %in% analysis_years) %>%
  mutate(
    cases_with_resistance =
      as.numeric(.data[[cases_res_col]]),
    
    cases_without_resistance =
      as.numeric(.data[[cases_nores_col]]),
    
    deaths_with_resistance =
      as.numeric(.data[[deaths_res_col]]),
    
    deaths_without_resistance =
      as.numeric(.data[[deaths_nores_col]]),
    
    additional_cases_due_to_resistance =
      cases_with_resistance - cases_without_resistance,
    
    additional_deaths_due_to_resistance =
      deaths_with_resistance - deaths_without_resistance
  )

# ---------------------------------------------------------
# 5. Summarise by posterior draw
# ---------------------------------------------------------

resistance_burden_by_draw <- resistance_burden_yearly %>%
  group_by(Draw) %>%
  summarise(
    total_additional_cases =
      sum(additional_cases_due_to_resistance, na.rm = TRUE),
    
    total_additional_deaths =
      sum(additional_deaths_due_to_resistance, na.rm = TRUE),
    
    .groups = "drop"
  )

# ---------------------------------------------------------
# 6. National summary with 95% credible intervals
# ---------------------------------------------------------

resistance_burden_summary <- resistance_burden_by_draw %>%
  summarise(
    mean_additional_cases =
      mean(total_additional_cases, na.rm = TRUE),
    
    lower_cases =
      quantile(total_additional_cases, 0.025, na.rm = TRUE),
    
    upper_cases =
      quantile(total_additional_cases, 0.975, na.rm = TRUE),
    
    mean_additional_deaths =
      mean(total_additional_deaths, na.rm = TRUE),
    
    lower_deaths =
      quantile(total_additional_deaths, 0.025, na.rm = TRUE),
    
    upper_deaths =
      quantile(total_additional_deaths, 0.975, na.rm = TRUE)
  )

print(resistance_burden_summary)

# ---------------------------------------------------------
# 7. Values to insert into the manuscript
# ---------------------------------------------------------

XXX_cases <- round(resistance_burden_summary$mean_additional_cases)
YYY_deaths <- round(resistance_burden_summary$mean_additional_deaths)

cat("\nValues to insert in the manuscript:\n")
cat("XXX =", comma(XXX_cases), "clinical cases\n")
cat("YYY =", comma(YYY_deaths), "deaths\n\n")

# ---------------------------------------------------------
# 8. Print manuscript sentence
# ---------------------------------------------------------

sentence_plain <- sprintf(
  paste0(
    "Results indicate that resistance has substantially undermined the protective effect ",
    "of pyrethroid-only nets, nationally leading to an additional %s clinical cases ",
    "and %s deaths due to the disease."
  ),
  comma(XXX_cases),
  comma(YYY_deaths)
)

cat(sentence_plain, "\n\n")

# ---------------------------------------------------------
# 9. Print Overleaf sentence with uncertainty
# ---------------------------------------------------------

sentence_latex <- sprintf(
  paste0(
    "Results indicate that resistance has substantially undermined the protective effect ",
    "of pyrethroid-only nets, nationally leading to an additional %s clinical cases ",
    "(95\\%% CrI: %s--%s) and %s deaths ",
    "(95\\%% CrI: %s--%s) due to the disease."
  ),
  comma(round(resistance_burden_summary$mean_additional_cases)),
  comma(round(resistance_burden_summary$lower_cases)),
  comma(round(resistance_burden_summary$upper_cases)),
  comma(round(resistance_burden_summary$mean_additional_deaths)),
  comma(round(resistance_burden_summary$lower_deaths)),
  comma(round(resistance_burden_summary$upper_deaths))
)

cat("Overleaf sentence:\n")
cat(sentence_latex, "\n")
###################### -----------------------------------------------------------------------------

# =========================================================
# PRINT % INCREASE IN CASES AND DEATHS DUE TO RESISTANCE
# =========================================================

library(dplyr)
library(scales)

# Choose period
analysis_years <- 2009:2024
# For only 2024, use:
# analysis_years <- 2024

# Use WHO-calibrated columns
cases_res_col    <- "Sum_WithNets_Resist_who"
cases_nores_col  <- "Sum_WithNets_NoResist_who"
deaths_res_col   <- "Sum_Deaths_WithNets_Resist_who"
deaths_nores_col <- "Sum_Deaths_WithNets_NoResist_who"

# Check columns
required_cols <- c(
  "Year", "Draw",
  cases_res_col, cases_nores_col,
  deaths_res_col, deaths_nores_col
)

missing_cols <- setdiff(required_cols, names(national_yearly_pw))

if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing:\n",
    paste(missing_cols, collapse = ", ")
  )
}

# Calculate additional cases/deaths and percentage increase
resistance_percent_by_draw <- national_yearly_pw %>%
  filter(Year %in% analysis_years) %>%
  mutate(
    cases_with_resistance = as.numeric(.data[[cases_res_col]]),
    cases_without_resistance = as.numeric(.data[[cases_nores_col]]),
    deaths_with_resistance = as.numeric(.data[[deaths_res_col]]),
    deaths_without_resistance = as.numeric(.data[[deaths_nores_col]])
  ) %>%
  group_by(Draw) %>%
  summarise(
    total_cases_with_resistance = sum(cases_with_resistance, na.rm = TRUE),
    total_cases_without_resistance = sum(cases_without_resistance, na.rm = TRUE),
    total_deaths_with_resistance = sum(deaths_with_resistance, na.rm = TRUE),
    total_deaths_without_resistance = sum(deaths_without_resistance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    percent_increase_cases =
      100 * (total_cases_with_resistance - total_cases_without_resistance) /
      total_cases_without_resistance,
    
    percent_increase_deaths =
      100 * (total_deaths_with_resistance - total_deaths_without_resistance) /
      total_deaths_without_resistance
  )

# Summary with 95% CrI
resistance_percent_summary <- resistance_percent_by_draw %>%
  summarise(
    mean_percent_cases = mean(percent_increase_cases, na.rm = TRUE),
    lower_percent_cases = quantile(percent_increase_cases, 0.025, na.rm = TRUE),
    upper_percent_cases = quantile(percent_increase_cases, 0.975, na.rm = TRUE),
    
    mean_percent_deaths = mean(percent_increase_deaths, na.rm = TRUE),
    lower_percent_deaths = quantile(percent_increase_deaths, 0.025, na.rm = TRUE),
    upper_percent_deaths = quantile(percent_increase_deaths, 0.975, na.rm = TRUE)
  )

print(resistance_percent_summary)

# Values to put in the text
XXX_percent_cases <- resistance_percent_summary$mean_percent_cases
YYY_percent_deaths <- resistance_percent_summary$mean_percent_deaths

cat("\nValues to insert:\n")
cat("XXX =", round(XXX_percent_cases, 1), "% increase in clinical cases\n")
cat("YYY =", round(YYY_percent_deaths, 1), "% increase in deaths\n\n")

# Manuscript sentence
sentence_plain <- sprintf(
  paste0(
    "Results indicate that resistance has substantially undermined the protective effect ",
    "of pyrethroid-only nets, nationally leading to a %.1f%% increase in clinical cases ",
    "and a %.1f%% increase in deaths due to the disease."
  ),
  XXX_percent_cases,
  YYY_percent_deaths
)

cat(sentence_plain, "\n\n")

# Overleaf sentence with uncertainty
sentence_latex <- sprintf(
  paste0(
    "Results indicate that resistance has substantially undermined the protective effect ",
    "of pyrethroid-only nets, nationally leading to a %.1f\\%% increase in clinical cases ",
    "(95\\%% CrI: %.1f--%.1f) and a %.1f\\%% increase in deaths ",
    "(95\\%% CrI: %.1f--%.1f) due to the disease."
  ),
  resistance_percent_summary$mean_percent_cases,
  resistance_percent_summary$lower_percent_cases,
  resistance_percent_summary$upper_percent_cases,
  resistance_percent_summary$mean_percent_deaths,
  resistance_percent_summary$lower_percent_deaths,
  resistance_percent_summary$upper_percent_deaths
)

cat("Overleaf sentence:\n")
cat(sentence_latex, "\n")


# =============================================================================
# End of script
# =============================================================================
