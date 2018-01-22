#' Combine data sets and save
#'
#' This script concatenates the Pinnacle survey data with the appropriate
#' environmental covariates and includes the resulting dataframe as a package
#' dataset.
#'
#' The survey and climate folders each contain pre-processing scripts to produce
#' the .csv files found here.
#'
#' Cover data does not contain absences and is aggregated as the mean cover
#' value for each species across all quadrats in a plot.
#'
#' @example
#'
#' source("data-raw/survey/survey_data.R")
#' source("data-raw/climate/climate_data.R")

# Soil and treatment data
env <- readr::read_csv("data-raw/pinnacle_covariates.csv") %>%
  dplyr::mutate(
    treatment = dplyr::case_when(
      grepl("control", treatment) ~ "Control",
      grepl("slash", treatment) ~ "Slashed",
      TRUE ~ treatment),
    fence = dplyr::case_when(
      grepl("fenced", fence) ~ "Fenced",
      grepl("open", fence) ~ "Grazed")) %>%
  dplyr::select(plot_id, treatment, fence, totalN_ppm)


# Rainfall data
rain <- readr::read_csv("data-raw/pinnacle_rainfall.csv")


# Cover data
cover <- readr::read_csv("data-raw/pinnacle_surveys.csv",
                         col_type = readr::cols(
                           cover = readr::col_double()))

# Join and filter
cover <- cover %>%
  dplyr::left_join(env) %>%
  dplyr::left_join(rain) %>%
  dplyr::filter(grepl("Control|Slashed", treatment))

# Aggregate to plot level
cover <- dplyr::group_by(cover, year, site, plot_id, species) %>%
  dplyr::mutate(cover = mean(cover)) %>%
  dplyr::select(-quadrat) %>%
  dplyr::distinct()

# Trait data
traits <- readr::read_csv("data-raw/pinnacle_traits.csv")


# Save datasets
devtools::use_data(cover, overwrite = T)
devtools::use_data(traits, overwrite = T)

rm(list = ls())
