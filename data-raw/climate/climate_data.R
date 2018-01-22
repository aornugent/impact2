#' Pre-processing script for rainfal data
#'
#' Data obtained from Australian Government Bureau of Meteorology
#'

# Load data
rain_1 <- readr::read_csv("data-raw/climate/IDCJAC0009_070093_1800_Data.csv")
rain_2 <- readr::read_csv("data-raw/climate/IDCJAC0009_070242_1800_Data.csv")
rain_3 <- readr::read_csv("data-raw/climate/IDCJAC0009_070277_1800_Data.csv")

stations <- readr::read_csv("data-raw/climate/weather_stations.csv") %>%
  dplyr::mutate(station_id = as.integer(`Bureau of Meteorology station number`))

# Join data
rain <- dplyr::bind_rows(rain_1, rain_2, rain_3) %>%
  magrittr::set_colnames(
    c("bom_code", "station_id", "year",
      "month", "day", "rainfall",
      "collection_period", "quality")) %>%
  dplyr::mutate(
    station_id = as.integer(station_id),
    month = as.integer(month),
    date = lubridate::dmy(paste(day, month, year)),
    unit = "mm") %>%
  dplyr::left_join(stations)

# Summarise spring rainfall
spring_rain <- dplyr::filter(rain, year %in% c(2010:2016),
                             month %in% c(8:11),
                             grepl("MELBA|ARANDA", `Station name`)) %>%
  dplyr::group_by(year, station_id) %>%
  dplyr::summarise(spring_total = sum(rainfall, na.rm = T)) %>%
  dplyr::summarise(rain_mm = round(mean(spring_total), 0))

write_csv(x = spring_rain, path = "data-raw/pinnacle_rainfall.csv")

rm(list = ls())
