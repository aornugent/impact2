#' Pre-processing script for Pinnacle surveys
#'
#' Data-sheet format varies by years:
#'   2010 - Percent cover of each 1 x 1m quadrat, with plots as columns
#'   2011-2013 - As 2010 plus 0.3 x 0.3m sub-quadrat counts, 9 total
#'   2015-2016 - Presence absence and percent cover, sub-quadrats as columns


# Convert all files into long format, select cover only
cover_10 <- readr::read_csv("data-raw/survey/pinnacle_2010_quadrats.csv") %>%
  tidyr::gather(plot, cover, `1_pc`:`10_pc`) %>%
  tidyr::separate(plot, c("plot", "measure"), sep = "_", convert = T) %>%
  dplyr::filter(measure == "pc") %>%
  dplyr::select(date, site, plot, quadrat, id, cover)

cover_11 <- readr::read_csv("data-raw/survey/pinnacle_2011_quadrats.csv") %>%
  tidyr::gather(plot, cover, `1_pc`:`10_9`) %>%
  tidyr::separate(plot, c("plot", "measure"), sep = "_", convert = T) %>%
  dplyr::filter(measure == "pc") %>%
  dplyr::select(date, site, plot, quadrat, id, cover)

cover_12 <- readr::read_csv("data-raw/survey/pinnacle_2012_quadrats.csv") %>%
  tidyr::gather(plot, cover, `1_pc`:`10_9`) %>%
  tidyr::separate(plot, c("plot", "measure"), sep = "_", convert = T) %>%
  dplyr::filter(measure == "pc") %>%
  dplyr::select(date, site, plot, quadrat, id, cover)

cover_13 <- readr::read_csv("data-raw/survey/pinnacle_2013_quadrats.csv") %>%
  tidyr::gather(plot, cover, `1_pc`:`10_9`) %>%
  tidyr::separate(plot, c("plot", "measure"), sep = "_", convert = T) %>%
  dplyr::filter(measure == "pc") %>%
  dplyr::select(date, site, plot, quadrat, id, cover)

cover_15 <- readr::read_csv("data-raw/survey/pinnacle_2015_quadrats.csv") %>%
  dplyr::select(date, site, plot, quadrat, id, cover)

cover_16 <- readr::read_csv("data-raw/survey/pinnacle_2016_quadrats.csv") %>%
  dplyr::select(date, site, plot, quadrat, id, cover)

# Bind together and filter out absences
cover <- bind_rows(cover_10, cover_11, cover_12, cover_13, cover_15, cover_16) %>%
  dplyr::filter(cover > 0) %>%
  dplyr::mutate(year = as.integer(format(as.Date(date, format="%d/%m/%Y"),"%Y")),
         site = tolower(site),
         plot_id = paste(site, plot, sep = "_"))


# Remove intermediary data frames
rm(list = grep("_", ls(), value = T))


# Resolve species names
species_id <- readr::read_csv("data-raw/survey/id_to_species.csv")

# Aggregate species complexes
cover <- dplyr::left_join(cover, species_id) %>%
  dplyr::mutate(
    species = dplyr::case_when(
      grepl("Vulpia",  species) ~ "Vulpia.sp",
      grepl("Rytidosperma", species) ~ "Rytidosperma.sp",
      grepl("Aira", species) ~ "Aira.sp",
      TRUE ~ as.character(species))) %>%
  dplyr::group_by(year, site, plot, plot_id, quadrat, species) %>%
  dplyr::summarise(cover = ifelse(sum(cover) > 100, 100, sum(cover))) %>%
  dplyr::ungroup()

# Remove non-plant IDs
cover <- dplyr::filter(cover, !grepl("Bare.ground|Rock|Log|Feather|Scat|Litter", species))

# Remove one maple tree seedling
cover <- dplyr::filter(cover, !grepl("Acer.negundo", species))

# Join species origin
tax <- readr::read_csv("data-raw/survey/pinnacle_taxonomy.csv") %>%
  dplyr::select(species, introduced)

cover <- dplyr::left_join(cover, tax)

write_csv(x = cover, path = "data-raw/pinnacle_surveys.csv")

# Do the same for traits

traits <- readr::read_csv("data-raw/survey/pinnacle_traits_2016.csv") %>%
  dplyr::filter(species %in% species_id$species) %>%
  dplyr::mutate(
    species = dplyr::case_when(
      grepl("Vulpia",  species) ~ "Vulpia.sp",
      grepl("Rytidosperma", species) ~ "Rytidosperma.sp",
      grepl("Aira", species) ~ "Aira.sp",
      TRUE ~ as.character(species))) %>%
  dplyr::group_by(species, trait) %>%
  dplyr::summarise(mean = mean(mean),
                   max = mean(max))

write_csv(x = traits, path = "data-raw/pinnacle_traits.csv")


rm(list = ls())
