#' Run models
#'
#' This analysis compares several models of increasing complexity. Models
#' are coded in Stan and can be found in the `/models` subdirectory.
#'
#' Model output is saved to .Rdata files for later inspection. These are
#' ~0.4-0.9Gb and are not distributed. Previously saved analysis can be
#' inspected using \code{check_models()} or loaded with \code{load_models()}.
#'
#' @param model M0-M6, defaults to all (~2hrs).
#' @param path  Directory to save model outputs, defaults to /models.
#'
#' @usage run_model(model = "m1")
#'
#' @export

run_models <- function(
  models = c("m0", "m1", "m2", "m3", "m4", "m5", "m6"),
  dir = "models/", ...) {

  for(model in models) {
    run_stan_model(model, dir, ...)
  }
}


#' Run Stan model
#'
#' Contains pre-processing and wraps \link[rstan] to run analyses using
#' adaptive Hamiltonian Monte Carlo. Sensible default settings are provided,
#' but can be tweaked as necessary.
#'
#' @param model one of c("m0", "m1", "m2", "m3", "m4", "m5", "m6")
#' @param ... tuning parameters for Stan sampler
#'
#' @export

run_stan_model <- function(model, dir, ...) {

  # Set parallel options
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores() - 1)

  # FIGURE OUT ...
  # if(cores = NULL){
  #   options(mc.cores = parallel::detectCores() - 1)
  # }

  # Do data prep
  data_list <- format_data(model, ...)

  # FIGURE OUT ...
  # Pass number of years, currently hard coded

  # Get appropriate file
  model_file <- switch(
    model,
    m0 = "models/m0_prop_logit",
    m1 = "models/m1_independent_species_tobit.stan",
    m2 = "models/m2_independent_species_tobit_varying_intercepts_slopes.stan",
    m3 = "models/m3_interacting_species_tobit.stan",
    m4 = "models/m4_interacting_species_tobit_varying_intercepts_slopes.stan",
    m5 = "models/m5_interacting_species_tobit_varying_interactions.stan",
    m6 = "models/m6_interacting_species_tobit_varying_intercepts_slopes_covariance.stan"
  )

  # Pre-compile model
  printf("\nCompiling model...")

  suppressWarnings(
    testfit <- rstan::stan(
      file = model_file,
      data = data_list,
      iter = 1,
      chains = 1,
      init_r = 1)
    )

  # Run model
  printf("Running model...")

  mod <- rstan::stan(
    fit = testfit,
    data = data_list,
    iter = 2000,
    chains = 3,
    init_r = 1,
    save_warmup = F
  )

  # FIGURE OUT ...
  # rstan::stan(control = ...)

  # Append data_list to model output
  model_output = list(stan_output = mod,
                      data_list = data_list)

  # Save model output
  filename = paste0(dir, model, "_output.Rdata")
  printf(c("Saving", filename))
  save(model_output, file = filename)

  return(NULL)
}


#' Format data for analysis
#'
#' Takes the preloaded Pinnacle vegetation dataset and prepares it for
#' analysis in Stan. The Stan models are written such that the are able
#' to take largely the same data structure, with the exception of M0.
#'
#' @param threshold restrict analysis to some proportion of total abundance.
#' @param years restrict analysis by years.
#' @param lkj_prior shape value for covariance matrices in joint models.

format_data <- function(model,
                        threshold = 0.9,
                        years = c(2013:2016),
                        lkj_prior = 25, ...) {

  # Load dataset
  data(cover)

  if(model == "m0")
    years = 2010:2016

  # Restrict by years, drop unidentified species.
  dat <- dplyr::filter(cover,
                       year %in% years,
                       !grepl("Unidentified", species)) %>%
         dplyr::ungroup()

  # Subset for common species
  species_list <- subset_species(dat, threshold)

  # Extract environmental covariates
  env_covariates <- dplyr::select(dat,
                                  year,
                                  site,
                                  plot_id,
                                  fence,
                                  treatment,
                                  totalN_ppm,
                                  rain_mm)  %>%
                    dplyr::distinct() %>%
                    dplyr::mutate(
                      rain_scaled = scale(rain_mm),
                      totalN_scaled = scale(totalN_ppm),
                      year_id = as.numeric(factor(year)),
                      treatment_id = as.numeric(factor(paste(fence, treatment)))
                    )

  if(model == "m0"){
    # Format data for logistic regression of proportion of exotic species.

    # Calculate proportion of exotic species
    prop <- dplyr::group_by(dat, year, plot_id, introduced) %>%
      dplyr::summarise(group_cover = sum(cover)) %>%
      dplyr::mutate(prop_cover = group_cover/sum(group_cover)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(introduced == 1)

    # Scaled response: 0 < y < 1, Smithson (2006)
    logit <- dplyr::mutate(prop,
                           prop_scaled = (prop_cover * (n() - 1) + 0.5)/n(),
                           logit = log(prop_scaled /(1 - prop_scaled)))

    y = logit$logit
  }  else {
    # Format data for tobit regression of cover values.

    # Reshape multivariate data
    cover_wide <- dplyr::select(dat, year, plot_id, species, cover) %>%
                  dplyr::filter(species %in% species_list$species) %>%
                  tidyr::spread(species, cover, fill = 0)

    y = as.matrix(cover_wide[, c(-1, -2)])
  }

  # Set up model matrix with intercept, fertility and rainfall (scaled)
  X <- model.matrix(~ 1 + totalN_scaled + rain_scaled, env_covariates)

  # Nest plots within sites
  plots <- as.numeric(factor(env_covariates$plot_id))
  sites <- nested(env_covariates$plot_id, env_covariates$site)

  data_list <- list(y_observed = y,
                    N = nrow(env_covariates),
                    N_censored = sum(y == 0),
                    S = nrow(species_list),
                    K = ncol(X),
                    X = X,
                    E = max(env_covariates$treatment_id),
                    treatment = env_covariates$treatment_id,
                    plot = plots,
                    site = sites,
                    year = env_covariates$year_id,
                    N_plots = max(plots),
                    N_sites = max(sites),
                    N_years = max(env_covariates$year),
                    shape_prior = lkj_prior,
                    species_list = species_list)

  return(data_list)
}


#' Generate proportional abundance lists
#'
#' Selects the most abundant species in the dataset up to some poportional
#' threshold. Species are ordered alphabetically and given numeric ids.
#'
#' 90% of observed cover is a typical default.
#'
#' @usage subset_species(cover, 0.9)

subset_species <- function(dat, threshold) {
  dplyr::group_by(dat, species) %>%
  dplyr::summarise(total = sum(cover)) %>%
  dplyr::mutate(prop = total / sum(total)) %>%
  dplyr::arrange(desc(prop)) %>%
  dplyr::mutate(common = ifelse(cumsum(prop) < threshold, 1, 0)) %>%
  dplyr::filter(common == 1) %>%
  dplyr::select(species) %>%
  dplyr::mutate(species_id = as.numeric(as.factor(species))) %>%
  dplyr::arrange(species_id)
}

#' Index nested factors
#'
#' Used for double-indexing in Stan models. Assigns numerical id of top-level
#' factors to each entry of lower-level factors.
#'
#' @usage nested(plot, site)


nested <- function(x, y){
  nestbl <- table(factor(x), factor(y))
  idx <- unname(apply(nestbl, 1, function(x) which(x > 0)))
  return(idx)
}

