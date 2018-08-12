#' Run models
#'
#' This analysis compares several models of increasing complexity. Models
#' are coded in Stan and can be found in the `/models` subdirectory.
#'
#' Model output is saved to .Rdata files for later inspection. These are
#' ~0.4-0.9Gb and are not distributed. Previously saved analysis can be
#' inspected using \code{check_models()} or loaded with \code{load_models()}.
#'
#' Models include:
#' \enumerate{
#'  \item{M0: Logistic regression of the proportion of nonnative species cover}
#'  \item{M1: Joint species model with varying intercepts and slopes, constant covariance}
#'  \item{M2: Joint species model with constant intercepts and slopes, varying covariance}
#'  \item{M3: Joint species model with varying intercepts, slopes and covariance}
#' }
#'
#' M0, M1, and M3 are presented in the manuscript (O'Reilly-Nugent, et al. 2018)
#'
#' @param model m0-m3, defaults to all (~2hrs).
#' @param path  Directory to save model outputs, defaults to /models.
#' @param check auto-run model checks, boolean (default = T).
#' @param ...   passes additional parameters to Stan models such as the LKJ shape prior and abundance threshold of the analysis.
#'
#' @usage run_model(model = "m1")
#'
#' @export

run_models <- function(
  models = c("m0", "m1", "m2", "m3"),
  path = "models/",
  check = T, ...) {

  for(model in models) {
    model_output <- run_stan_model(model, path, ...)

    if(check == T)
      check_models(model, path, model_output)
  }
}


#' Run Stan model
#'
#' Pre-processes data and wraps \link{rstan} to run analyses using
#' adaptive Hamiltonian Monte Carlo. Sensible default settings are provided,
#' but can be tweaked as necessary.
#'
#' @param model one of c("m0", "m1", "m2", "m3")
#' @param ... tuning parameters for Stan sampler
#'
#' @export

run_stan_model <- function(model, path, ...) {

  # Set parallel options
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores() - 1)

  # Do data prep
  data_list <- format_data(model, ...)

  # Get appropriate file
  model_file <- switch(
    model,
    m0 = "models/m0_proportional_logistic_regression.stan",
    m1 = "models/m1_interacting_species_tobit_varying_intercepts_slopes.stan",
    m2 = "models/m2_interacting_species_tobit_varying_interactions.stan",
    m3 = "models/m3_interacting_species_tobit_varying_intercepts_slopes_covariance.stan"
  )

  # Pre-compile model
  printf("\nCompiling model... (ignore non-fatal warnings)")

  suppressWarnings(
    testfit <- rstan::stan(
      file = model_file,
      data = data_list,
      iter = 1,
      chains = 1,
      init_r = 0.5)
    )

  # Run model
  printf("Running model...")

  mod <- rstan::stan(
    fit = testfit,
    data = data_list,
    iter = 2000,
    chains = 3,
    init_r = 0.5,
    save_warmup = F,
    control = list(max_treedepth = 15,
                   adapt_delta = 0.8)
  )

  summary <- rstan::summary(mod)$summary %>%
    as.data.frame()

  # Append data_list to model output
  model_output = list(model_summary = summary,
                      stan_output = mod,
                      data_list = data_list)

  # Save model output
  filename = paste0(path, model, "_output.Rdata")
  printf(c("Saving", filename))
  save(model_output, file = filename)

  return(model_output)
}


#' Logistic regression
#'
#' Run logit regression in JAGS
#'
#' @param model "m0" only
#' @param path path to jags model (defaults to "models/")
#'
#' @export

logistic_regression <- function(model = "m0", path = "models/"){

  m <- list()
  m$data_list <- format_data(model)

  pars <- c("B_int", "B_slope", "B_rain", "B_plot", "B_site",
            "mu_int", "mu_slope","sigma_plot", "sigma_site",
            "sigma_int", "sigma_slope", "sigma")

  m$jags_output <- jagsUI::jags(model = paste0(path, "logit_prop.jags"),
                                data = m$data_list[1:15],
                                parameters.to.save = pars,
                                n.iter = 50000,
                                n.burnin = 20000,
                                n.thin = 30,
                                n.chains = 4,
                                parallel = T)


  m$model_summary <- as.data.frame(m$jags_output$summary)

  filename <- paste0(path, model, "_output.Rdata")
  save(m, file = filename)
}


#' Format data for analysis
#'
#' Takes the preloaded Pinnacle vegetation dataset and prepares it for
#' analysis in Stan. The Stan models are written such that the are able
#' to take largely the same data structure, with the exception of M0.
#'
#' @param years restrict analysis by years.
#' @param threshold restrict analysis to some proportion of total abundance.
#' @param subset subset species by presence or abundance.
#' @param lkj_prior shape value for covariance matrices in joint models.

format_data <- function(model,
                        years = c(2013:2016),
                        threshold = 0.2,
                        subset = "presence",
                        lkj_prior = 25, ...) {

  # Load dataset
  data(cover)

  if(model == "m0")
    years = 2010:2016

  # Restrict by years, drop unidentified species.
  dat <- filter(cover,
           year %in% years,
           !grepl("Unidentified", species)) %>%
         ungroup() %>%
    mutate(quadrat_id = as.numeric(factor(paste(plot_id, quadrat))))

  # Subset for common species
  if(subset == "presence") {
    species_list <- subset_species_plots(dat, threshold)
  }

  if(subset == "abundance") {
    species_list <- subset_species_abundance(dat, threshold)
  }


  # Extract environmental covariates
  env_covariates <-
    select(
      dat,
      year,
      site,
      plot_id,
      quadrat_id,
      fence,
      treatment,
      totalN_ppm,
      rain_mm)  %>%
    unique() %>%
    arrange(year, plot_id, quadrat_id) %>%
    mutate(
      rain_scaled = scale(rain_mm),
      totalN_scaled = scale(totalN_ppm),
      year_id = as.numeric(factor(year)),
      treatment_id = as.numeric(factor(paste(fence, treatment)))
    )

  if(model == "m0"){
    # Format data for logistic regression of proportion of exotic species.
    env_covariates <- select(env_covariates, -quadrat_id) %>%
      unique()

    # Calculate proportion of exotic species
    prop <- group_by(dat, year, plot_id, introduced) %>%
      summarise(group_cover = sum(cover)) %>%
      mutate(prop_cover = group_cover/sum(group_cover)) %>%
      ungroup() %>%
      filter(introduced == 1)

    # Scaled response: 0 < y < 1, Smithson (2006)
    logit <- mutate(prop,
                     prop_scaled = (prop_cover * (n() - 1) + 0.5)/n(),
                     logit = log(prop_scaled /(1 - prop_scaled)))

    y = logit$logit

  }  else {
    # Format data for tobit regression of cover values.

    cover <- group_by(dat,
                      year,
                      plot_id,
                      quadrat_id,
                      species) %>%
      summarise(cover = mean(cover))

    # Reshape multivariate data
    cover_wide <- filter(cover,
                         species %in% unlist(species_list$species)) %>%
                  tidyr::spread(species, cover, fill = 0)

    # Drop unlabelled columns
    y = as.matrix(cover_wide[, c(-1, -2, -3)])
  }

  # Set up model matrix with intercept, fertility and rainfall (scaled)
  X <- model.matrix(~ 1 + totalN_scaled + rain_scaled, env_covariates)

  # Nest quadrats within plots, within sites
  plots <- as.numeric(as.factor(env_covariates$plot_id))
  sites <- nested(env_covariates$plot_id, env_covariates$site)

  data_list <-
    list(
      y_observed = y,
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
      N_years = max(env_covariates$year_id),
      shape_prior = lkj_prior,
      species_list = species_list,
      env_covariates = env_covariates)

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

subset_species_abundance <- function(dat, threshold) {

  # Aggregate abundance by species, then rank and filter.
  group_by(dat, species, introduced) %>%
    summarise(total_abun = sum(cover)) %>%
    ungroup() %>%
    mutate(prop_abun = total_abun / sum(total_abun)) %>%
    arrange(desc(prop_abun)) %>%
    filter(cumsum(prop_abun) < threshold) %>%
    mutate(species_id = as.numeric(as.factor(species))) %>%
    arrange(species_id)
}


#' Generate proportional presence lists
#'
#' Selects the most abundant species in the dataset by the proportion plots
#' a species was observed in. Discards species absent in one or more years.
#'
#' 15% of plots is a typical default.
#'
#' @usage subset_species(cover, 0.15)
#'
#' @importFrom tidyr gather

subset_species_plots <- function(dat, threshold) {

  # Get total number of plots
  n_total = n_distinct(paste(dat$year, dat$plot_id))

  # Aggregate presences by species, filter by threshold.
  group_by(dat, introduced, species, year) %>%
    summarise(n = n_distinct(plot_id)) %>%
    spread(year, n, fill = 0) %>%
    gather(year, n, -species, -introduced) %>%
    summarise(total_plots = sum(n),
              prop_plots = total_plots / n_total) %>%
    filter(prop_plots >= threshold) %>%
    ungroup() %>%
    mutate(species_id = as.numeric(as.factor(species))) %>%
    arrange(species_id)
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

