#' Create figures
#'
#' Creates a series of figures to aide the interpretation of model outputs.
#' The figures are model dependent and should only be used if the models have
#' converged correctly.
#'
#' Figures include:
#' \enumerate{
#'  \item{a highest-density interval summary of logistic regression (m0)}
#'  \item{linear tobit predictions of species cover along a fertility gradient
#'  (m1-6)}
#'  \item{heatmaps of between species covariances (m3-6)}
#'  \item{network diagrams of negative interactions between species (m3-6)}
#' }
#'
#' @param model m0-m6, defaults to all.
#' @param path  Directory to save model outputs, defaults to /models.
#' @param figs  numeric: 1-4, defaults to all.
#'
#' @usage create_figures_models(model = "m0", fig = 1)
#' @export

create_figures <-  function(
  models = c("m0", "m1", "m2", "m3", "m4", "m5", "m6"),
  path = "models/",
  figs = 1:4, ...) {

  for(fig in figs){
    for(model in models) {

      # Create HDI figure for logistic regression
      if(fig == 1 & model == "m0")
        hdi_interval_figure(path)

      # Create linear predictions for tobit models
      if(fig == 2 & !grepl("m0", model))
        linear_tobit_figure(model, path)

      # Create covariance heatmap for joint tobit models
      if(fig == 3 & !grepl("m0", model))
        covariance_heatmap_figure(model, path)

      # Create interaction network diagram for joint tobit models
      if(fig == 4 & !grepl("m0", model))
        interaction_network_figure(model, path)

      # Create interaction density figure for models with varying covariances
      if(fig == 5)
        covariance_density_figure(model, path)

      # Create trait covariance correlation plots
      if(fig == 6)
        trait_correlation_figure(model, path)

    }
  }
}

#' Highest density interval figure
#'
#' Summarises posterior intervals in a whisker plot. Only defined for logistic
#' regression of the proportion of exotic species (m0).
#'
#' @usage hdi_interval_figure(model = "m0")
#' @export

hdi_interval_figure <- function(model = "m0", path){

  # Load model
  model_output <- load_model(model, path)

  # Get summary
  pars <- extract_pars(model_output$stan_output,
                       names = c("diff_int", "diff_slope"),
                       index = c("year", "treatment"))

  # Label years, treatments
  pars <- mutate(pars,
    year = factor(year, labels = c(2010:2013, 2015, 2016)),
    treatment = factor(treatment, labels = c("Control", "Slash")),
    parameter = factor(parameter, labels = c("Intercept", "Slope")))

  # Set y-scale
  scale <- ceiling(max(pars$conf_high))

  p <- ggplot(pars, aes(x = year, y = mean)) +
    geom_point(size = 2) +
    geom_errorbar(
      aes(
        ymin = conf_low,
        ymax = conf_high),
      width = 0,
      size = 1,
      alpha = 0.3) +
    geom_hline(
      aes(yintercept = 0),
      linetype = "dashed") +
    facet_grid(parameter ~ treatment) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1)) +
    coord_cartesian(ylim = c(-scale, scale)) +
    labs(
      x = "Year",
      y = "Difference in posterior means (fence - open)",
      title = "Effect of fencing on exotic species dominance") +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = Inf, yend = -Inf)

  # Save plot
  filename <- paste0("figures/", model, "_hdi_interval.png")
  ggsave(filename = filename, plot = p, dpi = 600)

  # Display plot
  print(p)
}

#' Linear tobit predictions figure
#'
#' Predicts the latent values of species abundance along a fertility gradient.
#'
#' @usage linear_tobit_figure(model = "m1")
#' @export

linear_tobit_figure <- function(model, path) {

  # Load model
  model_output <- load_model(model, path)

  # Extract species coefficients
  if(grepl("m3|m5", model)) {
   species_coef <- extract_pars(model_output$stan_output,
      pars = c("B"),
      index = c("species", "coef")) %>%
     filter(parameter == "B") %>%
      mutate(coef = factor(coef, labels = c("int", "fert", "rain")),
           fence = "All",
           treatment = "All")
  }

  if(grepl("m4|m6", model)) {
    # Get environmental covariates
    env_covariates <- model_output$data_list$env_covariates %>%
      select(treatment_id, fence, treatment) %>%
      distinct()

    species_coef <- extract_pars(model_output$stan_output,
                                 pars = c("B"),
                                 index = c("treatment_id", "species", "coef")) %>%
      filter(parameter == "B") %>%
      mutate(coef = factor(coef, labels = c("int", "fert", "rain"))) %>%
      left_join(., env_covariates)
  }

  # Predict values
  pred_cover <- pred_cover(species_coef, range = 2.5)

  scale <- min(c(ceiling(max(pred_cover$mean)), 100))

  p <- ggplot(pred_cover, aes(x = x, group = species)) +
    geom_line(
      aes(y = mean),
      size = 1,
      alpha = 0.3) +
    geom_hline(
      aes(yintercept = 0),
      colour = "red") +
    facet_grid(treatment ~ fence) +
    coord_cartesian(ylim = c(-scale, scale)) +
    labs(
      x = "Total nitrogen (standardised)",
      y = "Latent performance (uncensored)") +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

  # Save plot
  filename <- paste0("figures/", model, "_linear_tobit_predictions.png")
  suppressMessages(ggsave(filename = filename, plot = p, dpi = 600))

  # Display plot
  print(p)
}


#' Gather and predict cover
#'
#' Long dataframe of coefficents from extract_pars(), return predictions.
#'
#' @param estimates   estimates of species coefficients
#' @param range       range to predict over, default is (-2, 2)
#' @param coef        two of c("int", "fert", "rain")
#'
#' @importFrom tidyr spread
#'
#' @export

pred_cover <- function(estimates,
                       range = 2,
                       coef = c("int", "fert")){

  estimates <- filter(estimates, coef %in% coef) %>%
    select(-conf_high, -conf_low) %>%
    spread(coef, mean)

  y <- expand_grid(x = seq(-range, range, length.out = 1000), estimates) %>%
    mutate(mean = int + fert * x)

  return(y)
}



#' Covariance heatmap figure
#'
#' Creates an S x S heatmap of covariances between species. Red is negative,
#' blue is positive.
#'
#' @usage covariance_heatmap_figure(model = "m3")
#'
#' @importFrom forcats fct_reorder
#' @export

covariance_heatmap_figure <- function(model, path) {

  # Load model
  model_output <- load_model(model, path)

  species_list <- model_output$data_list$species_list
  n = max(species_list$species_id)

  env_covariates <- model_output$data_list$env_covariates %>%
    select(treatment_id, fence, treatment) %>%
    distinct()

  # Models have different number of treatments
  if(grepl("m3|m4", model)) {
    Sigma <- extract_pars(model_output$stan_output,
                          c("Sigma"),
                          index = c("species_a", "species_b")) %>%
      filter(parameter == "Sigma") %>%
      left_join(., species_list, by = c("species_a" = "species_id")) %>%
      left_join(., species_list, by = c("species_b" = "species_id")) %>%
      mutate(treatment_name = "Average across all treatments")

    # Order species by magnitude of covariance
    Sigma_ordered <- mutate(Sigma,
        species.x = fct_reorder(species.x, mean, .desc = F),
        species.y = fct_reorder(species.y, mean, .desc = F),
        species_a = as.numeric(species.x),
        species_b = as.numeric(species.y),
        mean = ifelse(species_a == species_b, 0, mean)) %>%
      filter(species_b <= species_a)


  } else if(grepl("m5|m6", model)) {
    Sigma <- extract_pars(model_output$stan_output,
                          c("Sigma"),
                          index = c("treatment_id", "species_a", "species_b")) %>%
      filter(parameter == "Sigma") %>%
      left_join(., species_list, by = c("species_a" = "species_id")) %>%
      left_join(., species_list, by = c("species_b" = "species_id")) %>%
      left_join(.,  env_covariates, by = c("treatment_id" = "treatment_id")) %>%
      mutate(treatment_name = paste(fence, treatment))

    # Order species by magnitude of covariance
    Sigma_ordered <- mutate(Sigma,
        species.x = fct_reorder(species.x, mean, .desc = F),
        species.y = fct_reorder(species.y, mean, .desc = F),
        species_a = as.numeric(species.x),
        species_b = as.numeric(species.y),
        mean = ifelse(species_a == species_b, 0, mean)) %>%
      filter(
        ifelse(grepl("Control", treatment_name),
               as.numeric(as.factor(species_b)) < as.numeric(as.factor(species_a)),
               as.numeric(as.factor(species_b)) > as.numeric(as.factor(species_a))))

  } else {
    printf("Heatmaps only defined for models with interactions (m3-m6)")
    return(NULL)
  }


  p <- ggplot(Sigma_ordered,
              aes(
                x = species.x,
                y = species.y,
                fill = mean)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "blue") +
    labs(x = "",
         y = "",
         fill = "",
         title = "Covariance between species") +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0,
        hjust = 1),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.key.height = unit(1.6, "cm"))

  if(grepl("m5|m6", model)) {
    p <- p +
      facet_wrap(~ fence) +
      annotate("segment", x = 0.5, xend = n + 0.5,
              y = 0.5, yend = n + 0.5, linetype = "dashed") +
      annotate("text", x = 5, y = n - 2, label = "Removal") +
      annotate("text", x = n - 4, y = 3.5, label = "Control") +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)
  }

  # Save plot
  filename <- paste0("figures/", model, "_covariance_heatmap.png")
  ggsave(filename = filename, plot = p,
         width = 12, height = 9, dpi = 600)

  # Display plot
  print(p)
}

#' Interaction network figure
#'
#' Filters between species covariance matrices for significant (non-zero)
#' negative interactions and displays them as a network. Links between species
#' are weighted by the magnitude of the interaction, but are only relative
#' within a given network
#'
#' @usage interaction_network_figure(model = "m6")
#' @import circlize
#' @export

interaction_network_figure <- function(model, path) {

  # Load model
  model_output <- load_model(model, path)

  species_list <- model_output$data_list$species_list

  env_covariates <- model_output$data_list$env_covariates %>%
    select(treatment_id, fence, treatment) %>%
    distinct()

  # Models have different number of treatments
  if(grepl("m3|m4", model)) {
    E = 1

    Sigma <- extract_pars(model_output$stan_output,
                          c("Sigma"),
                          index = c("species_a", "species_b")) %>%
      mutate(treatment_id = 1) %>%
      filter(parameter == "Sigma") %>%
      left_join(., species_list, by = c("species_a" = "species_id")) %>%
      left_join(., species_list, by = c("species_b" = "species_id"))
  }
  else if (grepl("m5|m6", model)) {
    E = max(env_covariates$treatment_id)

    Sigma <- extract_pars(model_output$stan_output,
                        c("Sigma"),
                        index = c("treatment_id", "species_a", "species_b")) %>%
    filter(parameter == "Sigma") %>%
    left_join(., species_list, by = c("species_a" = "species_id")) %>%
    left_join(., species_list, by = c("species_b" = "species_id")) %>%
    left_join(.,  env_covariates, by = c("treatment_id" = "treatment_id"))
  }
  else{
    printf("Networks only defined for models with interactions (m3-m6)")
    return(NULL)
  }

  # Create a plot for each treatment
  n <- max(Sigma$species_a)
  factors <- 1:n

  # Scale covariances for figure
  Sigma <- mutate(Sigma, interaction = mean / min(mean) * 2)

  for(e in 1:E){

    # Get figure details
    treatment = filter(env_covariates, treatment_id == e)
    name = ifelse(E == 1, "all", paste0(treatment$fence, treatment$treatment))
    filename = paste0("figures/", model, "_network_", name, ".png")

    # Filter out significant negative interactions (upper limit below zero).
    interactions <-
      filter(Sigma,
             conf_high < 0,
             treatment_id == e)

    # Check that there are interactions to plot
    if(nrow(interactions) == 0){
      printf(paste(model, name, "has no significant negative interactions"))
      next
    }

    # Open device
    circos.clear()
    png(filename= filename,
        width = 7200,
        height = 7200,
        res = 600,
        pointsize = 18)

    # Initialise plot
    circos.par(
      track.height = 0.1,
      gap.after = 0,
      canvas.xlim = c(-1.5, 1.5),
      canvas.ylim = c(-1.6,
                      1.6),
      cell.padding = c(0.01, 0.01, 0.01, 0.01)
    )
    circos.initialize(factors = factors, xlim = c(0, 1))


    # Buffer between labels and links
    circos.track(
      ylim = c(0, 1),
      factors = factors,
      track.height = 0.05,
      bg.border = NA
    )


    # Tweak labels here --------------------------------------------------

    # eg. sorting, colouring

    # Colour labels by introduced status
    #labels <- data.frame(species = sort(species_list$species))
    #colour = ifelse(traits$introduced == 1, "black", "grey")

    labels = gsub("\\.", " ", species_list$species)

    # Add labels
    circos.trackText(
      x = rep(0.5, n),
      y = rep(1, n),
      labels = labels,
      cex = 0.8,
      factors = factors,
      #col = labels$colour,
      font = 1,
      adj = c(0, 0.5),
      facing = "clockwise",
      niceFacing = T
    )

    # Add links for interactions
    for (i in 1:nrow(interactions)){
      circos.link(
        sector.index1 = interactions$species_a[i],
        point1 = 0.5,
        sector.index2 = interactions$species_b[i],
        point2 = 0.5,
        col = "red",
        lwd = interactions$interaction[i],
        h.ratio = 0.8)
    }

    # Save file
    dev.off()
  }
}

#' Covariance change between treatments
#'
#' Creates density plots of covariances of selected species to compare change
#' in interaction strength between treatments.
#'
#' @usage covaraiance_density_figure("m6",
#' species = c("Avena.fatua", "Bromus.diandrus", "Acetosella.vulgaris))
#' @export

covariance_density_figure <- function(model, path,
                                      species = c("Avena.fatua",
                                                  "Bromus.diandrus",
                                                  "Acetosella.vulgaris")) {

  # Load model
  model_output <- load_model(model, path)

  env_covariates <- model_output$data_list$env_covariates %>%
    select(treatment_id, fence, treatment) %>%
    distinct()

  # Models have different number of treatments
  if(!grepl("m5|m6", model)) {
    printf("Figure only defined for models with varying covariances")
    break
  }

  E = max(env_covariates$treatment_id)

  Sigma <- extract_pars(model_output$stan_output,
                        c("Sigma"),
                        index = c("treatment_id", "species_a", "species_b")) %>%
    filter(parameter == "Sigma") %>%
    left_join(., species_list, by = c("species_a" = "species_id")) %>%
    left_join(., species_list, by = c("species_b" = "species_id")) %>%
    left_join(.,  env_covariates, by = c("treatment_id" = "treatment_id"))

  sp_cov <- filter(Sigma, species.x %in% species,
                   species_a != species_b,
                   mean <= 0,
                   mean > min(mean)) %>%
    mutate(treatment2 = paste(fence, treatment))

  p <- ggplot(sp_cov, aes(x = mean,
                     group = treatment2,
                     linetype = fence)) +
    geom_density() +
    geom_vline(aes(xintercept = 0), size = 1) +
    facet_wrap(species.x ~ treatment, scales = "free_y", ncol = 2) +
    coord_cartesian(expand = F) +
    labs(x = "Between species covariance",
         y = "Density",
         linetype = "") +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

  # Save plot
  filename <- paste0("figures/", model, "changing_covariance_density.png")
  ggsave(filename = filename, plot = p,
         width = 12, height = 9, dpi = 600)

  # Display plot
  print(p)

}

#' Changes in covariance with relative trait values
#'
#' Plots correlations between the pairwise differences in species traits and
#' between species covariances.
#'
#' The trait dataset includes:
#' \enumerate{
#'  \item{Leaf.length}
#'  \item{Leaf.width}
#'  \item{Max.height}
#'  \item{SLA}
#'  \item{Vegetative.height}
#'  \item{Vegetative.width}
#' }
#'
#' @usage trait_correlation_figure("m4", trait = "SLA")
#' @export

trait_correlation_figure <- function(model, path, trait_id = "Max.height"){

  # Load model
  model_output <- load_model(model, path)

  species_list <- model_output$data_list$species_list

  env_covariates <- model_output$data_list$env_covariates %>%
    select(treatment_id, fence, treatment) %>%
    distinct()

  # Models have different number of treatments
  if(!grepl("m5|m6", model)) {
    printf("Figure only defined for models with varying covariances")
    break
  }

  Sigma <- extract_pars(model_output$stan_output,
                        c("Sigma"),
                        index = c("treatment_id", "species_a", "species_b")) %>%
    filter(parameter == "Sigma") %>%
    left_join(., species_list, by = c("species_a" = "species_id")) %>%
    left_join(., species_list, by = c("species_b" = "species_id")) %>%
    left_join(.,  env_covariates, by = c("treatment_id" = "treatment_id")) %>%
    rename(mean_cov = mean)

  trait <- filter(traits, grepl(paste0(trait_id, collapse = "|"), trait)) %>%
    tidyr::gather(quantile, value, mean:max) %>%
    filter(quantile == "max")

  trait_diff <- filter(Sigma, species.x %in% species_list$species,
           species_a != species_b) %>%
    left_join(trait, by = c("species.x" = "species")) %>%
    left_join(trait, by = c("species.y" = "species")) %>%
    filter(quantile.x == quantile.y) %>%
    mutate(difference = abs(value.x - value.y))

  p <- ggplot(trait_diff, aes(x = difference,
                              y = mean_cov)) +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_point(alpha = 0.8) +
    geom_smooth(method = "lm", se = F, color = "red") +
    facet_grid(fence ~ treatment) +
    coord_cartesian(expand = F, ylim = c(-300, 300)) +
    labs(x = paste("Absolute", trait_id, "difference"),
         y = "Between species covariance") +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

  # Save plot
  filename <- paste0("figures/", model, "_covariance_trait_correlation.png")
  ggsave(filename = filename, plot = p,
         width = 12, height = 9, dpi = 600)

  # Display plot
  print(p)

}


#' Extract parameter estimates
#'
#' Stan output provides summaries of parameter estimates, this can be
#' faster than working directly with samples
#'
#' @param fit output of a Stan model
#' @param pars vector of parameter names
#' @param index vector of labels for parameter indexes
#'
#' @usage pars <- extract_par(model_output, c("B"),
#' c("species", "covariate", "treatment"))


extract_pars <- function(fit, pars, index) {

  # Get model summary
  summary <- rstan::summary(fit)$summary %>%
    as.data.frame()

  # Grep for desired parameter estimates
  suppressWarnings(estimates <- summary %>%
    mutate(id = rownames(.)) %>%
    filter(grepl(paste(pars, collapse = "|"), id)
                  & !grepl("raw", id)) %>%
    tidyr::separate(id,
                    into = c("parameter", index),
                    sep = "\\[|,|\\]", extra = "drop") %>%
    mutate_at(vars(one_of(index)), as.numeric) %>%
    select(parameter, index, mean,
                  conf_low = `2.5%`, conf_high = `97.5%`))

  return(estimates)
}
