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
#'  (m1-3)}
#'  \item{heatmaps of between species covariances (m2 & m3)}
#'  \item{network diagrams of negative interactions between species (m1-3)}
#'  \item{density plot of changes in negative interactions between treatments (m2 & m3)}
#'  \item{trait-covariance correlation plots, varying by treatment (m3)}
#' }
#'
#' @param model m0-m6, defaults to all.
#' @param path  Directory to save model outputs, defaults to /models.
#' @param figs  numeric: 1-4, defaults to all.
#'
#' @usage create_figures_models(model = "m0", fig = 1)
#' @export

create_figures <-  function(
  models = c("m0", "m1", "m2", "m3"),
  path = "models/",
  figs = 1:6, ...) {

  dir.create("figures/", showWarnings = F)

  for(model in models) {

    # Load model
    model_output <- load_model(model, path)

    for(fig in figs){

      # Create HDI figure for logistic regression
      if(fig == 1 & model == "m0") {
        nonnative_dominance_figures(model_output)
        figs = 1
      }

      # Create linear predictions for tobit models
      if(fig == 2)
        linear_tobit_figure(model, model_output)

      # Create covariance heatmap for joint tobit models
      if(fig == 3 & model == "m1")
        correlation_corrplot_figure(model, model_output)

      # Create interaction network diagram for joint tobit models
      if(fig == 4)
        interaction_network_figure(model, model_output)

      # Create interaction density figure for models with varying covariances
      if(fig == 5 & model == "m3")
        covariance_distribution_figure(model, model_output)

      # Create trait covariance correlation plots
      if(fig == 6 & model == "m3")
        trait_correlation_figure(model)

      # if(fig == 7)
      #   impact_heatmap_figure(model, model_output)

    }
  }
}

#' Nonnative dominance
#'
#' Plots mean regression lines for each treatment and annual change in
#' nonnative-fertility relationship following fencing.
#'
#' @usage dominance_figures(model = "m0")
#'
#' @importFrom tidyr gather spread unite
#'
#' @export

nonnative_dominance_figures <- function(model_output){

  # Set up covariates
  years <- data.frame(y = 1:6,
                      year = c(2010, 2011, 2012, 2013, 2015, 2016))

  treatments <- data.frame(t = 1:4,
                           treatment = c("Unslashed Fenced", "Slashed Fenced",
                                         "Unslashed Grazed", "Slashed Grazed"),
                           slash = c("Unslashed", "Slashed",
                                     "Unslashed", "Slashed"),
                           fence = c("Fenced", "Fenced", "Grazed", "Grazed"))

  # Extract raw data
  raw_data <- data.frame(x = model_output$data_list$X[, 2],
                         x2 = model_output$data_list$X[, 3],
                   y_raw = model_output$data_list$y_observed,
                   y_scaled = unlogit(model_output$data_list$y_observed),
                   t = model_output$data_list$treatment,
                   y = model_output$data_list$year) %>%
    left_join(., treatments) %>%
    left_join(., years)


  # Extract mean regression coefficients
  mean_slope <- extract_pars(model_output$model_summary,
                             pars = c("mu_int", "mu_slope"),
                             index = "t")

  # Generate predicted regression lines
  pred <- select(mean_slope, -conf_low, -conf_high) %>%
    spread(parameter, mean) %>%
    expand_grid(x = seq(-2.2, 2.2, len = 1000), .) %>%
    mutate(y_raw = mu_int + mu_slope * x,
           y_scale = unlogit(y_raw)) %>%
    left_join(., treatments)

  # Plot against data
  p1 <- ggplot(data = raw_data,
               aes(x = x,
                   y = y_raw,
                   group = treatment)) +
    geom_point(aes(shape = treatment,
                   colour = x),
               fill = "white",
               size = 3.5) +
    geom_path(data = pred,
              aes(x = x,
                  y = y_raw,
                  linetype = treatment),
              size = 2) +
    coord_cartesian(ylim = c(-4, 6), expand = F) +
    viridis::scale_colour_viridis() +
    scale_shape_manual(values = c(19, 17, 21, 25)) +
    scale_linetype_manual(values = c(1, 5, 4, 3)) +
    labs(x = "Fertility (scaled)",
         y = "Proportion of nonnative species (logit scale)",
         color = "",
         shape = "",
         linetype = "") +
    guides(shape = guide_legend(nrow = 4,
                                keywidth = 5),
           color = F) +
    theme(legend.position = c(1, 0),
          legend.justification = c(1, 0))

  # Extract annual regression coefficients
  slope <- extract_pars(model_output$model_summary,
                        pars = c("B_int", "B_slope"),
                        index = c("t", "y")) %>%
    left_join(., years) %>%
    left_join(., treatments) %>%
    gather(quantile, value, mean:conf_high) %>%
    unite(parameter, parameter, quantile, remove = T) %>%
    spread(parameter, value)

  # Plot change in intercepts overtime, coloured by slope
  p2 <- ggplot(slope,
               aes(x = year,
                   y = B_int_mean,
                   group = treatment)) +
    geom_line(aes(color = B_slope_mean),
              size = 1.8) +
    geom_line(aes(y = B_int_conf_high,
                  linetype = fence))+
    geom_line(aes(y = B_int_conf_low,
                  linetype = fence)) +
    geom_point(aes(shape = treatment,
                   color = B_slope_mean),
               fill = "white",
               size = 4.5) +
    geom_ribbon(aes(ymin = B_int_conf_low,
                    ymax = B_int_conf_high,
                    color = NULL),
                alpha = 0.05) +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    coord_cartesian(xlim = c(2010, 2016.2), ylim = c(-2.4, 5), expand = F) +
    viridis::scale_colour_viridis() +
    scale_shape_manual(values = c(19, 17, 21, 25)) +
    facet_wrap(~ slash, ncol = 1) +
    labs(x = "Year",
         y = "Logitstic intercept",
         color = "Logistic slope",
         linetype = "") +
    guides(shape = F,
           color = guide_colourbar(title.position = "top",
                                   title.hjust = 0.5,
                                   direction = "horizontal",
                                   barwidth = 8,
                                   barheight = 0.8),
           linetype = guide_legend(nrow = 2)) +
    theme(legend.box = "horizontal",
          legend.position = c(1, 0),
          legend.justification = c(1, 0),
          legend.title = element_text(size = 14))

  # Save and print
  p <- grid.arrange(label(p1, "A)"), label(p2, "B)"),
               nrow = 1, widths = c(1.5, 1))


  ggsave(p, filename = "figures/S3_logistic_regression.png",
         height = 8, width = 14, device = "png", dpi = 600)


  # Plot effect of rainfall
  B <- extract_pars(model_output$model_summary,
                    pars = c("mu_int", "B_rain"),
                    index = c("year")) %>%
    group_by(parameter) %>%
    summarise(mean = mean(mean),
              conf_low = min(conf_low),
              conf_high = max(conf_high)) %>%
    gather(quantile, value, mean:conf_high) %>%
    spread(parameter, value) %>%
    mutate

  p4 <- ggplot(raw_data, aes(x = x2, y_raw)) +
    geom_jitter(aes(color = x), size = 3, width = 0.02) +
    geom_abline(data = B, aes(intercept = mu_int, slope = B_rain,
                              linetype = grepl("conf", quantile)),
                size = 1.4) +
    coord_cartesian(xlim = c(-2, 2)) +
    scale_color_viridis_c(limits = c(-2.2, 2.2)) +
    labs(x = "Spring rainfall (scaled)",
         y = "Proportion of nonnative species (logit scale)",
          color = "Fertility (scaled)") +
    guides(linetype = F,
           colour = guide_colourbar(title.position = "top",
                                    barwidth = 8,
                                    barheight = 0.8)) +
    theme(legend.position = c(0.8, 0.1),
          legend.direction = "horizontal")

  ggsave(label(p4, "A"), filename = "figures/S5a_rainfall_figure.png",
         device = "png", dpi = 600, height = 6, width = 6)

  pred_year <-
    extract_pars(model_output$model_summary,
      pars = c("B_int", "B_slope"),
      index = c("t", "y")) %>%
    left_join(., years) %>%
    left_join(., treatments) %>%
    select(-conf_low, -conf_high) %>%
    spread(parameter, mean) %>%
    expand_grid(x = seq(-2.2, 2.3, len = 1000), .) %>%
    mutate(y_raw = B_int + B_slope * x,
      y_scale = unlogit(y_raw))

  # Alternative plot with slopes over time
  p3 <- ggplot(raw_data, aes(x = x, y = y_raw)) +
    geom_line(data = pred_year,
                aes(group = t,
                    linetype = fence),
                size = 1.4) +
    geom_point(aes(fill = fence),
               shape = 21,
               size = 4) +
    facet_grid(slash ~ year) +
    scale_fill_manual(values = c("black", "white")) +
    coord_cartesian(xlim = c(-2.2, 2.3), ylim = c(-4.1, 6.9), expand = F) +
    labs(x = "Fertility (scaled)",
         y = "Proportion of nonnative species (logit scale)",
         fill = "",
         linetype = "") +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    theme(legend.position = c(1, 0.51),
          legend.justification = c(1, 0),
          legend.key.width = unit(1.2, "cm"))


  # Save and print
  ggsave(p3, filename = "figures/F1_logistic_regression_by_year.png",
         height = 5, width = 14, device = "png", dpi = 600)
  print(p3)



}


#' Linear tobit predictions figure
#'
#' Predicts the latent values of species abundance along a fertility gradient.
#'
#' @usage linear_tobit_figure(model = "m1")
#'
#' @export

linear_tobit_figure <- function(model, model_output) {

  # Get species origin
  species <- model_output$data_list$species_list %>%
    mutate(origin = ifelse(introduced == 1, "Nonnative", "Native"))

  # Extract species coefficients
  if(grepl("m2", model)) {
   species_coef <- extract_pars(model_output$model_summary,
      pars = c("B"),
      index = c("species", "coef")) %>%
     filter(parameter == "B") %>%
      mutate(coef = factor(coef, labels = c("int", "fert", "rain")),
           fence = "All",
           treatment = "All")
  } else if(grepl("m1|m3", model)) {
    # Get environmental covariates
    env_covariates <- model_output$data_list$env_covariates %>%
      select(treatment_id, fence, treatment) %>%
      distinct()

    species_coef <- extract_pars(model_output$model_summary,
                                 pars = c("B"),
                                 index = c("treatment_id", "species", "coef")) %>%
      filter(parameter == "B") %>%
      mutate(coef = factor(coef, labels = c("int", "fert", "rain"))) %>%
      left_join(., env_covariates)
  } else {
    printf("Linear predictions only defined for tobit models (m1-3)")
    return(NULL)
  }


  # Predict values, highlight species with positive slopes
  pred_fert <- pred_cover(species_coef, range = 2.5) %>%
    left_join(., species, by = c("species" = "species_id"))

  scale <- min(c(ceiling(max(pred_fert$mean)), 1))

  p <- ggplot(pred_fert,
              aes(x = x,
                  group = desc(species))) +
    geom_line(aes(y = mean,
                  color = as.factor(fert),
                  linetype = origin),
              size = 1.5) +
    geom_hline(aes(yintercept = 0),
               colour = "black") +
    facet_grid(treatment ~ fence) +
    coord_cartesian(xlim = c(-2, 2),
                    ylim = c(-100, 100),
                    expand = F) +
    scale_linetype_manual(values = c("twodash", "solid")) +
    viridis::scale_colour_viridis(begin = 0.2, discrete = T) +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    labs(x = "Fertility (scaled)",
      y = "Latent habitat suitability",
      linetype = "",
      color = "Tobit slope") +
    guides(color = F,
           linetype = guide_legend(nrow = 2, keywidth = 4)) +
    theme(legend.box = "horizontal",
          legend.position = c(1, 0),
          legend.justification = c(1, 0),
          legend.title = element_text(size = 14))

  # Save plot
  filename <- paste0("figures/F2_", model, "_linear_tobit_predictions.png")
  ggsave(filename = filename, plot = p,
         height = 7, width = 7, device = "png", dpi = 600)

  # Display plot
  print(p)

  # Get rain effect quantiles from full posterior
  b_rain <- as.data.frame(model_output$stan_output, pars = "B") %>%
    gather(par, val) %>%
    separate(par, c("par", "treatment_id", "species_id", "coef"), sep = "\\[|,") %>%
    mutate(species_id = as.numeric(species_id),
           treatment_id = as.numeric(treatment_id),
           coef = factor(coef, labels = c("int", "fert", "rain"))) %>%
    group_by(species_id, treatment_id, coef) %>%
    summarise(mean = mean(val),
              conf_low = quantile(val, 0.1),
              conf_high = quantile(val, 0.90)) %>%
    left_join(species) %>%
    left_join(env_covariates) %>%
    filter(coef == "rain")

  p2 <- ggplot(b_rain, aes(x = reorder(species, mean), shape = origin)) +
    geom_hline(aes(yintercept = 0)) +
    geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0, size = 1) +
    geom_point(aes(y = mean), size  = 2, fill = "white") +
    facet_grid(treatment ~ fence) +
    labs(x = "", y = "Rainfall effect", shape = "") +
    scale_shape_manual(values = c(21, 16)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
          legend.position = c(0.9, 0.06))

  filename <- paste0("figures/S5b_", model, "_rainfall_effect.png")
  ggsave(filename = filename, plot = label(p2, "B"),
         height = 12, width = 10, device = "png", dpi = 600)

  pred_rain <- filter(species_coef, grepl("int|rain", coef)) %>%
    select(-conf_high, -conf_low) %>%
    spread(coef, mean) %>%
    expand_grid(x = seq(-2.2, 2, length.out = 1000), .) %>%
    mutate(mean = int + rain * x) %>%
    left_join(., species, by = c("species" = "species_id"))

  p3 <- ggplot(pred_rain,
              aes(x = x,
                  group = desc(species))) +
    geom_line(aes(y = mean,
                  color = as.factor(rain),
                  linetype = origin),
              size = 1.5) +
    geom_hline(aes(yintercept = 0),
               colour = "black") +
    facet_grid(treatment ~ fence) +
    coord_cartesian(xlim = c(-2, 2),
                    ylim = c(-100, 100),
                    expand = F) +
    scale_linetype_manual(values = c("twodash", "solid")) +
    viridis::scale_colour_viridis(begin = 0.2, discrete = T) +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    labs(x = "Spring rainfall (scaled)",
         y = "Latent habitat suitability",
         linetype = "",
         color = "Tobit slope") +
    guides(color = F,
           linetype = guide_legend(nrow = 2, keywidth = 4)) +
    theme(legend.box = "horizontal",
          legend.position = c(1, 0),
          legend.justification = c(1, 0),
          legend.title = element_text(size = 14))

  filename <- paste0("figures/S5b_", model, "_rainfall_tobit_prediction.png")
  ggsave(filename = filename, plot = label(p3, "B"),
         height = 7, width = 7, device = "png", dpi = 600)

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

covariance_heatmap_figure <- function(model, model_output) {

  species_list <- model_output$data_list$species_list
  n = max(species_list$species_id)

  env_covariates <- model_output$data_list$env_covariates %>%
    select(treatment_id, fence, treatment) %>%
    distinct()

  # Models have different number of treatments
  if(grepl("m1", model)) {
    Sigma <- extract_pars(model_output$model_summary,
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


  } else if(grepl("m2|m3", model)) {
    Sigma <- extract_pars(model_output$model_summary,
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
    printf("Heatmaps only defined for models with interactions (m1-3)")
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

  if(grepl("m2|m3", model)) {
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
         width = 10, height = 9, dpi = 600)

  # Display plot
  print(p)
}


#' Correlation heatmap figure
#'
#' Creates an S x S heatmap of covariances between species. Red is negative,
#' blue is positive.
#'
#' @usage correlation_heatmap_figure(model = "m3")
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom forcats fct_reorder
#' @export

correlation_heatmap_figure <- function(model, model_output) {

  species_list <- model_output$data_list$species_list %>%
    mutate(species = ifelse(introduced == 1,
                            paste0(gsub("\\.", " ", species), "*"),
                            paste0(gsub("\\.", " ", species), " ")))

  n = max(species_list$species_id)

  env_covariates <- model_output$data_list$env_covariates %>%
    select(treatment_id, fence, treatment) %>%
    distinct()

  # Models have different number of treatments
  if(grepl("m1", model)) {
    sigma <- extract_pars(model_output$model_summary,
                          c("sigma"),
                          index = "species") %>%
      filter(parameter == "sigma") %>%
      mutate(variance = mean^2) %>%
      select(species, variance)

    Omega <- extract_pars(model_output$model_summary,
                          c("Omega"),
                          index = c("species_a", "species_b")) %>%
      filter(parameter == "Omega") %>%
      left_join(., sigma, by = c("species_a" = "species")) %>%
      left_join(., sigma, by = c("species_b" = "species")) %>%
      left_join(., species_list, by = c("species_a" = "species_id")) %>%
      left_join(., species_list, by = c("species_b" = "species_id"))


    # Order species by magnitude of Correlation
    Omega_ordered <- mutate(Omega,
                            species.x = fct_reorder(species.x, variance.x, .desc = T),
                            species.y = fct_reorder(species.y, variance.y, .desc = T),
                            species_a = as.numeric(species.x),
                            species_b = as.numeric(species.y),
                            mean = ifelse(species_a == species_b, 0, mean)) %>%
      filter(species_b <= species_a)

  } else if(grepl("m2|m3", model)) {

    sigma <- extract_pars(model_output$model_summary,
                          c("sigma"),
                          index = c("treatment_id", "species")) %>%
      filter(parameter == "sigma")


    Omega <- extract_pars(model_output$model_summary,
                          c("Omega"),
                          index = c("treatment_id", "species_a", "species_b")) %>%
      filter(parameter == "Omega") %>%
      left_join(., species_list, by = c("species_a" = "species_id")) %>%
      left_join(., species_list, by = c("species_b" = "species_id")) %>%
      left_join(.,  env_covariates, by = c("treatment_id" = "treatment_id")) %>%
      mutate(treatment_name = paste(fence, treatment)) %>%
      group_by(species_a) %>%
      mutate(total_correlation_a = sum(mean)) %>%
      group_by(species_b) %>%
      mutate(total_correlation_b = sum(mean)) %>%
      ungroup()

    # Order species by magnitude of Correlation
    Omega_ordered <- mutate(Omega,
                            species.x = fct_reorder(species.x, total_correlation_a, .desc = F),
                            species.y = fct_reorder(species.y, total_correlation_b, .desc = F),
                            species_a = as.numeric(species.x),
                            species_b = as.numeric(species.y),
                            mean = ifelse(species_a == species_b, 0, mean)) %>%
      filter(
        ifelse(grepl("Control", treatment_name),
               as.numeric(as.factor(species_b)) <= as.numeric(as.factor(species_a)),
               as.numeric(as.factor(species_b)) > as.numeric(as.factor(species_a))))

  } else {
    printf("Heatmaps only defined for models with interactions (m1-3)")
    return(NULL)
  }

  p <- ggplot(Omega_ordered,
              aes(x = species.x,
                y = species.y,
                fill = mean)) +
    geom_tile() +
    scale_fill_gradient2(low = "red",
                         mid = "white",
                         high = "blue",
      limits = c(-0.51, 0.51)) +
    labs(x = "",
         y = "",
         fill = "Between species correlation") +
    guides(fill = guide_colorbar(title.position = "top",
                                 title.hjust = 0.5,
                                 direction = "horizontal",
                                 barwidth = 16)) +
    theme(axis.text = element_text(face = "italic"),
      axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
      axis.ticks = element_blank(),
      legend.position = c(0, 1),
      legend.justification = c(-0.1, 1),
      plot.margin = unit(c(0, 0, 0, 0), "in"))

  if(grepl("m2|m3", model)) {
    p <- p +
      facet_wrap(~ fence) +
      annotate("segment", x = 0.5, xend = n + 0.5,
               y = 0.5, yend = n + 0.5, linetype = "dashed") +
      annotate("text", x = 5, y = n - 2, label = "Removal") +
      annotate("text", x = n - 4, y = 3.5, label = "Control") +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)
  }

  # # Save plot
  # filename <- paste0("figures/", model, "_correlation_heatmap.png")
  # ggsave(filename = filename, plot = p,
  #        width = 10, height = 10, dpi = 600)
  #
  # # Display plot
  # print(p)

  # Plot variance only
  variance <- filter(Omega_ordered, species_a == species_b) %>%
    mutate(species.x = "Variance")

  # Or with barplot
  v2 <- ggplot(sigma, aes(x = species_b, y = variance.x)) +
    geom_bar(stat = "identity") +
    coord_flip(expand = F) +
    labs(x = "",
         y = "Variance") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          aspect.ratio = 1.25,
          plot.margin = unit(c(0.4, 0.3, 2, 0), "in"))

  if(grepl("m2|m3", model)) {
    v2 <- v2 +
      facet_wrap(~ fence) +
      annotate("segment", x = 0.5, xend = n + 0.5,
               y = 0.5, yend = n + 0.5, linetype = "dashed") +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)
  }

  pv <- grid.arrange(label(p, "B)", 0),
                     label(v2, "C)", 0),
                     nrow = 1,
                     widths = c(1.5, 1))

  filename <- paste0("figures/F3b_", model, "_covariance_barplot.png")
  ggsave(filename = filename, plot = pv,
         width = 14, height = 9, dpi = 600)

}

#' Corrplot
#'
#' @export

correlation_corrplot_figure <- function(model, model_output) {

  species_list <- model_output$data_list$species_list %>%
    mutate(species = ifelse(introduced == 1,
                            paste0(gsub("\\.", " ", species), "*"),
                            paste0(gsub("\\.", " ", species), " ")))

  Omega <- extract_pars(model_output$model_summary, "Omega", c("sp_a", "sp_b")) %>%
    left_join(., species_list, by = c("sp_a" = "species_id")) %>%
    left_join(., species_list, by = c("sp_b" = "species_id")) %>%
    mutate(mean = ifelse(sp_a == sp_b, 0, mean)) %>%
    select(species.x, species.y, mean) %>%
    spread(species.y, mean) %>%
    select(-species.x)

  rownames(Omega) <- colnames(Omega)

  filename <- paste0("figures/", model, "_correlation_corrplot.png")
  png(filename = filename,
      width = 9,
      height = 9,
      units = "in",
      res = 600)

  m <- as.matrix(Omega)
  o <- rev(corrplot::corrMatOrder(m, "FPC"))

  p <- corrplot::corrplot(m[o, o],
           is.corr = F,
           diag = F,
           type = "lower",
           method = "circle",
           tl.srt = 50,
           tl.offset = 1,
           tl.col = "black",
           tl.cex = 0.8,
           font = 3,
           cl.lim = c(-0.5, 0.5))

  dev.off()

  sigma <- extract_pars(model_output$model_summary,
                        c("sigma"),
                        index = "species_id") %>%
    filter(parameter == "sigma") %>%
    mutate(variance = mean^2) %>%
    select(species_id, variance) %>%
    left_join(., species_list)

  sigma_ordered <- sigma[o,] %>%
    mutate(species_id = n():1)

  v2 <- ggplot(sigma_ordered, aes(x = species_id, y = variance)) +
    geom_bar(stat = "identity") +
    coord_flip(expand = F) +
    scale_y_reverse() +
    labs(x = "",
         y = "Variance") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          aspect.ratio = 1.25)
         # plot.margin = unit(c(0.4, 0.3, 2, 0), "in"))

  filename <- paste0("figures/F3c_", model, "_covariance_barplot.png")
  ggsave(filename = filename, plot = v2,
    width = 14, height = 9, dpi = 600)
}

#' Interaction network figure
#'
#' Filters between species covariance matrices for significant (non-zero)
#' negative interactions and displays them as a network. Links between species
#' are weighted by the magnitude of the interaction, but are only relative
#' within a given network
#'
#' @param subset select positive, negative or specific interactions
#'
#' @usage interaction_network_figure(model = "m3")
#' @export

interaction_network_figure <- function(model, model_output, subset = c("negative")) {

  species_list <- model_output$data_list$species_list

  env_covariates <- model_output$data_list$env_covariates %>%
    select(treatment_id, fence, treatment) %>%
    distinct()

  # Models have different number of treatments
  if(grepl("m1", model)) {
    E = 1

    Sigma <- extract_pars(model_output$model_summary,
                          c("Sigma"),
                          index = c("species_a", "species_b")) %>%
      mutate(treatment_id = 1) %>%
      filter(parameter == "Sigma",
             species_a != species_b) %>%
      left_join(., species_list, by = c("species_a" = "species_id")) %>%
      left_join(., species_list, by = c("species_b" = "species_id"))
  }
  else if (grepl("m2|m3", model)) {
    E = max(env_covariates$treatment_id)

    Sigma <- extract_pars(model_output$model_summary,
                        c("Sigma"),
                        index = c("treatment_id", "species_a", "species_b")) %>%
    filter(parameter == "Sigma",
           species_a != species_b) %>%
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
  Sigma <- mutate(Sigma, interaction = as.numeric(abs(scale(mean))) / 1.5)


  for(e in 1:E){

    # Get figure details
    treatment = filter(env_covariates, treatment_id == e)
    name = ifelse(E == 1, "all", paste0(treatment$fence, treatment$treatment))
    filename = paste0("figures/F3a_", model, "_", subset, "_network_", name, ".png")

    # Filter out significant negative interactions (upper limit below zero).
    if(subset == "negative"){
      interactions <-
        filter(Sigma,
               conf_high < 0,
               treatment_id == e)

      int_colour = "red"

      printf("Displaying negative interactions")
    } else if(subset == "positive") {
      interactions <-
          filter(Sigma,
                 conf_low > 0,
                 treatment_id == e)

      int_colour = "blue"

      printf("Displaying positive interactions")
    } else {
      printf("Interaction figure for individual species not yet defined")
    }

    # Check that there are interactions to plot
    if(nrow(interactions) == 0){
      printf(paste(model, name, "has no significant interactions"))
      next
    }

    # Open device
    circlize::circos.clear()
    png(filename= filename,
        width = 6000,
        height = 5400,
        res = 600,
        pointsize = 14)

    # Initialise plot
    circlize::circos.par(
      track.height = 0.1,
      gap.after = 0,
      canvas.xlim = c(-1.5, 1.5),
      canvas.ylim = c(-1.6,
                      1.6),
      cell.padding = c(0.01, 0.01, 0.01, 0.01)
    )
    circlize::circos.initialize(factors = factors, xlim = c(0, 1))


    # Buffer between labels and links
    circlize::circos.track(
      ylim = c(0, 1),
      factors = factors,
      track.height = 0.05,
      bg.border = NA
    )


    # Tweak labels here --------------------------------------------------

    # eg. sorting, colouring

    # Colour labels by introduced status
    #labels <- data.frame(species = sort(species_list$species))
    colour = ifelse(species_list$introduced == 1, "black", "grey")

    labels = gsub("\\.", " ", species_list$species)

    # Add labels
    circlize::circos.trackText(
      x = rep(0.5, n),
      y = rep(1, n),
      labels = labels,
      cex = 1,
      factors = factors,
      col = colour,
      font = 3,
      adj = c(0, 0.5),
      facing = "clockwise",
      niceFacing = T
    )

    # Add links for interactions
    for (i in 1:nrow(interactions)){
      circlize::circos.link(
        sector.index1 = interactions$species_a[i],
        point1 = 0.5,
        sector.index2 = interactions$species_b[i],
        point2 = 0.5,
        col = int_colour,
        lwd = interactions$interaction[i],
        h.ratio = 0.8)
    }

    # Save file
    dev.off()
  }
}

#' Combine covariance heatmap and network
#'
#' @importFrom png readPNG
#' @importFrom grid rasterGrob
#' @export

heatmap_network_figure <- function(model, model_output) {

  if(model != "m1"){
    printf("Combined figure only defined for m1")
  }

  # Create figures
  correlation_heatmap_figure(model, model_output)
  interaction_network_figure(model, model_output)

  # Load as images
  p1 <- readPNG(paste0("figures/", model, "_covariance_barplot.png"))
  p2 <- readPNG(paste0("figures/", model, "_negative_network_all.png"))

  p <- arrangeGrob(rasterGrob(p1),
               label(rasterGrob(p2), "C)"),
               nrow = 2,
               widths = c(2, 1))

  filename <- paste0("figures/", model, "_covariance_heatmap_network.png")
  ggsave(filename = filename, plot = p,
         width = 6, height = 9, dpi = 600)

}

#' Covariance change between treatments
#'
#' Creates density plots of covariances of selected species to compare change
#' in interaction strength between treatments.
#'
#' @usage covaraiance_density_figure("m6", species = c("Avena.fatua", "Bromus.diandrus", "Acetosella.vulgaris"))
#' @export


covariance_distribution_figure <- function(model, model_output,
                                      species = c("Avena.fatua",
                                                  "Bromus.diandrus",
                                                  "Acetosella.vulgaris")) {

  species_list <- model_output$data_list$species_list

  env_covariates <- model_output$data_list$env_covariates %>%
    select(treatment_id, fence, treatment) %>%
    distinct()

  # Models have different number of treatments
  if(!grepl("m2|m3", model)) {
    printf("Covariance distribution figure only defined for models with varying covariances")
    return(NULL)
  }

  E = max(env_covariates$treatment_id)

  Sigma <- extract_pars(model_output$model_summary,
                        c("Sigma"),
                        index = c("treatment_id", "species_a", "species_b")) %>%
    filter(parameter == "Sigma") %>%
    left_join(., species_list, by = c("species_a" = "species_id")) %>%
    left_join(., species_list, by = c("species_b" = "species_id")) %>%
    left_join(.,  env_covariates, by = c("treatment_id" = "treatment_id"))

  sp_cov <- Sigma %>%
    filter(species.x %in% species,
           species_a != species_b) %>%
    group_by(species.x, fence, treatment) %>%
    summarise(ucl = quantile(mean, probs = 0.975),
              ucr = quantile(mean, probs = 0.75),
              med = quantile(mean, probs = 0.50),
              lcr = quantile(mean, probs = 0.25),
              lcl = quantile(mean, probs = 0.025)) %>%
    ungroup() %>%
    mutate(species.x = paste0(gsub("\\.", " ", species.x), "*"))

  s <- ggplot(sp_cov, aes(x = treatment,
                          fill = fence,
                          shape = fence)) +
    geom_hline(aes(yintercept = 0)) +
    geom_errorbar(aes(ymin = lcl, ymax = ucl),
                  position = position_dodge(width = 0.3),
                  width = 0,
                  size = 1.5) +
    geom_errorbar(aes(ymin = lcr, ymax = ucr),
                  position = position_dodge(width = 0.3),
                  width = 0,
                  size = 3.2,
                  color = "dodgerblue3") +
    geom_point(aes(y = med),
               position = position_dodge(width = 0.3),
               size = 4.5) +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    scale_fill_manual(values = c("black", "white")) +
    scale_shape_manual(values = c(21, 21)) +
    facet_grid(~ species.x) +
    labs(y = "Covariance",
         x = "",
         fill = "",
         colour = "",
         shape = "") +
    guides(linetype = guide_legend(nrow = 2,
                                   keywidth = 3)) +
    theme(strip.text = element_text(face = "italic"),
          legend.position = c(1, 0),
          legend.justification = c(1, 0),
          legend.box = "horizontal")

  filename <- paste0("figures/F4_", model, "_covariance_distributions.png")
  ggsave(filename = filename, plot = s,
         width = 12, height = 5, dpi = 600)

  # Display plot
  print(s)

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
#'
#' @importFrom tidyr gather
#' @export

trait_correlation_figure <- function(model){

  model_output <- load_model("m3_trait_regression")

  tr <- data.frame(tr = 1:nlevels(factor(traits$trait)),
    trait = levels(factor(traits$trait))) %>%
    mutate(trait = gsub("\\.", " ", trait),
      trait = gsub("Vegetative", "Canopy", trait))

  t <- data.frame(t = 1:model_output$data_list$N_treatment,
    treatment = c("Unslashed", "Slashed", "Unslashed", "Slashed"),
    fence = c("Fenced", "Fenced", "Grazed", "Grazed"))

  mu <- extract_pars(model_output$model_summary, "B_slope", c("tr", "t")) %>%
    left_join(., tr) %>%
    left_join(., t)

  p <- ggplot(mu, aes(x = trait,
    fill = fence,
    shape = fence)) +
    geom_hline(aes(yintercept = 0)) +
    geom_errorbar(aes(ymin = conf_low, ymax = conf_high),
      size = 1.5,
      position = position_dodge(width = 0.5),
      width = 0) +
    geom_point(aes(y = mean),
      size = 4.5,
      position = position_dodge(width = 0.5)) +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    scale_fill_manual(values = c("black", "white")) +
    scale_shape_manual(values = c(21, 21)) +
    facet_grid( ~ treatment) +
    labs(x = "",
      y = "Slope coefficient",
      shape = "",
      fill = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.box = "horizontal")

  # Save plot
  filename <- paste0("figures/F5_", model, "_covariance_trait_correlation.png")
  ggsave(filename = filename, plot = p,
         width = 12, height = 9, dpi = 600)

  # Display plot
  print(p)

}

tobit_example_figure <- function(n = 25){

  tobit <- data.frame(x = seq(-2, 2, len = n)) %>%
    mutate(Uncensored = x * 50 + rnorm(n, 0, 10),
      Censored = ifelse(Uncensored < 0 , 0, Uncensored)) %>%
    gather(obs, val, -x)

  p <- ggplot(tobit, aes(x, val, shape = obs)) +
    geom_abline(aes(intercept = 0, slope = 50),
                linetype = "dashed",
                color = "red",
                size = 1.5) +
    geom_hline(aes(yintercept = 0)) +
    geom_point(size = 4, fill = "white") +
    coord_cartesian(xlim = c(-2, 2.2),
                    ylim = c(-100, 100),
                    expand = F) +
    scale_shape_manual(values = c(16, 21)) +
    scale_y_continuous(sec.axis = dup_axis(name = "Observed cover",
                                            breaks = c(0, 50, 100))) +
    annotate("segment", x = Inf, xend = Inf, y = -0.5, yend = Inf) +
    labs(x = "Environmental gradient",
         y = "Latent habitat suitability",
         shape = "") +
    theme(legend.position = c(1, 0),
          legend.justification = c(1, 0),
          axis.line.y.right = element_blank())

  ggsave(filename = "figures/S1_tobit_example.png", plot = p,
         height = 6, width = 6, device = "png", dpi = 600)

  print(p)
}

species_abundance_distribution_figure <- function(years = 2013:2016) {

  prop <- filter(cover, year %in% years) %>%
    subset_species_plots(., threshold = 0) %>%
    mutate(common = ifelse(prop_plots > .19, "Common", "Uncommon"),
           species = ifelse(introduced == 1,
                            paste0(gsub("\\.", " ", species), "*"),
                            paste0(gsub("\\.", " ", species), " "))) %>%
    filter(!is.na(species))

  p <- ggplot(prop, aes(x = reorder(species, prop_plots),
                   y = prop_plots,
                   fill = common)) +
    geom_bar(stat = "identity") +
    coord_cartesian(expand = F) +
    labs(x = "",
         y = "Proportion of plots occupied",
         fill = "") +
    scale_fill_manual(values = c("black", "grey")) +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size = 10),
          aspect.ratio = 0.2,
          legend.position = c(0, 1),
          legend.justification = c(-0.5, 1))

  ggsave(filename = "figures/S2_species_abundance_distribution.png", plot = p,
         height = 6, width = 18, device = "png", dpi = 400)

}

#' Extract parameter estimates
#'
#' Stan output provides summaries of parameter estimates, this can be faster than working directly with samples
#'
#' @param summary output of a Stan model summary
#' @param pars vector of parameter names
#' @param index vector of labels for parameter indexes
#'
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr separate
#'
#' @examples
#' pars <- extract_par(fit = model_output$model_summary,
#'                           pars = c("B"),
#'                           index = c("species", "covariate", "treatment"))


extract_pars <- function(summary, pars, index) {

  # Grep for desired parameter estimates
  suppressWarnings(estimates <- summary %>%
    rownames_to_column(var = "id") %>%
    filter(grepl(paste0("^", pars, collapse = "|"), id)
                  & !grepl("raw", id)) %>%
    separate(id,
             into = c("parameter", index),
             sep = "\\[|,|\\]", extra = "drop") %>%
    mutate_at(vars(one_of(index)), as.numeric) %>%
    select(parameter, index, mean,
                  conf_low = `2.5%`, conf_high = `97.5%`))

  return(estimates)
}
