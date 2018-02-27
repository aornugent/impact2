#' Check models
#'
#' Having used \code{run_models()}, we can now evaluate the convergence of
#' the model and a variety of posterior predictive checks. These results
#' are output to the console.
#'
#' @param model m0-m6, defaults to all.
#' @param path  Directory to save model outputs, defaults to /models.
#'
#' @usage check_models(model = "m1")
#'
#' @export

check_models <- function(
  models = c("m0", "m1", "m2", "m3", "m4", "m5", "m6"),
  path = "models/",
  model_output = NA, ...) {

  for(model in models) {

    # Load model if need be
    if(is.na(model_output[1]))
      model_output <- load_model(model, path, ...)

    # Check for convergence
    convergence_check(model_output)

    # Run posterior checks, not defined for m0
    if(model != "m0")
      posterior_check(model_output)

    # Delete loaded output to save space
    rm(model_output)

  }
}


#' Load model
#'
#' Loads model output from .Rdata file.
#'
#' @param model m0-m6, defaults to all.
#' @param path  Directory to save model outputs, defaults to /models.
#'
#' @usage model_output <- load_models(model = "m1", path = "models/")
#'
#' @export

load_model <- function(model, path, ...) {

  # Format filename
  filename = paste0(path, model, "_output.Rdata")

  printf(paste("Loading:", filename))

  # Check that this fails nicely
  tryCatch(
    model_output <- get(load(filename)),
    error = function(e) {
      printf(paste(filename, "does not exist"))
    }
  )

  # Check model type
  if(class(model_output$stan_output) != "stanfit"){
    printf(paste(filename, "does not contain the output of a Stan model"))
    return(NULL)
  }

  return(model_output)
}


#' Convergence checks
#'
#' Returns the model run time, the number of parameters with Rubin-Gelman
#' statistics greater than 1.1 (indicating that they have not converged),
#' and the number of effective samples of the worst sampled parameter.
#'
#' @usage convergence_check(model_output)
#'
#' @export

convergence_check <- function(model_output) {

  # Extract model fit, summary and time
  fit <- model_output$stan_output
  summary <- data.frame(rstan::summary(fit)$summary)
  time <- max(rowSums(rstan::get_elapsed_time(fit)))

  # Print results
  printf("Time elapsed: %.1f", time)

  Rhat = sum(summary$Rhat > 1.1, na.rm = T)
  n_eff = floor(min(summary$n_eff))

  if (Rhat > 0) {
    printf("Model has **not** converged")
  }
  else {
    printf("Model has converged")
  }

  printf("Rhats > 1.1 =  %d", Rhat)
  printf("Minimum number of effective samples = %d, (%.2f/sec)",
         n_eff, n_eff/time)
}


#' Posterior check
#'
#' Runs posterior checks on model output by comparing the observed data to
#' posterior predictions under the same conditions. Presently, these are only
#' defined for tobit models (m1-m6) and include the root square mean error
#' (RMSE) of predictions after censoring, the total Euclidean distance between
#' all censored predictions and the data, the accuracy of predicted presences
#' and absences given as proportion of correct predictions, and the R-squared
#' statistic of predicted positive abundances against data.
#'
#' @usage posterior_check(model_output)
#'
#' @export

posterior_check <- function(model_output) {

  # Extract observed data
  y_obs <- model_output$data_list$y_observed

  # Extract posterior predictions
  y_pred <- rstan::extract(model_output$stan_output,
                           pars = "y_pred")$y_pred %>%
            apply(., c(2, 3), mean)

  # Censor predictions
  y_pred_censored = ifelse(y_pred < 0, 0, y_pred)

  # Set up dataframe, calculate presences/absences
  df <- data.frame(obs = matrix(y_obs, ncol = 1),
                   pred = matrix(y_pred, ncol = 1),
                   pred_censored = matrix(y_pred_censored, ncol = 1)) %>%
        mutate(obs_pa = ifelse(obs > 0, 1, 0),
               pred_pa = ifelse(pred > 0, 1, 0))

  # Root square mean error of censored predictions
  error <- df$obs - df$pred_censored
  rmse_censored <- round(sqrt(mean(error^2)), 2)
  printf("Root square mean error of censored predictions = %.3f", rmse_censored)

  # Total Euclidean distance between data and censored predictions
  dist <- sqrt(sum((df$obs - df$pred_censored)^2))
  printf("Distance between data and censored predictions = %.3f", dist)

  # Predicted presences
  pres_acc <-  filter(df, pred_pa == 1) %>%
    summarise(sum(obs_pa == 1) / n()) %>%
    round(2)

  printf("Accuracy of predicted presences = %.2f", pres_acc)

  # Predicted absences
  abs_acc <- filter(df, pred_pa == 0) %>%
    summarise(sum(obs_pa == 0) / n()) %>%
    round(2)

  printf("Accuracy of predicted absences = %.2f", abs_acc)

  # R-squared of predicted abundances
  abun_Rsq <- filter(df, pred > 0) %>%
    summarise(R_squared = 1 - (sum((obs - pred)^2) / sum((obs - mean(obs))^2)))

  printf("R-squared of predicted abundances = %.2f", abun_Rsq)
}

#' Extract parameters
#'
#' Extract samples from model output and calculate a summary statistic for
#' each parameter
#'
#' @param model_output output of a Stan model.
#' @param par vector of parameter names.
#' @param stat pass a summary statistic like mean, median or quantile.
#' @param ... extra arguments to summary statistic, eg probs = c(0.2, 0.8)
#'
#' @usage parameter <- extract_pars(model_output, par = "y_pred", stat = mean)

# extract_pars <- function(model_output,
#                          par = "y_pred",
#                          stat = mean) {
#
#   # Extract samples
#   samples <- rstan::extract(model, pars = par) %>%
#     unlist()
#
#   # Calculate point summary of posterior
#   parameters <- apply(pred, c(2, 3), mean)
# }

