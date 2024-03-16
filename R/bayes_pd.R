#' Partial Dependence for Bayesian Posterior Distributions
#'
#' Compute partial dependence statistics with credible intervals.
#'
#' @param x Data frame of covariates used in the model-fitting process.
#' @param f_hat Function that returns a posterior distribution of predictions given covariates \code{x}. The
#' function should return a matrix with \code{nrow(x)} rows and \code{num_save} columns.
#' @param vars Character vector of covariate names to vary in partial dependence calculations.
#' @param k Numeric vector of unique values to evaluate each covariate in \code{vars} at. By default this value will be decreased to the number of unique values if necessary.
#' @param k_vals Alternatively, a list of \code{length(vars)} numeric vectors, each representing the unique values to evaluate the corresponding variable at.
#' @param limits Numeric vector of length 2 indicating the quantiles used for the lower and upper bounds of the credible interval.
#' @param f (Optional) function that returns the true prediction for a set of covariates. Useful when conducting simulation studies.
#' @param num_burn Number of posterior samples to skip when making predictions.
#' @param num_thin Thinning interval for posterior samples following \code{num_burn} samples.
#'
#' @return A data frame of posterior mean estimates and credible intervals for
#' a predictor grid.
#' @export

bayes_pd <- function (x, f_hat, vars,
                      k, k_vals = NULL,
                      limits = c(0.025, 0.975),
                      num_burn = 0, num_thin = 1,
                      f = NULL) {

  # Check limits
  if ((length(limits) != 2) | !is.numeric(limits)) {
    stop("limits must be a numeric vector of length 2.")
  }

  if (!all(limits < 1) | !all(limits > 0) | (limits[1] >= limits[2])) {
    stop("limits must lie in the interval (0,1), and the second element must be larger than the first.")
  }

  if (num_burn %% 1 != 0 | num_burn < 0) stop("num_burn must be a non-negative integer.")
  if (num_thin %% 1 != 0 | num_thin < 1) stop("num_burn must be a positive integer.")


  # Identify samples to retain in prediction
  if (num_burn == 0 & num_thin == 1) {
    f_hat_post <- f_hat
  } else {
    num_samples <- ncol(f_hat(x[1:2,])) # fix this later so it works for 1 observation
    if (num_burn >= num_samples) {
      stop("Requested burn-in period exceeds the number of available samples.")
    } else if (num_burn + num_thin >= num_samples) {
      stop("Requested burn-in period and thinning exceeds the number of available samples.")
    } else {
      # Samples to retain
      good_samples <- seq.int(from = num_burn + num_thin, to = num_samples, by = num_thin)

      # Updated prediction function
      f_hat_post <- function (x) f_hat(x)[,good_samples]
    }
  }

  # Process input data
  x <- as.data.frame(x)
  p <- length(vars)

  # Determine grid intervals
  if (!is.null(k_vals)) {
    if ((!is.list(k_vals) & !is.numeric(k_vals)) |
        (is.list(k_vals) & !(length(k_vals) %in% c(1, p)))) {
      stop("k_vals must be supplied as a list of numeric vectors for each variable supplied in vars.")
    }
    # Recycle k_vals if only one set are given
    if (length(k_vals) == 1 | is.numeric(k_vals)) {
      if (p > 1) message("Using the same k_vals for each variable in vars.")
      if (is.numeric(k_vals)) k_vals <- replicate(p, k_vals, simplify = FALSE)
      if (is.list(k_vals)) k_vals <- rep(k_vals, p)
    }

    # Send a warning if supplied intervals extrapolate from the data
    for (j in seq_len(p)) {
      v <- x[, vars[j], drop = TRUE]
      if (any(k_vals[[j]] < min(v)) | any(k_vals[[j]] > max(v))) {
        message(paste0("Interval requested for variable ", vars[j], " extrapolates from the data."))
      }
    }

  } else {

    # Create equally spaced intervals if only number of intervals is provided
    if (k %% 1 != 0) stop("K must be a positive integer.")
    if (length(k) == 1) k <- rep(k, p)
    k_vals <- list()
    for (j in seq_len(p)) {
      var <- vars[j]
      v <- x[, var, drop = TRUE]
      n_unique <- min(length(unique(v)), k[j])
      k_vals[[var]] <- seq(min(v), max(v), length.out = n_unique)
    }
  }

  # Create named grid
  names(k_vals) <- colnames(x[vars])
  k.grid <- expand.grid(k_vals)

  # Calculate partial dependence for each grid cell
  for (row in seq_len(nrow(k.grid))) {
    pdata <- x
    for (var in vars) pdata[[var]] <- k.grid[[var]][row]
    preds <- colMeans(f_hat_post(pdata))
    k.grid$est[row] <- mean(preds)
    k.grid$lcl[row] <- stats::quantile(preds, limits[1])
    k.grid$ucl[row] <- stats::quantile(preds, limits[2])

    # Calculate the true partial dependence if provided
    if (!is.null(f)) k.grid$truth[row] <- mean(f(pdata))
  }

  return(k.grid)

}
