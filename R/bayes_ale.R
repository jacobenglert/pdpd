#' Accumulated Local Effects for Bayesian Posterior Distributions
#'
#' Compute first and second-order accumulated local effects (ALE) with credible
#' intervals.
#'
#' @param x Data frame of covariates used in the model-fitting process.
#' @param f_hat Function that returns a posterior distribution of predictions given covariates \code{x}. The
#' function should return a matrix with \code{nrow(x)} rows and \code{num_save} columns.
#' @param vars Character vector of covariate names to vary in partial dependence calculations.
#' @param k Numeric vector of unique values to evaluate each covariate in \code{vars} at. By default this value will be decreased to the number of unique values if necessary.
#' @param limits Numeric vector of length 2 indicating the quantiles used for the lower and upper bounds of the credible interval.
#' @param f (Optional) function that returns the true prediction for a set of covariates. Useful when conducting simulation studies.
#' @param num_burn Number of posterior samples to skip when making predictions.
#' @param num_thin Thinning interval for posterior samples following \code{num_burn} samples.
#' @param include_main_effects Logical value indicating whether main effect ALEs should be added to
#' second-order ALE calculations. Has no effect when \code{length(vars)} is 1.
#' @param center Either a logical value indicating whether or not to vertically center the ALE
#' with respect to the average prediction, or 'median' to indicate centering should be done with
#' respect to the median x value(s).
#'
#'
#' @return A data frame of posterior mean estimates and credible intervals for
#' a predictor grid. For second-order ALE, plotting intervals are also returned
#' for custom plotting similar to \code{image}.
#' @export

bayes_ale <- function (x, f_hat, vars, k = 40, limits = c(0.025, 0.975),
                       num_burn = 0, num_thin = 1,
                       include_main_effects = TRUE, center = TRUE,
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

  if (!is.logical(center) & center != 'median') {
    stop("center must be a logical value or 'median'.")
  }

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

  # Compute ALE
  if (p == 1) {

    # Check variable type
    if (!(class(x[[vars[1]]]) %in% c("numeric","integer"))) {
      stop("bayes_ale can only take numeric or integer variables.")
    }

    fJ <- bayes_ale_1D(x, f_hat, vars, k, center)
    # x <- fJ$x
    fJ.post   <- fJ$ale
    fJ.mean   <- rowMeans(fJ.post)
    fJ.lower  <- rowQuants(fJ.post, limits[1])
    fJ.upper  <- rowQuants(fJ.post, limits[2])

    # Predict truth if f is provided
    if (!is.null(f)) {
      fJ.true <- ale_1D(x, f, vars, k, center)
    } else fJ.true <- rep(NA, length(fJ$x))

    # Return as nice data frame
    out <- data.frame(k = 1:length(fJ$x),
                      var = rep(vars, length(fJ$x)),
                      x = fJ$x,
                      est = fJ.mean,
                      lcl = fJ.lower,
                      ucl = fJ.upper,
                      truth = fJ.true,
                      n = fJ$n)

    return(out)

  } else if (p == 2) {

    # Check variable types
    if (!all(c(class(x[[vars[1]]]), class(x[[vars[2]]])) %in% c("numeric","integer"))) {
      stop("bayes_ale can only take numeric or integer variables.")
    }

    # Process input into k x k quantile intervals (k+1 x k+1 total points)
    z1 <- c(min(x[[vars[1]]]), as.numeric(stats::quantile(x[[vars[1]]], seq(1 / k, 1, length.out = k), type = 1)))
    z1 <- unique(z1)
    k1 <- length(z1) - 1
    a1 <- as.numeric(cut(x[[vars[1]]], breaks = z1, include.lowest = TRUE))
    z2 <- c(min(x[[vars[2]]]), as.numeric(stats::quantile(x[[vars[2]]], seq(1 / k, 1, length.out = k), type = 1)))
    z2 <- unique(z2)
    k2 <- length(z2) - 1
    a2 <- as.numeric(cut(x[[vars[2]]], breaks = z2, include.lowest = TRUE))

    # Create counterfactual datasets
    x11 <- x12 <- x21 <- x22 <- x
    x11[vars] <- cbind(z1[a1], z2[a2])
    x12[vars] <- cbind(z1[a1], z2[a2 + 1])
    x21[vars] <- cbind(z1[a1 + 1], z2[a2])
    x22[vars] <- cbind(z1[a1 + 1], z2[a2 + 1])

    # Make predictions
    y.hat11.pd <- f_hat(x11)
    y.hat12.pd <- f_hat(x12)
    y.hat21.pd <- f_hat(x21)
    y.hat22.pd <- f_hat(x22)

    if (!is.null(f)) {

      # Make predictions
      y11 <- f(x11)
      y12 <- f(x12)
      y21 <- f(x21)
      y22 <- f(x22)

    }

    rm(x11, x12, x21, x22)

    # Compute local effects
    delta.pd <- (y.hat22.pd - y.hat12.pd) - (y.hat21.pd - y.hat11.pd)
    rm(y.hat11.pd, y.hat12.pd, y.hat21.pd, y.hat22.pd)
    delta.pd.list <- as.list(asplit(delta.pd, 2))
    avg.delta.pd.list <- lapply(delta.pd.list, \(x) as.matrix(tapply(x, list(a1, a2), mean)))
    rm(delta.pd.list)

    # For non-existent regions, use nearest neighbor
    na.delta <- is.na(avg.delta.pd.list[[1]])
    na.ind <- which(na.delta, arr.ind = TRUE, useNames = FALSE)
    if(nrow(na.ind) > 0){
      not.na.ind <- which(!na.delta, arr.ind = TRUE, useNames = FALSE)
      range1 <- max(z1) - min(z1)
      range2 <- max(z2) - min(z2)
      z.na <- cbind((z1[na.ind[, 1]] + z1[na.ind[, 1] + 1]) / 2 / range1,
                    (z2[na.ind[, 2]] + z2[na.ind[, 2] + 1]) / 2 / range2)
      z.not.na <- cbind((z1[not.na.ind[, 1]] + z1[not.na.ind[, 1] + 1]) / 2 / range1,
                        (z2[not.na.ind[, 2]] + z2[not.na.ind[, 2] + 1]) / 2 / range2)
      nbrs <- yaImpute::ann(z.not.na, z.na, k = 1, verbose = F)$knnIndexDist[, 1]
      avg.delta.pd.list <- lapply(avg.delta.pd.list, \(x){
        x[na.ind] <- x[matrix(not.na.ind[nbrs,], ncol = 2)]
        return(x)
      })
    }

    # Accumulate local effects
    fJ <- lapply(avg.delta.pd.list, \(x) apply(t(apply(x, 1, cumsum)), 2, cumsum))
    fJ <- lapply(fJ, \(x) rbind(rep(0, k2), x))
    fJ <- lapply(fJ, \(x) cbind(rep(0, k1 + 1), x))

    # Center ALE
    b <- as.matrix(table(a1, a2)) # cell counts
    b1 <- apply(b, 1, sum) # row cell counts
    b2 <- apply(b, 2, sum) # column cell counts

    # Compute first-order ALEs
    delta     <- lapply(fJ, \(x) x[2:(k1 + 1), ] - x[1:k1, ]) # local effects
    b.delta   <- lapply(delta, \(x) b * (x[, 1:k2] + x[, 2:(k2 + 1)]) / 2) # weighted average numerator
    avg.delta <- lapply(b.delta, \(x) rowSums(x) / b1) # weighted average calculation
    fJ1       <- lapply(avg.delta, \(x) c(0, cumsum(x))) # accumulate

    delta     <- lapply(fJ, \(x) x[, 2:(k2 + 1)] - x[, 1:k2])
    b.delta   <- lapply(delta, \(x) b * (x[1:k1, ] + x[2:(k1 + 1), ]) / 2)
    avg.delta <- lapply(b.delta, \(x) colSums(x) / b2)
    fJ2       <- lapply(avg.delta, \(x) c(0, cumsum(x)))

    # Subtract first-order ALEs from second-order ALE
    fJ <- mapply(\(x, y, z) x - outer(y, rep(1, k2 + 1)) - outer(rep(1, k1 + 1), z),
                 x = fJ, y = fJ1, z = fJ2, SIMPLIFY = FALSE)

    # Optional centering of the effects
    if (is.logical(center) & center) {
      fJ <- lapply(fJ, \(x) x - sum(b * (x[1:k1, 1:k2] +
                                           x[1:k1, 2:(k2 + 1)] +
                                           x[2:(k1 + 1), 1:k2] +
                                           x[2:(k1 + 1), 2:(k2 + 1)]) / 4) / sum(b))
    } else if (center == 'median') {
      m1 <- max(which(stats::median(x[[vars[1]]]) > z1))
      m2 <- max(which(stats::median(x[[vars[1]]]) > z2))
      fJ <- lapply(fJ, \(x) x - ((x[m1,m2] + x[m1,m2+1] + x[m1+1,m2] + x[m1+1,m2+1]) / 4))
    }

    # If main effects requested, add back the marginally centered main effects
    if (include_main_effects) {

      # Compute centered ALE main effects
      ale1 <- bayes_ale_1D(x, f_hat, vars[1], k, center)$ale
      ale2 <- bayes_ale_1D(x, f_hat, vars[2], k, center)$ale

      # Convert to list format
      fJ1 <- fJ2 <- vector("list", ncol(ale1))
      for (i in 1:ncol(ale1)) {
        fJ1[[i]] <- ale1[,i]
        fJ2[[i]] <- ale2[,i]
      }

      # Add ALE main effects back to the centered second-order predictions
      fJ <- mapply(\(x, y, z) x + outer(y, rep(1, k2 + 1)) + outer(rep(1, k1 + 1), z),
                   x = fJ, y = fJ1, z = fJ2, SIMPLIFY = FALSE)
    }

    # Predict truth if f is provided
    if (!is.null(f)) {

      # Compute local effects
      delta <- (y22 - y12) - (y21 - y11)
      avg.delta <- as.matrix(tapply(delta, list(a1, a2), mean))
      if (nrow(na.ind) > 0) {
        avg.delta[na.ind] <- avg.delta[matrix(not.na.ind[nbrs, ], ncol = 2)]
      }

      # Accumulate local effects
      fJ.true <- apply(t(apply(avg.delta, 1, cumsum)), 2, cumsum)
      fJ.true <- rbind(rep(0, k2), fJ.true)
      fJ.true <- cbind(rep(0, k1 + 1), fJ.true)

      # Center ALE
      # Compute first-order true ALEs
      delta <- fJ.true[2:(k1 + 1), ] - fJ.true[1:k1, ]
      b.delta <- b * (delta[, 1:k2] + delta[, 2:(k2 + 1)]) / 2
      avg.delta <- rowSums(b.delta) / b1
      fJ1.true <- c(0, cumsum(avg.delta))

      delta     <- fJ.true[, 2:(k2 + 1)] - fJ.true[, 1:k2]
      b.delta   <- b * (delta[1:k1, ] + delta[2:(k1 + 1), ]) / 2
      avg.delta <- colSums(b.delta) / b2
      fJ2.true <- c(0, cumsum(avg.delta))

      # Subtract first-order true ALEs from second-order true ALE
      fJ.true <- fJ.true - outer(fJ1.true, rep(1, k2 + 1)) - outer(rep(1, k1 + 1), fJ2.true)

      # Optional centering
      if (is.logical(center) & center) {
        fJ.true <- fJ.true - sum(b * (fJ.true[1:k1, 1:k2] +
                                        fJ.true[1:k1, 2:(k2 + 1)] +
                                        fJ.true[2:(k1 + 1), 1:k2] +
                                        fJ.true[2:(k1 + 1), 2:(k2 + 1)]) / 4) / sum(b)
      } else if (center == 'median') {
        fJ.true <- fJ.true - (fJ.true[m1,m2] + fJ.true[m1,m2+1] + fJ.true[m1+1,m2] + fJ.true[m1+1,m2+1]) / 4
      }

      # If main effects requested, add back the marginally centered main effects
      if (include_main_effects) {

        # Get true ALE main effect functions
        fJ1.true <- ale_1D(x, f, vars[1], k, center)
        fJ2.true <- ale_1D(x, f, vars[2], k, center)

        # Add true ALE main effects back to true ALE second-order effects
        fJ.true <- fJ.true + outer(fJ1.true, rep(1, k2 + 1)) + outer(rep(1, k1 + 1), fJ2.true)
      }

    } else fJ.true <- NULL

    # Return as nice data frame
    fJ.array <- array(unlist(fJ), dim = c(dim(fJ[[1]]), length(fJ)))
    fJ.mean <- apply(fJ.array, 2, rowMeans)
    fJ.lower <- apply(fJ.array, 2, rowQuants, q = limits[1])
    fJ.upper <- apply(fJ.array, 2, rowQuants, q = limits[2])

    # Lower and upper bounds for plotting
    #   - Plotting coordinates: ((z1 - w1, z1 + w2, z2 - h1, z2 + h1))
    w1 <- c(diff(z1)[1], diff(z1)) / 2
    w2 <- c(rev(abs(diff(rev(z1)))), rev(abs(diff(z1)))[1]) / 2
    h1 <- c(diff(z2)[1], diff(z2)) / 2
    h2 <- c(rev(abs(diff(rev(z2)))), rev(abs(diff(z2)))[1]) / 2

    out <- data.frame(matrix(nrow = length(z1) * length(z2), ncol = 15))
    colnames(out) <- c('k1','k2','var1','var2','x1','x2',
                       'est','lcl','ucl','truth','n',
                       'w1','w2','h1','h2')
    index <- 1
    for (i in 1:length(z1)) {
      for (j in 1:length(z2)) {
        out$k1[index] <- i
        out$k2[index] <- j
        out$x1[index] <- z1[i]
        out$x2[index] <- z2[j]
        out$est[index] <- fJ.mean[i,j]
        out$lcl[index] <- fJ.lower[i,j]
        out$ucl[index] <- fJ.upper[i,j]
        if (!is.null(fJ.true)) out$truth[index] <- fJ.true[i,j]
        out$w1[index] <- w1[i]
        out$w2[index] <- w2[i]
        out$h1[index] <- h1[j]
        out$h2[index] <- h2[j]

        # Record sample size used
        # We always start accumulating from 0 in both directions, so the first
        # point doesn't use data. But we will just input the sample size for the
        # second point.
        out$n[index] <- b[max(i-1, 1), max(j-1, 1)]

        index <- index + 1
      }
    }
    out$var1 <- vars[1]
    out$var2 <- vars[2]

    return(out)

  } else {
    stop("bayes_ale can only take 1 or 2 variables.")
  }
}


bayes_ale_1D <- function (x, f_hat, var, k = 40, center = TRUE) {

  # Process input into k quantile intervals (k+1 total points)
  z <- c(min(x[[var]]), as.numeric(stats::quantile(x[[var]], seq(1 / k, 1, length.out = k), type = 1)))
  z <- unique(z)
  k <- length(z) - 1 # update k if necessary
  a1 <- as.numeric(cut(x[[var]], breaks = z, include.lowest = TRUE))

  # Create counterfactual datasets
  x1 <- x2 <- x
  x1[[var]] <- z[a1]
  x2[[var]] <- z[a1 + 1]

  # Make predictions
  y.hat1.pd <- f_hat(x1)
  y.hat2.pd <- f_hat(x2)
  rm(x1, x2)

  # Compute individual local effects
  delta.pd <- y.hat2.pd - y.hat1.pd
  rm(y.hat1.pd, y.hat2.pd)

  # Compute average local effects
  avg.delta.pd <- apply(delta.pd, 2, \(x) tapply(x, a1, mean))

  # Accumulate average local effects
  fJ <- apply(avg.delta.pd, 2, \(x) c(0, cumsum(x)))

  # Center ALE
  b1 <- as.numeric(table(a1))
  if (is.logical(center)) {
    if (center) {
      fJ <- sweep(fJ, 2, apply(fJ, 2, \(x) sum((x[1:k] + x[2:(k + 1)]) / 2 * b1) / sum(b1)), '-')
    }
  } else if (center == 'median') {
    m <- max(which(stats::median(x[[var]]) > z))
    fJ <- sweep(fJ, 2, apply(fJ, 2, \(x) (x[m] + x[m+1]) / 2), '-')
  } else {
    stop("center must be a logical value or 'median'.")
  }

  return (list(ale = fJ, x = z, n = c(b1[1], b1)))
}

ale_1D <- function (x, f, var, k = 40, center = TRUE) {

  # Process input into k quantile intervals (k+1 total points)
  z <- c(min(x[[var]]), as.numeric(stats::quantile(x[[var]], seq(1 / k, 1, length.out = k), type = 1)))
  z <- unique(z)
  k <- length(z) - 1 # update k if necessary
  a1 <- as.numeric(cut(x[[var]], breaks = z, include.lowest = TRUE))

  # Create counterfactual datasets
  x1 <- x2 <- x
  x1[[var]] <- z[a1]
  x2[[var]] <- z[a1 + 1]

  # Make predictions
  y1 <- f(x1)
  y2 <- f(x2)
  rm(x1, x2)

  # Compute individual local effects
  delta <- y2 - y1
  rm(y1, y2)
  avg.delta <- tapply(delta, a1, mean)

  # Accumulate local effects
  fJ <- c(0, cumsum(avg.delta))

  # Center ALE
  b1 <- as.numeric(table(a1))
  if (is.logical(center)) {
    if (center) {
      fJ <- fJ - sum((fJ[1:k] + fJ[2:(k + 1)]) / 2 * b1) / sum(b1)
    }
  } else if (center == 'median') {
    m <- max(which(stats::median(x[[var]]) > z))
    fJ <- fJ - ((fJ[m] + fJ[m+1]) / 2)
  } else {
    stop("center must be a logical value or 'median'.")
  }

  return (as.numeric(fJ))
}

