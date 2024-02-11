
#' Accumulated Local Effects for a Bayesian Posterior
#'
#' @description
#' This function provides pointwise first- and second-order ALE with uncertainty
#' intervals using the posterior samples of predictions from a Bayesian model. The
#' function currently only supports computing ALEs for numeric predictors.
#'
#'
#' @param X The data frame of predictor variables to which the model was fit.
#' @param model The fitted Bayesian model object from which predictions can be made.
#' @param pred_fun A user-supplied function that will be used to predict the response
#' for \code{model} for some specified input. \code{pred_fun} has two arguments. The first argument
#' is named \code{model} and the second argument is named \code{newdata}. The output of \code{pred_fun}
#' when applied to \code{model} and \code{X} must be a numeric matrix \code{nrow(X)} rows and number of
#' columns equal to the number of posterior samples drawn (this is the main difference
#' between this function and \code{ALEPlot} from the \href{https://cran.r-project.org/web//packages//ALEPlot/ALEPlot.pdf}{\code{ALEPlot}} package).
#' @param J same as in \href{https://cran.r-project.org/web//packages//ALEPlot/ALEPlot.pdf}{\code{ALEPlot}}.
#' @param K same as in \href{https://cran.r-project.org/web//packages//ALEPlot/ALEPlot.pdf}{\code{ALEPlot}}.
#' @param CrI a numeric vector of length two indicating the pointwise quantiles
#' to return in addition to the point estimate.
#' @param f_true an optional function which maps \code{X} to its true predicted value
#' (useful for comparing results in a simulation study).
#' @param include_main_effects logical determining whether the centered main effect
#' ALEs should be added back to the second-order ALE (has no effect when \code{dim(J) = 1})
#' @param center centring method to be used. Must be one of 'avg_pred' or 'median'.
#' By default, ALE is centered with respect to the predictor value(s) that corresponds
#' to the average (mean) ALE prediction. Setting to 'median' will have the effect
#' of centering the ALE about the ALE prediction for the median predictor value(s).
#'
#' @details
#' See the Apley (2020) reference below for a background on ALE.
#'
#' @references
#' Daniel W. Apley, Jingyu Zhu, Visualizing the Effects of Predictor Variables in
#' Black Box Supervised Learning Models, Journal of the Royal Statistical Society
#' Series B: Statistical Methodology, Volume 82, Issue 4, September 2020, Pages
#' 1059–1086, \url{https://doi.org/10.1111/rssb.12377}
#'
#' @export
#'
mcmc_ale <- function (X, model, pred_fun, J, K = 40,
                      CrI = c(0.025, 0.975), f_true = NULL, include_main_effects = TRUE,
                      center = 'avg_pred') {

  # Get dimensions of data
  N <- dim(X)[1]
  d <- dim(X)[2]

  if(length(J) == 1){

    stopifnot("mcmc_ale currently only works with numeric variables" = class(X[, J[1]]) %in% c("numeric", "integer"))

    fJ <- get_ale_post_1D(X, model, pred_fun, J, K, center)
    x <- fJ$x
    fJ.post   <- fJ$ale
    fJ.mean   <- rowMeans(fJ.post)
    fJ.lower  <- apply(fJ.post, 1, \(x) stats::quantile(x, CrI[1]))
    fJ.upper  <- apply(fJ.post, 1, \(x) stats::quantile(x, CrI[2]))
    if (!is.null(f_true)) fJ.true <- get_ale_true_1D(X, f_true, J, K, center)
    else fJ.true <- rep(NA, length(x))

    # Return as nice data frame
    out <- data.frame(k = 1:length(x),
                      var = rep(colnames(X)[J], length(x)),
                      x = x,
                      est = fJ.mean,
                      lcl = fJ.lower,
                      ucl = fJ.upper,
                      truth = fJ.true)
    return(out)


  } else if(length(J) == 2) {

    stopifnot("mcmc_ale currently only works with numeric variables" = class(X[, J[1]]) %in% c("numeric", "integer") & class(X[, J[2]]) %in% c("numeric", "integer"))

    # Process inputs into intervals
    z1 <- c(min(X[, J[1]]), as.numeric(stats::quantile(X[, J[1]], seq(1 / K, 1, length.out = K), type = 1)))
    z1 <- unique(z1)
    K1 <- length(z1) - 1
    a1 <- as.numeric(cut(X[, J[1]], breaks = z1, include.lowest = TRUE))
    z2 <- c(min(X[, J[2]]), as.numeric(stats::quantile(X[, J[2]], seq(1 / K, 1, length.out = K), type = 1)))
    z2 <- unique(z2)
    K2 <- length(z2) - 1
    a2 <- as.numeric(cut(X[, J[2]], breaks = z2, include.lowest = TRUE))

    # Create counterfactual datasets
    X11 <- X12 <- X21 <- X22 <- X
    X11[, J] <- cbind(z1[a1], z2[a2])
    X12[, J] <- cbind(z1[a1], z2[a2 + 1])
    X21[, J] <- cbind(z1[a1 + 1], z2[a2])
    X22[, J] <- cbind(z1[a1 + 1], z2[a2 + 1])

    # Make predictions
    y.hat11.mat <- pred_fun(model = model, newdata = X11)
    y.hat12.mat <- pred_fun(model = model, newdata = X12)
    y.hat21.mat <- pred_fun(model = model, newdata = X21)
    y.hat22.mat <- pred_fun(model = model, newdata = X22)

    # Compute local effects
    Delta.mat <- (y.hat22.mat - y.hat21.mat) - (y.hat12.mat - y.hat11.mat)
    Delta.list <- asplit(Delta.mat, 2) |> as.list()
    Delta.mat.list <- lapply(Delta.list, \(x) as.matrix(tapply(x, list(a1, a2), mean)))

    # For non-existent regions, use nearest neighbor
    NA.Delta <- is.na(Delta.mat.list[[1]])
    NA.ind <- which(NA.Delta, arr.ind = TRUE, useNames = FALSE)
    if(nrow(NA.ind) > 0){
      notNA.ind <- which(!NA.Delta, arr.ind = TRUE, useNames = FALSE)
      range1 <- max(z1) - min(z1)
      range2 <- max(z2) - min(z2)
        Z.NA <- cbind((z1[NA.ind[, 1]] + z1[NA.ind[, 1] + 1]) / 2 / range1,
                    (z2[NA.ind[, 2]] + z2[NA.ind[, 2] + 1]) / 2 / range2)
      Z.notNA <- cbind((z1[notNA.ind[, 1]] + z1[notNA.ind[, 1] + 1]) / 2 /range1,
                       (z2[notNA.ind[, 2]] + z2[notNA.ind[, 2] + 1])/2/range2)
        nbrs <- yaImpute::ann(Z.notNA, Z.NA, k = 1, verbose = F)$knnIndexDist[, 1]
      Delta.mat.list <- lapply(Delta.mat.list,
                               \(x){
                                 x[NA.ind] <- x[matrix(notNA.ind[nbrs,], ncol = 2)]
                                 return(x)
                                 })
    }

    # Accumulate local effects
    fJ <- lapply(Delta.mat.list, \(x) apply(t(apply(x, 1, cumsum)), 2, cumsum))
    fJ <- lapply(fJ, \(x) rbind(rep(0, K2), x))
    fJ <- lapply(fJ, \(x) cbind(rep(0, K1 + 1), x))

    # Center ALE
    b <- as.matrix(table(a1, a2))
    b1 <- apply(b, 1, sum)
    b2 <- apply(b, 2, sum)

    Delta <- lapply(fJ, \(x) x[2:(K1 + 1), ] - x[1:K1, ])
    b.Delta <- lapply(Delta, \(x) b * (x[, 1:K2] + x[, 2:(K2 + 1)]) / 2)
    Delta.Ave <- lapply(b.Delta, \(x) apply(x, 1, sum) / b1)
    fJ1 <- lapply(Delta.Ave, \(x) c(0, cumsum(x)))
    Delta <- lapply(fJ, \(x) x[, 2:(K2 + 1)] - x[, 1:K2])
    b.Delta <- lapply(Delta, \(x) b * (x[1:K1, ] + x[2:(K1 + 1), ]) / 2)
    Delta.Ave <- lapply(b.Delta, \(x) apply(x, 2, sum) / b2)
    fJ2 <- lapply(Delta.Ave, \(x) c(0, cumsum(x)))
    fJ <- mapply(\(x, y, z) x - outer(y, rep(1, K2 + 1)) - outer(rep(1, K1 + 1), z),
                 x = fJ, y = fJ1, z = fJ2, SIMPLIFY = FALSE)


    if (center == 'avg_pred') {
      fJ <- lapply(fJ, \(x) x - sum(b * (x[1:K1, 1:K2] + x[1:K1, 2:(K2 + 1)] + x[2:(K1 + 1), 1:K2] + x[2:(K1 + 1), 2:(K2 + 1)]) / 4) / sum(b))
    }
    else if (center == 'median') {
      M1 <- max(which(stats::median(X[, J[1]]) > z1))
      M2 <- max(which(stats::median(X[, J[2]]) > z2))
      fJ <- lapply(fJ, \(x) x - ((x[M1,M2] + x[M1,M2+1] + x[M1+1,M2] + x[M1+1,M2+1]) / 4))
    }

    # If main effects requested, add back the marginally centered main effects
    if (include_main_effects) {

      # Compute centered ALE main effects
      ale1 <- get_ale_post_1D(X, model, pred_fun, J[1], K, center)$ale
      ale2 <- get_ale_post_1D(X, model, pred_fun, J[2], K, center)$ale

      # Convert to list format
      fJ1 <- fJ2 <- vector("list", ncol(ale1))
      for (i in 1:ncol(ale1)) {
        fJ1[[i]] <- ale1[,i]
        fJ2[[i]] <- ale2[,i]
      }

      # Add ALE main effects back to the centered second-order predictions
      fJ <- mapply(\(x, y, z) x + outer(y, rep(1, K2 + 1)) + outer(rep(1, K1 + 1), z),
                   x = fJ, y = fJ1, z = fJ2, SIMPLIFY = FALSE)
    }

    # Return true ALE function if true function is provided
    if (!is.null(f_true)) {

      # Make predictions
      y.hat11.true <- f_true(X11)
      y.hat12.true <- f_true(X12)
      y.hat21.true <- f_true(X21)
      y.hat22.true <- f_true(X22)

      # Compute local effects
      Delta.true <- (y.hat22.true - y.hat21.true) - (y.hat12.true - y.hat11.true)
      Delta.true <- as.matrix(tapply(Delta.true, list(a1, a2), mean))
      if (nrow(NA.ind) > 0) {
        Delta.true[NA.ind] <- Delta.true[matrix(notNA.ind[nbrs, ], ncol = 2)]
      }

      # Accumulate local effects
      fJ.true <- apply(t(apply(Delta.true, 1, cumsum)), 2, cumsum)
      fJ.true <- rbind(rep(0, K2), fJ.true)
      fJ.true <- cbind(rep(0, K1 + 1), fJ.true)

      # Center ALE
      Delta.true <- fJ.true[2:(K1 + 1), ] - fJ.true[1:K1, ]
      b.Delta.true <- b * (Delta.true[, 1:K2] + Delta.true[, 2:(K2 + 1)]) / 2
      Delta.Ave.true <- apply(b.Delta.true, 1, sum) / b1
      fJ1.true <- c(0, cumsum(Delta.Ave.true))
      Delta.true <- fJ.true[, 2:(K2 + 1)] - fJ.true[, 1:K2]
      b.Delta.true <- b * (Delta.true[1:K1, ] + Delta.true[2:(K1 + 1), ]) / 2
      Delta.Ave.true <- apply(b.Delta.true, 2, sum) / b2
      fJ2.true <- c(0, cumsum(Delta.Ave.true))
      fJ.true <- fJ.true - outer(fJ1.true, rep(1, K2 + 1)) - outer(rep(1, K1 + 1), fJ2.true)

      if (center == 'avg_pred') {
        fJ0.true <- sum(b * (fJ.true[1:K1, 1:K2] +
                             fJ.true[1:K1, 2:(K2 + 1)] +
                             fJ.true[2:(K1 + 1), 1:K2] +
                             fJ.true[2:(K1 + 1), 2:(K2 + 1)]) / 4) / sum(b)
      }
      else if (center == 'median') {
        fJ0.true <- (fJ.true[M1,M2] + fJ.true[M1,M2+1] + fJ.true[M1+1,M2] + fJ.true[M1+1,M2+1]) / 4
      }

      fJ.true <- fJ.true - fJ0.true

      # If main effects requested, add back the marginally centered main effects
      if (include_main_effects) {

        # Get true ALE main effect functions
        fJ1.true <- get_ale_true_1D(X, f_true, J[1], K, center)
        fJ2.true <- get_ale_true_1D(X, f_true, J[2], K, center)

        # Add true ALE main effects back to true ALE second-order effects
        fJ.true <- fJ.true + outer(fJ1.true, rep(1, K2 + 1)) + outer(rep(1, K1 + 1), fJ2.true)
      }

    } else fJ.true <- NULL


    # Format output
    # x <- list(z1, z2)
    # K <- c(K1, K2)

    fJ.array <- array(unlist(fJ), dim = c(dim(fJ[[1]]), length(fJ)))
    fJ.mean  <- apply(fJ.array, c(1, 2), mean)
    fJ.lower <- apply(fJ.array, c(1, 2), \(x) stats::quantile(x, CrI[1]))
    fJ.upper <- apply(fJ.array, c(1, 2), \(x) stats::quantile(x, CrI[2]))


    # Return as nice data frame
    # Lower and upper bounds for plotting
    #   - Plotting coordinates: ((z1 - w1, z1 + w2, z2 - h1, z2 + h1))
    w1 <- c(diff(z1)[1], diff(z1)) / 2
    w2 <- c(rev(abs(diff(rev(z1)))), rev(abs(diff(z1)))[1]) / 2
    h1 <- c(diff(z2)[1], diff(z2)) / 2
    h2 <- c(rev(abs(diff(rev(z2)))), rev(abs(diff(z2)))[1]) / 2

    out <- data.frame(matrix(nrow = length(z1) * length(z2), ncol = 14))
    colnames(out) <- c('k1','k2','var1','var2','x1','x2',
                       'est','lcl','ucl','truth',
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
        index <- index + 1
      }
    }
    out$var1 <- colnames(X)[J[1]]
    out$var2 <- colnames(X)[J[2]]

    return(out)

  } else{
    print("error:  J must be a vector of length one or two")
  }

}
