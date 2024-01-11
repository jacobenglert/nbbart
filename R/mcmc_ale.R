
#' @export
mcmc_ale <- function(X, model, pred_fun, J, K = 40, CrI = c(0.025, 0.975), f_true = NULL){

  N <- dim(X)[1]
  d <- dim(X)[2]

  if(length(J) == 1){

    stopifnot("mcmc_ale currently only works with numeric variables" = class(X[, J[1]]) %in% c("numeric", "integer"))

    # Process input into intervals
    z <- c(min(X[, J]), as.numeric(stats::quantile(X[, J], seq(1 / K, 1, length.out = K), type = 1)))
    z <- unique(z)
    K <- length(z) - 1
    a1 <- as.numeric(cut(X[, J], breaks = z, include.lowest = TRUE))

    # Create counterfactual datasets
    X1 <- X2 <- X
    X1[, J] <- z[a1]
    X2[, J] <- z[a1 + 1]

    # Make predictions
    y.hat1.mat <- pred_fun(model = model, newdata = X1)
    y.hat2.mat <- pred_fun(model = model, newdata = X2)

    # Compute local effects
    Delta.mat <- y.hat2.mat - y.hat1.mat
    Delta.mat <- apply(Delta.mat, 2, \(x) tapply(x, a1, mean))

    # Accumulate local effects
    fJ.mat <- apply(Delta.mat, 2, \(x) c(0, cumsum(x)))

    # Center ALE
    b1 <- as.numeric(table(a1))
    fJ.mat.cen <- sweep(fJ.mat, 2, apply(fJ.mat, 2, \(x) sum((x[1:K] + x[2:(K + 1)]) / 2 * b1) / sum(b1)), '-')

    # Format output
    x <- z
    fJ.mean <- rowMeans(fJ.mat.cen)
    fJ.lower <- apply(fJ.mat.cen, 1, \(x) stats::quantile(x, CrI[1]))
    fJ.upper <- apply(fJ.mat.cen, 1, \(x) stats::quantile(x, CrI[2]))

    if(!is.null(f_true)){
      # Make predictions
      y.hat1.true <- f_true(X1)
      y.hat2.true <- f_true(X2)

      # Compute local effects
      Delta.true <- y.hat2.true - y.hat1.true
      Delta.true <- tapply(Delta.true, a1, mean)

      # Accumulate local effects
      fJ.true <- c(0, cumsum(Delta.true))

      # Center ALE
      fJ.true <- fJ.true - sum((fJ.true[1:K] + fJ.true[2:(K + 1)])/2 * b1) / sum(b1)

    } else fJ.true <- NULL


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
    NA.ind <- which(NA.Delta, arr.ind = T, useNames = F)
    if(nrow(NA.ind) > 0){
      notNA.ind <- which(!NA.Delta, arr.ind = T, useNames = F)
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
    fJ <- lapply(fJ, \(x) x - sum(b * (x[1:K1, 1:K2] + x[1:K1, 2:(K2 + 1)] + x[2:(K1 + 1), 1:K2] + x[2:(K1 + 1), 2:(K2 + 1)]) / 4) / sum(b))

    # Format output
    x <- list(z1, z2)
    K <- c(K1, K2)

    fJ.array <- array(unlist(fJ), dim = c(dim(fJ[[1]]), length(fJ)))
    fJ.mean <- apply(fJ.array, c(1, 2), mean)
    fJ.lower <- apply(fJ.array, c(1, 2), \(x) stats::quantile(x, CrI[1]))
    fJ.upper <- apply(fJ.array, c(1, 2), \(x) stats::quantile(x, CrI[2]))

    # contour(x[[1]], x[[2]], fJ, add = TRUE, drawlabels = TRUE)
    # if (NA.plot == FALSE) {
    #   if (nrow(NA.ind) > 0) {
    #     rect(xleft = z1[NA.ind[, 1]],
    #          ybottom = z2[NA.ind[,2]],
    #          xright = z1[NA.ind[, 1] + 1],
    #          ytop = z2[NA.ind[,2] + 1],
    #          col = "black")
    #   }
    # }

    if(!is.null(f_true)){
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
      fJ.true <- apply(t(apply(Delta.true, 1, cumsum)), 2, cumsum)
      fJ.true <- rbind(rep(0, K2), fJ.true)
      fJ.true <- cbind(rep(0, K1 + 1), fJ.true)
      Delta.true <- fJ.true[2:(K1 + 1), ] - fJ.true[1:K1, ]
      b.Delta.true <- b * (Delta.true[, 1:K2] + Delta.true[, 2:(K2 + 1)]) / 2
      Delta.Ave.true <- apply(b.Delta.true, 1, sum) / b1
      fJ1.true <- c(0, cumsum(Delta.Ave.true))
      Delta.true <- fJ.true[, 2:(K2 + 1)] - fJ.true[, 1:K2]
      b.Delta.true <- b * (Delta.true[1:K1, ] + Delta.true[2:(K1 + 1), ]) / 2
      Delta.Ave.true <- apply(b.Delta.true, 2, sum) / b2
      fJ2.true <- c(0, cumsum(Delta.Ave.true))
      fJ.true <- fJ.true - outer(fJ1.true, rep(1, K2 + 1)) - outer(rep(1, K1 + 1), fJ2.true)
      fJ0.true <- sum(b * (fJ.true[1:K1, 1:K2] +
                           fJ.true[1:K1, 2:(K2 + 1)] +
                           fJ.true[2:(K1 + 1), 1:K2] +
                           fJ.true[2:(K1 + 1), 2:(K2 + 1)]) / 4) / sum(b)
      fJ.true <- fJ.true - fJ0.true

    } else fJ.true <- NULL

  } else{
    print("error:  J must be a vector of length one or two")
  }

  return(list(K = K, x = x, est = fJ.mean, lcl = fJ.lower, ucl = fJ.upper, true = fJ.true))
}
