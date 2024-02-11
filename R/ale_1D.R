
get_ale_post_1D <- function (X, model, pred_fun, J, K = 40, center = 'avg_pred') {

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
  if (center == 'avg_pred') fJ.mat.cen <- sweep(fJ.mat, 2, apply(fJ.mat, 2, \(x) sum((x[1:K] + x[2:(K + 1)]) / 2 * b1) / sum(b1)), '-')
  else if (center == 'median') {
    M <- max(which(stats::median(X[,J]) > z))
    fJ.mat.cen <- sweep(fJ.mat, 2, apply(fJ.mat, 2, \(x) (x[M] + x[M+1]) / 2), '-')
  }

  return (list(ale = fJ.mat.cen, x = z))
}

get_ale_true_1D <- function (X, f_true, J, K = 40, center = 'avg_pred') {

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
  y.hat1.true <- f_true(X1)
  y.hat2.true <- f_true(X2)

  # Compute local effects
  Delta.true <- y.hat2.true - y.hat1.true
  Delta.true <- tapply(Delta.true, a1, mean)

  # Accumulate local effects
  fJ.true <- c(0, cumsum(Delta.true))

  # Center ALE
  b1 <- as.numeric(table(a1))
  if (center == 'avg_pred') fJ.true <- fJ.true - sum((fJ.true[1:K] + fJ.true[2:(K + 1)])/2 * b1) / sum(b1)
  else if (center == 'median') {
    M <- max(which(stats::median(X[, J]) > z))
    fJ.true <- fJ.true - ((fJ.true[M] + fJ.true[M+1]) / 2)
  }

  return (as.numeric(fJ.true))
}
