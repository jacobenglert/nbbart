# Program Name: spanbbart.R
# Author:       Jacob Englert
# Date:         31 October 2023
# Purpose:      Fit a Bayesian negative binomial regression with spatial random
#               effects under a proper CAR prior and nonparametric relative risk
#               mean function under a BART prior (using dbarts)

spanbbart <- function(x1, x2, y,
                      s = 1:length(y), t = rep(1, length(y)), # Data should be ordered by space then time
                      offset = numeric(length(y)), geo, seed = 2187,
                      m = 200, k = 2, base = 0.95, power = 2,
                      num_iter = 5000, num_burn = 2500, num_thin = 5,
                      light_store = TRUE,
                      b = rep(0, ncol(x1)), B = diag(ncol(x1)) * 1e4, # Fixed effect hyperparameters
                      c = 0.1, d = 0.1, # spatial random effect marginal variance (tau2) prior hyperparameters
                      s_xi = 0.1 # proposal distribution variance for xi (only used if xi_update_method = 'MH')
){

  set.seed(seed)

  # Sample size
  N <- length(y)
  nT <- length(unique(t))
  nS <- length(unique(s))

  # Spatial random effects design matrix
  s_unique <- match(unique(s), s)
  x_nu <- spam::spam(0, nrow = N, ncol = nS)
  x_nu[cbind(1:N, s)] <- 1

  # Compute spatial matrices
  W <- spdep::poly2nb(geo[s_unique]) |> spdep::nb2mat(style = 'B')
  D <- diag(rowSums(W))

  # Fit a non-spatial fixed effects NB model to get initial estimates
  init.fit <- MASS::glm.nb(y ~ 1 + as.matrix(x1) + offset(offset))

  # Initialize parameters
  beta  <- stats::coef(init.fit)[2:(ncol(x1) + 1)]
  nu    <- mvtnorm::rmvnorm(n = 1, sigma = diag(1, nS))[1,]
  tau2  <- 1
  rho   <- 0.5
  xi    <- init.fit$theta
  G     <- stats::coef(init.fit)[1]

  # Linear predictor
  fixeff <- as.numeric(x1 %*% beta)
  ranef <- as.numeric(x_nu %*% nu)
  eta <- offset + fixeff + G + ranef

  # Pre-calculate discrete prior distribution for rho
  lambda    <- eigen(solve(D) %*% W, only.values = TRUE)$values
  rho_vals  <- stats::qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
  rho_ll0   <- sapply(rho_vals, \(x) 0.5 * sum(log(1 - x * lambda)), simplify = TRUE)

  # Pre-calculate fixed effect precision matrix
  B_inv <- B; diag(B_inv) <- 1 / diag(B)

  # Create dbarts sampler object
  control <- dbarts::dbartsControl(n.trees = m,
                                   n.burn = 0,
                                   n.samples = num_iter - num_burn,
                                   n.thin = 1,
                                   n.chains = 1L,
                                   keepTrees = FALSE,
                                   keepTrainingFits = TRUE,
                                   updateState = TRUE,
                                   verbose = FALSE)

  # Generate an example auxillary outcome to initialize sampler with
  omega <- pg::rpg_hybrid(y + xi, eta)[,1]
  z <- (y - xi) / (2 * omega)
  r <- z - offset - fixeff - ranef
  split.probs <- rep(1 / ncol(x2), ncol(x2))
  sampler <- dbarts::dbarts(x2, r, weights = omega,
                            control = control,
                            resid.prior = fixed(1), sigma = 1,
                            tree.prior = cgm(power = power,
                                             base = base,
                                             split.probs = split.probs),
                            node.prior = normal(k))

  # Allocate posterior storage
  K <- (num_iter - num_burn) / num_thin
  k_keep <- seq(num_burn + 1, num_iter, by = num_thin)

  post <- list(alpha = rep.int(NA_real_, K),  # fixed intercept
               beta = matrix(NA_real_, K, ncol(x1)), # fixed effects
               nu = matrix(NA_real_, K, nS),   # spatial random effects
               xi = rep.int(NA_real_, K),     # dispersion
               rho = rep.int(NA_real_, K),    # spatial autocorrelation
               tau2 = rep.int(NA_real_, K),   # spatial random effects variance
               var_counts = matrix(0, nrow = K, ncol = ncol(x2),
                                   dimnames = list(NULL, colnames(x2))),
               logmean = rep.int(NA_real_, K),
               sigma = rep.int(NA_real_, K)
  )

  # If light storage is not selected, then return linear and BART predictors
  # These matrices are K x N, so require a lot of memory and may crash R
  if(!light_store){
    post$eta <- matrix(NA_real_, K, N)
    post$G <- matrix(NA_real_, K, N)
  }


  # MCMC
  xi_acc <- 0
  pb <- progress::progress_bar$new(
    format = "[:bar] Iteration :current/:total. Total time elapsed: :elapsedfull",
    total = num_iter, clear = FALSE, width = 100)
  for(k in seq_len(num_iter)){

    if(k > num_burn){
      control@keepTrees <- TRUE
      # control@updateState <- FALSE
      sampler$setControl(control)
    }

    # Sample latent Polya-Gamma RV
    # Attempt hybrid sampler unless approximation fails, then use truncated sum of gammas
    if(length(utils::capture.output(omega <- pg::rpg_hybrid(y + xi, eta)[,1])) > 0){
      omega <- pg::rpg_gamma(y + xi, eta)[,1]
    }

    # Convert to Gaussian form
    z <- (y - xi) / (2 * omega)  # z ~ N(eta, diag(1 / omega))

    # Update spatial weight matrices
    D_rho_W <- spam::as.spam(D - rho * W)

    # Update spatial random effects
    r <- z - offset - fixeff - G
    nu_Sigma <- spam::solve(spam::crossprod(x_nu * sqrt(omega)) + (1 / tau2) * (D_rho_W))
    nu_mu <- nu_Sigma %*% spam::crossprod(x_nu, omega * r)
    nu <- mvtnorm::rmvnorm(n = 1, mean = nu_mu, sigma = nu_Sigma)[1,]
    nu <- nu - mean(nu)
    ranef <- as.numeric(x_nu %*% nu)

    # Update spatial random effects variance
    tau2 <- 1 / stats::rgamma(1, c + nS / 2, d + (nu %*% (D_rho_W) %*% nu) / 2)

    # Update spatial random effects correlation
    rho_ll <- rho_ll0 + rho_vals / (2 * tau2) * as.numeric(nu %*% W %*% nu)
    rho <- sample(rho_vals, size = 1, prob = exp(rho_ll - max(rho_ll)))

    # Update fixed effects
    r <- z - offset - G - ranef
    beta_Sigma <- solve(B_inv + crossprod(x1 * sqrt(omega)))
    beta_mu <- beta_Sigma %*% (B_inv %*% b + crossprod(x1, omega * r))
    beta <- mvtnorm::rmvnorm(n = 1, mean = beta_mu, sigma = beta_Sigma)[1,]
    fixeff <- x1 %*% beta

    # Update BART
    r <- z - offset - fixeff - ranef
    sampler$setResponse(r)
    sampler$setWeights(omega)
    samples <- sampler$run(0L, 1L)
    G <- samples$train

    # Update linear predictor
    eta <- offset + fixeff + G + ranef

    # Update dispersion parameter
    # (Metropolis-Hastings proposed from centered truncated normal)
    q <- 1 / (1 + exp(eta)) # 1 - Pr(success)
    xi_prop <- msm::rtnorm(1, mean = xi, sd = s_xi, lower = 0)
    r_xi <- sum(stats::dnbinom(y, size = xi_prop, prob = q, log = TRUE)) -
      sum(stats::dnbinom(y, size = xi, prob = q, log = TRUE)) +
      msm::dtnorm(xi, mean = xi_prop, sd = s_xi, lower = 0, log = TRUE) -
      msm::dtnorm(xi_prop, mean = xi, sd = s_xi, lower = 0, log = TRUE)

    if(log(stats::runif(1)) < r_xi){
      xi <- xi_prop
      xi_acc <- xi_acc + 1
    }

    # Update tuning parameter
    if(k <= num_burn & k %% 100 == 0){
      xi_acc_rate <- xi_acc / k
      if(xi_acc_rate > 0.6) s_xi <- 1.1 * s_xi
      if(xi_acc_rate < 0.2) s_xi <- 0.8 * s_xi
      cat(s_xi)
    }
    #invisible(sampler$state)
    if(k %in% k_keep){

      Kk <- which(k_keep == k) # ID posterior sample index

      if(!light_store){
        post$eta[Kk,]   <- eta
        post$G[Kk,]     <- G
      }

      post$alpha[Kk]        <- mean(G)
      post$beta[Kk,]        <- beta
      post$nu[Kk,]          <- nu
      post$xi[Kk]           <- xi
      post$rho[Kk]          <- rho
      post$tau2[Kk]         <- tau2
      post$var_counts[Kk,]  <- samples$varcount
      post$sigma[Kk]        <- samples$sigma
      post$logmean[Kk]      <- log(xi) + post$alpha[Kk]
    }

    pb$tick()
  }

  control@updateState <- TRUE
  sampler$setControl(control)

  sampler$storeState()
  post$bart <- sampler

  return(post)
}
