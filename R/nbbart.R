#' @export
nbbart <- function(x, y, offset = rep.int(0, length(y)),
                   m = 200, k = 2, base = 0.95, power = 2, s_xi = 0.1,
                   num_iter = 5000, num_burn = 2500, num_thin = 5,
                   light_store = TRUE, seed = 2187){

  set.seed(seed)

  # Sample size
  N <- length(y)

  # Fit a fixed effects NB model to get initial estimates
  init.fit <- MASS::glm.nb(y ~ 1 + as.matrix(x) + offset(offset))

  # Initialize parameters
  xi <- init.fit$theta
  G  <- as.numeric(x %*% stats::coef(init.fit)[1:ncol(x)]) + stats::coef(init.fit)[1]

  # Linear predictor
  eta <- offset + G

  # Create dbarts sampler object
  control <- dbarts::dbartsControl(n.trees = m,
                                   n.burn = 0, n.samples = num_iter - num_burn, n.thin = 1,
                                   n.chains = 1L, keepTrees = FALSE, keepTrainingFits = TRUE,
                                   updateState = TRUE, verbose = FALSE)

  # Generate an example auxillary outcome to initialize sampler with
  omega <- pg::rpg_hybrid(y + xi, eta)[,1]
  z <- (y - xi) / (2 * omega)
  split.probs <- rep(1 / ncol(x), ncol(x))
  sampler <- dbarts::dbarts(x, z, offset = offset, weights = omega,
                            control = control,
                            resid.prior = fixed(1), sigma = 1,
                            tree.prior = cgm(power = power,
                                             base = base,
                                             split.probs = split.probs),
                            node.prior = normal(k))

  # Allocate posterior storage
  K <- (num_iter - num_burn) / num_thin
  k_keep <- seq(num_burn + 1, num_iter, by = num_thin)

  post <- list(G.mean = rep.int(NA_real_, K),
               xi = rep.int(NA_real_, K),
               var_counts = matrix(0, nrow = K, ncol = ncol(x),
                                   dimnames = list(NULL, colnames(x))),
               logmean = rep.int(NA_real_, K)
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
      sampler$setControl(control)
    }

    # Sample latent Polya-Gamma RV
    # Attempt hybrid sampler unless approximation fails, then use truncated sum of gammas
    if(length(utils::capture.output(omega <- pg::rpg_hybrid(y + xi, eta)[,1])) > 0){
      omega <- pg::rpg_gamma(y + xi, eta)[,1]
    }

    # Convert to Gaussian form
    z <- (y - xi) / (2 * omega)  # z ~ N(eta, diag(1 / omega))

    # Update BART
    sampler$setResponse(z)
    sampler$setWeights(omega)
    samples <- sampler$run(0L, 1L)
    G <- samples$train

    # Update linear predictor
    eta <- offset + G

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

      post$G.mean[Kk]       <- mean(G)
      post$xi[Kk]           <- xi
      post$var_counts[Kk,]  <- samples$varcount
      post$logmean[Kk]      <- log(xi) + post$G.mean[Kk]
    }

    pb$tick()
  }

  control@updateState <- TRUE
  sampler$setControl(control)

  sampler$storeState()
  post$bart <- sampler

  return(post)
}
