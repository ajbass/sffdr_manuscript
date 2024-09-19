# Generate informative summary statistic data
#
# @param m: total number of tests
# @param num.z: number of informative summary statistics
# @param num.null.z: number of null informative summary statistics
# @param prior.coverage: proportion of non-null tests upper bound
# @param signal.density.z: signal density of the informative statistics
#
# @return informative summary statistics
generate_z <- function(m,
                       num.z,
                       num.null.z,
                       prior.coverage,
                       signal.density.z) {
  # initializations
  p <- matrix(nrow = m, ncol = num.z + num.null.z)
  status <- matrix(0, nrow = m, ncol = num.z + num.null.z)
  mu <- rep(0, num.z + num.null.z)
  Sigma <- matrix(0,
                  nrow = num.z + num.null.z,
                  ncol = num.z + num.null.z)
  diag(Sigma) <- 1
  B <- 5

  if (signal.density.z == "High") {
    A <- rep(0.2, num.z + num.null.z)
  } else if (signal.density.z == "Medium") {
    A <- rep(0.3, num.z + num.null.z)
  } else if (signal.density.z == "Low") {
    A <- rep(0.4, num.z + num.null.z)
  }

  # Generate proportion of non-null tests for each informative study
  prop_null <- 1 - runif(n = num.z + num.null.z,
                         min = prior.coverage / 2,
                         max = prior.coverage)

  # generate p-values
  for (i in 1:(num.z + num.null.z)) {
    num.null <- round(m * prop_null[i])
    num.alt <- m - num.null
    p.null <- runif(n = num.null)
    p.alt <- rbeta(num.alt, shape1 = A[i], shape2 = B)
    p[1:num.alt, i] <- p.alt
    p[(num.alt+1):m, i] <- p.null
    status[1:num.alt,i] <- 1
  }

  # randomly permute null statistics so it is unassociated with primary
  if (num.null.z != 0) {
    p[, (num.z + 1):(num.z + num.null.z)] <- p[sample(1:nrow(p), replace = F),  (num.z + 1):(num.z + num.null.z)]
  }

  # rank transform
  z <- apply(p, 2, FUN = function(x) rank(x) / length(x))
  return(list(z = z,
              p = p,
              pi0 = prop_null,
              alpha = A,
              status = status))
}

# Simulate a functional pi0 relationship using the informative variables
#
# @return functional pi0
fpi0 <- function(z,
                 pi0z,
                 w,
                 status,
                 prior.strength) {
  # initializations
  if (prior.strength == "None") {
    pi0 <- matrix(0.98, nrow = nrow(z), ncol = 1)
    return(pi0)
  }

  tpi0 <- matrix(0.98, nrow = nrow(z), ncol = ncol(z))
  if (prior.strength == "Large") {
    beta <- 0.6
  } else {
    beta <- 0.3
  }

  # generate functional pi0
  for (i in 1:ncol(z)) {
    tpi0[z[, i] < pi0z[i] & status[, i] == 1, i] <- 0.98 * (z[z[, i] < pi0z[i] & status[, i] == 1, i] / pi0z[i]) ^ beta
  }
  pi0 <- rowSums(sapply(1:ncol(z), FUN = function(i) w[i] * tpi0[, i]))

  pi0[pi0 > 1] <- 1
  pi0[pi0 < 0] <- 1e-3
  return(pi0)
}

fden <- function(z,
                 pi0,
                 w,
                 status,
                 prior.strength = "Large",
                 signal.density = "High") {
  # Simulate a functional density that depends on the informative variable, z.
  # Output:
  #   functional density
  beta <- 5
  if (signal.density == "High") {
    alpha <- 0.3
  } else if (signal.density == "Medium") {
    alpha <- 0.4
  } else {
    alpha <- 0.5
  }

  ww <- matrix(0, nrow = nrow(z), ncol = ncol(z))
  for (i in 1:ncol(z)) {
    ww[z[, i] < pi0[i] & status[, i] == 1, i] <- (pi0[i] - z[z[, i] < pi0[i] & status[,i] == 1, i])  / pi0[i]
  }
  altw <- rowSums(sapply(1:ncol(z), FUN = function(i) w[i] * (1 - ww[,i])))
  if (prior.strength == "Large") {
    c <- alpha / 2
  } else if (prior.strength == "Moderate") {
    c <- alpha / 4
  } else {
    c <- 0
  }

  alpha <- alpha  - c * (1 - altw)
  p <- rbeta(n = nrow(z), shape1 = alpha, shape2 = beta)

  df <- list(p = p,
             beta = beta,
             alpha = alpha)
  return(df)
}

# Generate p-values for primary study
#
# @param z: the informative summary statistics
# @param pi0z: the pi0 of the informative studies
# @param status_z: the true status of each test
# @param w: weights
# @param prior.strength: the prior strength of informative studies
# @param signal.density: signal density of the informative studies
#
# @return data frame of p-values and oracle values for simulation
generate_p <- function(z,
                       pi0z,
                       status_z,
                       w,
                       prior.strength,
                       signal.density) {
  # initialization
  m <- nrow(z)
  p <- vector(length = m)

  # generate functional proportion of null tests
  pi0 <- fpi0(z = z,
              pi0z = pi0z,
              w = w,
              prior.strength = prior.strength,
              status = status_z)

  # generate alternative density of p-values
  alt <- fden(z = as.matrix(z),
              w = w,
              pi0 = pi0z,
              status = status_z,
              prior.strength = prior.strength,
              signal.density = signal.density)

  # generate true status of primary study p-values
  status <- (rbinom(n = m, size = 1, prob = 1 - pi0) == 1)

  # null + alternative p-values
  p[!status] <- runif(n = sum(!status))
  p[status] <- alt$p[status]

  # calculate oracle values
  fden <- (pi0 + (1 - pi0) * dbeta(p, shape1 = alt$alpha, shape2 = alt$beta))
  lfdr <- pi0 / fden
  fp <- sffdr:::fpvalues_raw(lfdr)
  out <- sort(lfdr, index.return = TRUE)
  out$x[out$x > 1] <- 1
  fdr <- cumsum(out$x) / (1:length(lfdr))
  fdr <- fdr[order(out$ix)]

  return(data.frame(p = p,
                    fp = fp,
                    oracle_lfdr = lfdr,
                    oracle_q = fdr,
                    oracle_pi0 = pi0,
                    oracle_fden = fden,
                    oracle = status))
}

weights <- function(a, b) {
  w <- (1 - a) / (b - a)
  w <- w / sum(w)
  return(w)
}
