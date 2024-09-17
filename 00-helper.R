# Wraper of sffdr function
#
# @param p: p-values of primary study
# @param z: summary statistics of informative studies
# @param indep_snps: vector of LD-independent SNPs (boolean)
# @param knots: location of knots for functional pi0
# @param epsilon: truncation for sffdr density estimation
# @param nn: nearest neighborhood paramter for density estimation
# @param lambda: lambda values for pi0 estimation
# @param method: pi0 estimation GAM model only
#
# @return functional FDR quantities
apply_sffdr <- function(p,
                        z,
                        indep_snps,
                        knots = c(0.005, 0.0025, 0.01, 0.05, 0.1),
                        epsilon = 1e-15,
                        nn = NULL,
                        lambda = seq(0.05, 0.9, 0.05),
                        method = "gam") {
  # initialization
  fqvalues <- fpi0 <- flfdr <- fpvalues <- rep(NA, length(p))

  # run sffdr
  create_model <- sffdr:::pi0_model(as.matrix(z),
                                    knots = knots)

  fpi0 <- sffdr:::fpi0est(p = as.matrix(p),
                          z = create_model$zt,
                          indep_snps = indep_snps,
                          pi0_model = create_model$fmod,
                          lambda = lambda,
                          method = method)

  fpi0 <- fpi0$fpi0
  sffdr.out <- sffdr:::sffdr(p,
                             fpi0 = fpi0,
                             indep_snps = indep_snps,
                             epsilon = epsilon,
                             nn = nn)

  fqvalues  <- sffdr.out$fqvalues
  fpvalues  <- sffdr.out$fpvalues
  flfdr  <- sffdr.out$flfdr
  fpi0  <- fpi0

  return(list(fpvalues = fpvalues,
              fqvalues = fqvalues,
              flfdr = flfdr,
              fpi0 = fpi0))
}

type1error <- function(p,
                       oracle,
                       threshold = c(5e-8, 5e-7, 5e-6, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1 )) {
  out <- NULL
  for (i in 1:length(threshold)) {
    alpha <- threshold[i]
    sig <- p <= alpha
    tot <- sum(sig)
    fpr <- mean(sig[!oracle])
    tpr <- mean(sig[oracle])
    tmp <- data.frame(threshold = alpha,
                      discoveries = tot,
                      emptdr = tpr,
                      empfdr = fpr)
    out <- rbind(tmp, out)
  }
  return(out)
}

create_bool <- function(z, thres = 1e-3) {
  boo <- rep(FALSE, nrow(z))
  zmin <- apply(z, 1, min)
  if (min(zmin) < thres) {
    id <- which.min(zmin)
  } else {
    id <- sample(1:nrow(z), size = 1)
  }
  boo[id] <- TRUE
  return(boo)
}

sumfunc <- function(q, oracle, fdr) {
  out <- NULL
  for (i in 1:length(fdr)) {
    fdr.level <- fdr[i]
    sig <- q <= fdr.level
    fp <- sum(sig[!oracle])
    tp <- sum(sig[oracle])
    tmp <- data.frame(threshold = fdr.level,
                      discoveries = fp + tp,
                      emptdr = tp /  sum(oracle),
                      empfdr = fp / max(1, (fp + tp)))
    out <- rbind(tmp, out)
  }
  return(out)
}
