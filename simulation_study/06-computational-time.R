source("./00-generate_data.R")
source("../00-helper.R")

summarise_study <- function(qvalues, oracle, fdr = seq(0.001, 0.01, 0.001)) {
  return(sumfunc(q = qvalues, oracle = oracle, fdr = fdr))
}

run_t1e <- function(m,
                    num.z,
                    num.null.z,
                    prior.strength,
                    signal.density,
                    signal.density.z,
                    transformation,
                    fmodel,
                    prior.coverage,
                    pi0.method,
                    seed) {
  # generate informative studies
  ivars <- generate_z(m,
                      num.z,
                      num.null.z,
                      signal.density.z = signal.density.z,
                      prior.coverage = prior.coverage)
  z <- ivars$z
  
  # generate primary study p-values
  study <- generate_p(z = z[, 1:num.z, drop = F],
                      pi0z = 1 - ivars$pi0,
                      status = ivars$status,
                      w =  rep(1/num.z, num.z),
                      prior.strength = prior.strength,
                      signal.density = signal.density)
  
  # rank transform
  z <- apply(ivars$z, 2, FUN = function(x) rank(x)/length(x))
  colnames(z) <- paste0("z", 1:ncol(z))
  
  # apply sffdr
  create_model <- sffdr:::pi0_model(z, knots = c(0.005, 0.01, 0.025, 0.05, 0.1))
  t1 <- proc.time()[3]
  fpi0 <- sffdr:::fpi0est(p = study$p,
                          z = create_model$zt,
                          pi0_model = create_model$fmod,
                          lambda = seq(0.05, 0.9, 0.05),
                          constrained.p = TRUE)$fpi0
  
  tmp <- sffdr:::sffdr(study$p,
                       surrogate = fpi0,
                       fpi0 = fpi0)
  t2 <- proc.time()[3] - t1

  df <- data.frame(type = "sfFDR",
                   time = t2) 
  return(df)
}

library(sffdr)
library(tidyverse)
library(qvalue)
library(digest)
library(gam)
library(splines)

# simulation design
design <- expand.grid(fmodel = "both",
                      m = c(1e5, 2.5e5, 5e5, 7.5e5, 1e6),
                      rep = 1:10,
                      pi0.method = "gam",
                      transformation = "probit",
                      prior.strength = c("Large"),
                      signal.density = c("High"),
                      signal.density.z = c("High"),
                      prior.coverage = 0.025,
                      num.z = 1:5,
                      num.null.z = 0)
 
design <- design %>%
  group_by(num.z, m, fmodel, transformation, prior.coverage, num.null.z, signal.density.z, prior.strength, pi0.method, signal.density, rep) %>%
  mutate(seed = readBin(digest(c(num.z, m, transformation, fmodel, prior.coverage, signal.density.z, pi0.method, num.null.z, prior.strength, signal.density, rep), raw = TRUE), "integer"))

out <- NULL
for (ii in 1:nrow(design)) {
  print(ii)
  set.seed(design[ii,]$seed)
  df0 <- run_t1e(design[ii,]$m,
                 design[ii,]$num.z,
                 design[ii,]$num.null.z,
                 design[ii,]$prior.strength,
                 design[ii,]$signal.density,
                 design[ii,]$signal.density.z,
                 design[ii,]$transformation,
                 design[ii,]$fmodel,
                 design[ii,]$prior.coverage,
                 design[ii,]$pi0.method,
                 design[ii,]$seed)
  df2 <- cbind(design[ii,], df0)
  out <- dplyr::bind_rows(df2, out)
}

saveRDS(out, paste0("./data-fdr/06-time.rds"))
