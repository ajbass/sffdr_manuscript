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
  
  # generate primary study
  study <- generate_p(z = z[, 1:num.z, drop = F],
                      pi0z = 1 - ivars$pi0,
                      status = ivars$status,
                      w =  weights(a = ivars$alpha[1:num.z], b = 5),
                      prior.strength = prior.strength,
                      signal.density = signal.density,
                      fmodel = fmodel)

  # rank transform summary statistics
  z <- apply(ivars$z, 2, FUN = function(x) rank(x)/length(x))
  colnames(z) <- paste0("z", 1:ncol(z))
  
  # LD blocks: assign cluster sizes to independent SNPs 
  cluster_size <- readRDS("./main-sffdr/LDblocks/ukbb_set_sizes.rds")
  cluster_size <- sample(cluster_size, replace = T, size = length(study$p))
  zp <- ivars$p[rep(1:nrow(ivars$p), times = cluster_size),]
  zp <- apply(zp, 2, FUN = function(x) rank(x) / length(x))
  np <- rep(study$p, times = cluster_size)
  id <- cumsum(cluster_size)
  indp_status <- rep(FALSE, length(np))
  indp_status[id] <- TRUE

  neworacle <- rep(study$oracle,
                   times = cluster_size)

  # apply sffdr
  create_model <- sffdr:::fmodel(zp, knots = c(0.005, 0.01, 0.025, 0.05, 0.1))

  fpi0 <- sffdr:::efpi0(p = np,
                        dt = create_model$zt,
                        fm = create_model$fmod,
                        indep_snps = indp_status,
                        method = pi0.method,
                        lambda = seq(0.05, 0.9, 0.05))$table$fpi0

  tmp <- sffdr:::sffdr(np,
                       surrogate = fpi0,
                       monotone.window = NULL,
                       transformation = transformation,
                       fpi0 = fpi0)
  ffdr.out <- tmp$fqvalues
  lfdr <- tmp$flfdr
  fp <- tmp$fpvalues

  # Summarize - Functional P-values
  oracle_t1e <- type1error(rep(study$fp,
                               times = cluster_size), neworacle)
  sffdr_t1e <- type1error(fp, neworacle)
  raw_t1e <- type1error(np, neworacle)
  sffdr_t1e$method = "sffdr"
  oracle_t1e$method = "oracle"
  raw_t1e$method = "raw"
  dt2 <- rbind(raw_t1e,
               oracle_t1e,
               sffdr_t1e)
  dt2$quantity <- "fpvalues"

  # Summarize - Q-values
  oracle_q <- summarise_study(rep(study$oracle_q,
                                  times = cluster_size), neworacle)
  sffdr_q <- summarise_study(ffdr.out, neworacle)
  raw_q <- summarise_study(qvalue(np)$qv, neworacle)
  sffdr_q$method = "sffdr"
  oracle_q$method = "oracle"
  raw_q$method = "raw"
  dt3 <- rbind(raw_q,
               oracle_q,
               sffdr_q)
  dt3$quantity <- "fqvalues"

  dt3 <- rbind(dt2, dt3)

  # Summarize - Pi0-values
  oracle_pi0 <- mean(!neworacle)
  sffdr_pi0 <- mean(tmp$pi0)
  raw_pi0 <- qvalue(np)$pi0
  dt4 <- data.frame(method = c("sffdr", "oracle", "raw"),
                    pi0 = c(sffdr_pi0, oracle_pi0, raw_pi0))
  dt4 <- dt3 %>% left_join(dt4)
  return(dt4)
}

library(adaptMT)
library(sffdr)
library(tidyverse)
library(swfdr)
library(qvalue)
library(digest)
library(locfit)
library(gam)
library(splines)
library(CAMT)

# study design
design <- expand.grid(fmodel = "both",
                      m = 150000,
                      rep = 1:100,
                      pi0.method = "gam",
                      transformation = "probit",
                      prior.strength = c("Large", "Moderate" ),
                      signal.density = c("High", "Medium", "Low"),
                      signal.density.z = c("High", "Medium", "Low"),
                      prior.coverage = 0.025,
                      num.z = 3,
                      num.null.z = 0)

# null setting (prior.strength = "None")
design2 <- expand.grid(fmodel = "both",
                       m = 150000,
                       rep = 1:100,
                       pi0.method = "gam",
                       transformation = "probit",
                       prior.strength = "None",
                       signal.density = c("High", "Medium", "Low"),
                       signal.density.z = "High",
                       prior.coverage = 0.025,
                       num.z = 3,
                       num.null.z = 0)

design <- rbind(design, design2)

design <- design %>%
  group_by(num.z, m, fmodel, transformation, prior.coverage, num.null.z, signal.density.z, prior.strength, pi0.method, signal.density, rep) %>%
  mutate(seed = readBin(digest(c(num.z, m, transformation, fmodel, prior.coverage, signal.density.z, pi0.method, num.null.z, prior.strength, signal.density, rep), raw=TRUE), "integer"))

args <- commandArgs(trailingOnly = TRUE)
id2 <- args[1]
id <- (1:nrow(design) %% 300) == id2
design <- design[id,]

out <- out2 <- NULL
for (ii in 1:nrow(design)) {
  print(ii)
  set.seed(design[ii,]$seed)
  df0 <- tryCatch(run_t1e(design[ii,]$m,
                          design[ii,]$num.z,
                          design[ii,]$num.null.z,
                          design[ii,]$prior.strength,
                          design[ii,]$signal.density,
                          design[ii,]$signal.density.z,
                          design[ii,]$transformation,
                          design[ii,]$fmodel,
                          design[ii,]$prior.coverage,
                          design[ii,]$pi0.method,
                          design[ii,]$seed),
                  error = function(x) data.frame(alpha = NA, empfdr = NA, total = NA, method = NA))

  df2 <- cbind(design[ii,], df0)
  out2 <- dplyr::bind_rows(df2, out2)
}

saveRDS(out2, paste0("./data/04-simulations-t1e-LD-", id2, ".rds"))
