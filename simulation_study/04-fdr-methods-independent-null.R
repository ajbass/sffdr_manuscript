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
  study$p <- runif(m)
  study$fp <- study$p
  study$oracle_lfdr <- rep(1, m)
  study$oracle_q <- rep(1, m)
  study$oracle_pi0 <- rep(1, m)
  study$oracle_fden <- rep(1, m)
  study$oracle <- rep(0, m)
  # rank transform
  z <- apply(ivars$z, 2, FUN = function(x) rank(x)/length(x))
  colnames(z) <- paste0("z", 1:ncol(z))
  
  create_model <- sffdr:::pi0_model(z, knots = c(0.005, 0.01, 0.025, 0.05, 0.1))
  
  # apply CAMT
  camt.obj.fdr <- camt.fdr(pvals = study$p,
                           pi0.var = create_model$fmod,
                           f1.var = create_model$fmod,
                           data = create_model$zt)
  
  # apply AdaPT
  adapt.fdr <- adapt_gam(pvals = study$p,
                         x = create_model$zt,
                         pi_formulas = as.character(create_model$fmod)[2],
                         mu_formulas = as.character(create_model$fmod)[2],
                         alphas = c(0.001, 0.005, 0.01))
  
  # apply Boca-Leek
  design <- model.matrix(as.formula(create_model$fmod), data = create_model$zt)
  bocaleek.out <- swfdr::lm_qvalue(study$p, X = design)
  
  # Summarize - Functional P-values
  oracle_t1e <- type1error(study$fp, study$oracle)
  raw_t1e <- type1error(study$p, study$oracle)
  oracle_t1e$method = "oracle"
  raw_t1e$method = "raw"
  dt2 <- rbind(raw_t1e,
               oracle_t1e)
  dt2$quantity <- "fpvalues"
  
  # Summarize - Q-values
  oracle_q <- summarise_study(study$oracle_q, study$oracle)
  adapt.sum <- summarise_study(adapt.fdr$q, study$oracle)
  camt.sum <- summarise_study(camt.obj.fdr$fdr, study$oracle)
  bl.sum <- summarise_study(bocaleek.out$qvalues, study$oracle)
  raw_q <- summarise_study(qvalue(study$p)$qv, study$oracle)
  oracle_q$method <- "oracle"
  camt.sum$method <- "CAMT"
  adapt.sum$method <- "adapt"
  bl.sum$method <- "Boca-Leek"
  raw_q$method <- "raw"
  dt3 <- rbind(raw_q,
               oracle_q,
               adapt.sum,
               camt.sum,
               bl.sum)
  dt3$quantity <- "fqvalues"
  dt3 <- rbind(dt2, dt3)
  
  lid <- length(adapt.fdr$params)
  raw_pi0 <- qvalue(study$p)$pi0
  
  # Density
  k <- camt.obj.fdr$k
  f1 <- (1 - k) * study$p^(-k)
  camt.f <- (1 - camt.obj.fdr$pi0) * f1 + camt.obj.fdr$pi0
  pix <- adapt.fdr$params[[lid]]$pix
  mux <- adapt.fdr$params[[lid]]$mux
  adapt.f <- (pix *  adapt.fdr$dist$h(study$p, mux) + 1 - pix)
  # Output pi0 values
  adapt <- 1 - adapt.fdr$params[[lid]]$pix
  df <- data.frame(method = c("adapt", "CAMT", "Boca-Leek"),
                   mse_pi0 = c(mean((adapt - study$oracle_pi0)^2), 
                               mean((camt.obj.fdr$pi0 - study$oracle_pi0)^2),
                               mean((bocaleek.out$pi0 - study$oracle_pi0)^2)),
                   mse_den = c(mean((-log(adapt.f) + log(study$oracle_fden))^2),
                               mean((-log(camt.f) + log(study$oracle_fden))^2),
                               NA)) 
  
  # global pi0 values
  oracle_pi0 <- mean(study$oracle_pi0)
  adapt_pi0 <- 1 - mean(adapt.fdr$params[[lid]]$pix)
  camt_pi0 <- mean(camt.obj.fdr$pi0)
  bl_pi0 <- mean(bocaleek.out$pi0)
  dt4 <- data.frame(method = c("oracle", "raw", "adapt", "CAMT", "Boca-Leek"),
                    pi0 = c(oracle_pi0, raw_pi0, adapt_pi0, camt_pi0, bl_pi0))
  dt4 <- dt3 %>% left_join(dt4) %>% left_join(df)
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

# simulation design (null setting: prior.strength = "None")
design <- expand.grid(fmodel = "both",
                       m = 150000,
                       rep = 1:500,
                       pi0.method = "gam",
                       transformation = "probit",
                       prior.strength = "None",
                       signal.density = c("High"),
                       signal.density.z = c("High", "Medium", "Low"),
                       prior.coverage = 0.025,
                       num.z = 3,
                       num.null.z = 0)

design <- design

design <- design %>%
  group_by(num.z, m, fmodel, transformation, prior.coverage, num.null.z, signal.density.z, prior.strength, pi0.method, signal.density, rep) %>%
  mutate(seed = readBin(digest(c(num.z, m, transformation, fmodel, prior.coverage, signal.density.z, pi0.method, num.null.z, prior.strength, signal.density, rep), raw = TRUE), "integer"))

out <- NULL
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
  out <- dplyr::bind_rows(df2, out)
}

saveRDS(out, paste0("./data-fdr/04-fdr-methods-independent-null.rds"))
