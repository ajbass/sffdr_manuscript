##################################################
### Apply sfFDR to EGPA
##################################################
set.seed(12345)

library(sffdr)
library(tidyverse)

source("../00-helper.R")

# LD blocks
df <- readRDS("../ukbb_ldblocks_01.rds")
df$BP38 <- as.numeric(unlist(str_split(df$alt38, "[:_]", simplify = T))[,2])
df <- as_tibble(df) %>% dplyr::select(chr, BP38, group)
colnames(df) <- c("CHR38", "BP38", "LD")

primary <- c("EGPA")
for (trait in primary) {
  out_EGPA <- read.delim(paste0("./processed/", trait, ".tsv"), sep = "\t")

  # Remove any missing data
  out_EGPA <- out_EGPA %>% left_join(df)
  out_train <- out_EGPA %>%
    dplyr::select(SNPID, P, ASTAO, ASTCO, EOSC, LD, CHR38, BP38, INFO, ALT) %>%
    dplyr::filter(!is.na(ASTAO) & !is.na(ASTCO) & !is.na(EOSC))

  # Select representative SNPs randomly and based on informative variables
  out_tmp2 <- out_train %>%
    filter(is.na(LD)) %>%
    mutate(indep_snp_rand = FALSE,
           indep_snp_inform = FALSE)
  out_tmp <- out_train %>% filter(!is.na(LD)) %>%
    group_by(LD) %>%
    mutate(indep_snp_rand = sample(c(rep(FALSE, length(LD) - 1), TRUE)),
           indep_snp_inform = create_bool(cbind(ASTAO, ASTCO, EOSC)))
  out_train <- rbind(out_tmp, out_tmp2)

  z <- out_train[, 3:5]
  id <-  which(rowSums(is.na(as.matrix(z))) == 0)
  p <-  out_train$P[id]
  z <- z[id,]
  indep_snps_inform <- out_train$indep_snp_inform[id]
  indep_snps_rand <- out_train$indep_snp_rand[id]

  # apply sffdr w/ randomly selected SNP and informative-selected
  knots <-  c(0.005, 0.01, 0.025, 0.05, 0.1)
  out_sffdr_indp_r <- apply_sffdr(p, z,
                                  indep_snps_rand,
                                  epsilon = 1e-15,
                                  method =  "gam",
                                  lambda = seq(0.05, 0.9, 0.05),
                                  knots = knots)
  t1 <- proc.time()
  out_sffdr_indp_i <- apply_sffdr(p, z,
                                  indep_snps_inform,
                                  epsilon = 1e-15,
                                  method =  "gam",
                                  lambda = seq(0.05, 0.9, 0.05),
                                  knots = knots)
  t2 <- proc.time()[3] - t1[3]

  fp_i <- out_sffdr_indp_i$fpvalues
  out_train2 <- out_train[id,]

  # Save output
  out_train3 <- out_train2
  out_train3$fPi <-  fp_i
  out_train3$flfdri <-  out_sffdr_indp_i$flfdr
  out_train3$fQi <-  out_sffdr_indp_i$fqvalues
  out_train3$fPIi <-  out_sffdr_indp_i$fpi0
  out_train3$indep_snps <- indep_snps_inform
  out_train3$elapsed <- t2
  saveRDS(out_train3, file = paste0("./sffdr/", trait))
}
