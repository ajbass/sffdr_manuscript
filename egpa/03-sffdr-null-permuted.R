##################################################
### Apply sfFDR to EGPA with null traits
##################################################
source("../00-helper.R")
set.seed(12345)

library(sffdr)
library(tidyverse)

# Gather null summary statistics
files <- list.files("../ukbb/data/assoc_null/", full.names = T)
cdf <- NULL
prop <- seq(0.1, 1, 0.1)
for (f in files) {
  out <- read_tsv(f)
  str <-  str_split(str_split(f, pattern = "subset")[[1]][2], "\\.")[[1]]
  i <- str[1]
  trait <- str[2]
  out$trait <- trait
  out$replicate <- as.numeric(i)
  out$downsample <- 1
  cdf <- dplyr::bind_rows(out, cdf)
}

cdf2 <- cdf %>% filter(downsample == 1)
bfp <-  cdf2 %>%
  filter(trait == "bfp")
cholesterol <-  cdf2 %>%
  filter(trait == "cholesterol")
triglycerides <-  cdf2 %>%
  filter(trait== "triglycerides")

# LD regions
df <- readRDS("../ukbb_ldblocks_01.rds")
df$BP38 <- as.numeric(unlist(str_split(df$alt38, "[:_]", simplify = T))[,2])
df <- as_tibble(df) %>% dplyr::select(chr, BP38, group)
colnames(df) <- c("CHR38", "BP38", "LD")

# Combine data
df_cond <- cbind(bfp[, c(1:8, 14)],
            data.frame(bfp = bfp$P,
                       cholesterol = cholesterol$P,
                       triglycerides = triglycerides$P))
df_cond <- df_cond %>% dplyr::rename(CHR38 ="#CHROM", BP = POS) %>% select(replicate, BP, CHR38, bfp, cholesterol, triglycerides)

# load EGPA data
out_raw <- read.delim(paste0("./primary/EGPA_Lyons_31719529_1-hg38.tsv.gz"), sep = "\t")
map <- out_raw %>% select(SNPID, BP38, BP)

primary <- c("EGPA")
ret_df <- NULL
for (rep in 1:10) {
  dfc <- df_cond %>% filter(replicate == rep)
  for (trait in primary) {
    # remove missing data and align UK biobank data w/ EGPA
    out_EGPA <- read.delim(paste0("./processed/", trait, ".tsv"), sep = "\t")
    out_EGPA <- out_EGPA %>%
      left_join(df) %>%
      left_join(map)

    out_EGPA2 <- dfc %>%
      left_join(out_EGPA, by = c("CHR38", "BP"))
    out_train <- out_EGPA2 %>%
      dplyr::filter(!is.na(P)) %>%
      dplyr::select(SNPID, P, bfp, cholesterol,triglycerides, LD) %>%
      dplyr::filter(!is.na(bfp) & !is.na(cholesterol) & !is.na(triglycerides))

    # select representative SNPs to train model
    out_tmp2 <- out_train %>%
      filter(is.na(LD)) %>%
      mutate(indep_snp_rand = FALSE,
             indep_snp_inform = FALSE)
    out_tmp <- out_train %>% filter(!is.na(LD)) %>%
      group_by(LD) %>%
      mutate(indep_snp_rand = sample(c(rep(FALSE, length(LD) - 1), TRUE)),
             indep_snp_inform = create_bool(cbind(bfp, cholesterol, triglycerides)))
    out_train <- rbind(out_tmp, out_tmp2)

    z <- out_train[, 3:5]
    id <-  which(rowSums(is.na(as.matrix(z))) == 0)
    p <-  out_train$P[id]
    z <- z[id,]
    indep_snps_inform =  out_train$indep_snp_inform[id]

    # apply sffdr
    knots <- c(0.005, 0.01, 0.025, 0.05, 0.1)
    out_sffdr_indp_i <- apply_sffdr(p, z,
                                    indep_snps_inform,
                                    epsilon = 1e-15,
                                    method =  "gam",
                                    lambda = seq(0.05, 0.9, 0.05),
                                    knots = knots)

    fp_i <- out_sffdr_indp_i$fpvalues
    ret <- data.frame(id = 1:length(p),
                      rep = rep,
                      p = p,
                      fp  = fp_i,
                      fq  = out_sffdr_indp_i$fqvalues,
                      fpi0 = out_sffdr_indp_i$fpi0,
                      indep_snps = indep_snps_inform,
                      trait = trait)
    ret_df <- rbind(ret, ret_df)
  }
}

saveRDS(ret_df, file = "./functional-pvalue-null-inform.rds")
