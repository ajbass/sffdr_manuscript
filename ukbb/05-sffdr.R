##################################################
### Applying sfFDR to UK Biobank data ###
##################################################
source("../00-helper.R")

library(tidyverse)
library(locfit)
library(qvalue)
library(splines)
library(fFDR)
library(sffdr)
library(tidyverse)

# load LD blocks
df <- readRDS("../ukbb_ldblocks_01.rds")
df$BP38 <- as.numeric(unlist(str_split(df$alt38, "[:_]", simplify = T))[,2])
df_LD <- as_tibble(df) %>% dplyr::select(rsid, BP38, group)
colnames(df_LD) <- c("ID", "POS", "LD")

# Load data: p-values from PLINK
pdf <- readRDS("./pcombined_pvals.rds")
cdf <- readRDS("./ccombined_pvals.rds")
fcdf <- readRDS("./cfcombined_pvals.rds")
pdf <- pdf %>% select(trait, downsample, `#CHROM`,
                      OBS_CT, POS, ID, indep_snps,
                      BETA, SE,  T_STAT, P, type)
cdf <- cdf %>% select(trait, downsample, `#CHROM`,
                      OBS_CT, POS, ID, indep_snps,
                      BETA, SE,  T_STAT,P, type)
fcdf <- fcdf %>% select(trait, downsample, `#CHROM`,
                        OBS_CT, POS, ID, indep_snps,
                        BETA, SE,  T_STAT,P, type)

# Reorganize data
fcdf_bmi <- fcdf %>% filter(trait == "bmi", downsample == 1)
cdf_bmi <- cdf %>% filter(trait == "bmi")
bfp <- fcdf %>%
  filter(trait == "bfp", downsample == 1)
cholesterol <- fcdf %>%
  filter(trait == "cholesterol", downsample == 1)
triglycerides <- fcdf %>%
  filter(trait== "triglycerides", downsample == 1)
df <- cbind(bfp[, 1:6],
            data.frame(bfp = bfp$P,
                       cholesterol = cholesterol$P,
                       triglycerides = triglycerides$P))
fcdf_cor <- df %>%
  filter(!is.na(bfp), !is.na(cholesterol), !is.na(triglycerides)) %>% distinct()

# Go through each downsample proportion and apply sffdr
prop <- seq(0.1, 1, 0.1)

for (i in 1:10) {
  set.seed(i)
  # randomly select representative SNP within each block
  fcdf_tmp <- fcdf_cor %>%
    left_join(df_LD, by = c("ID")) %>% filter(!is.na(LD))%>%
    group_by(LD) %>%
    mutate(indep_snp_rand = sample(c(rep(FALSE, length(LD) - 1), TRUE)))
  tmp_pdf <- pdf %>%
    filter(trait == "bmi", downsample == prop[i])

  fcdf_tmp <- fcdf_tmp %>%
    dplyr::rename(POS = POS.x) %>%
    left_join(tmp_pdf, by = c("POS", "ID"))

  ff <- is.na(fcdf_tmp$P) | is.na(fcdf_tmp$bfp)
  fcdf_tmp <- fcdf_tmp[!ff,]

  tmp_cdf <- cdf %>%
    filter(trait == "bmi", downsample == prop[i])
  tmp_cdf <- fcdf_tmp %>%
    select(-P)  %>% left_join(tmp_cdf, by = c("POS", "ID"))

  z <- fcdf_tmp[, 7:9]
  p <- fcdf_tmp$P

  indep_snps_rand <- fcdf_tmp$indep_snp_rand

  # apply sffdr
  out <- apply_sffdr(p,
                     z,
                     indep_snps_rand,
                     knots = c(0.005, 0.01, 0.025, 0.05, 0.1),
                     lambda = seq(0.05, 0.9, 0.05),
                     method = "gam")

  df_cor <- data.frame(type = "correlated",
                   CHR = fcdf_tmp$`#CHROM.x`,
                   POS = fcdf_tmp$POS,
                   ID = fcdf_tmp$ID,
                   p = p,
                   p_meta = tmp_cdf$P,
                   indep_snps = indep_snps_rand,
                   downsample = prop[i],
                   fq.rand = out$fqvalues,
                   fp.rand = out$fpvalues)

  save(df_cor, file = paste0("./data/subsample/ssfdr-results-fpvalues-", i, ".rds"))
}
