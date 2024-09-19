##################################################
### Apply sfFDR to EGPA with null traits
##################################################
source("../00-helper.R")
set.seed(12345)

library(sffdr)
library(tidyverse)
library(forcats)

# Load data: p-values from PLINK
fcdf <- readRDS("../ukbb/cfcombined_pvals.rds")
fcdf <- fcdf %>% select(trait, downsample, `#CHROM`, OBS_CT, POS, ID, indep_snps, BETA, SE,  T_STAT,P, type)

# Reorganize data: use original obesity-related traits
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
df_cond <- df %>%
  filter(!is.na(bfp), !is.na(cholesterol), !is.na(triglycerides)) %>% distinct()

df_cond <- df_cond %>% dplyr::rename(CHR38 ="#CHROM", BP = POS) %>%
  select(BP, CHR38, bfp, cholesterol, triglycerides)

# load EGPA
out_raw <- read.delim(paste0("./primary/EGPA_Lyons_31719529_1-hg38.tsv.gz"), sep = "\t")
map <- out_raw %>% select(SNPID, BP38, BP, CHR38) %>% mutate(CHR38 = as.numeric(CHR38))

# load LD blocks
df <- readRDS("../ukbb_ldblocks_01.rds")
df$BP38 <- as.numeric(unlist(str_split(df$alt38, "[:_]", simplify = T))[,2])
df_LD <- as_tibble(df) %>% dplyr::select(rsid, chr, BP38, group)
colnames(df_LD) <- c("ID", "CHR38", "BP38", "LD")

primary <- c("EGPA")
ret_df <- NULL
dfc <- df_cond
for (trait in primary) {
  # Merge data sets
  out_EGPA <- read.delim(paste0("./processed/", trait, ".tsv"), sep = "\t")
  out_EGPA <- out_EGPA %>%
    left_join(df_LD) %>%
    left_join(map)

  # remove missing data
  out_EGPA2 <- dfc %>%
    left_join(out_EGPA, by = c("CHR38", "BP"))
  out_train <- out_EGPA2 %>%
    dplyr::filter(!is.na(P)) %>%
    dplyr::select(SNPID, P, bfp, cholesterol, triglycerides, LD) %>%
    dplyr::filter(!is.na(bfp) & !is.na(cholesterol) & !is.na(triglycerides))

  # select representative SNPs
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
  id <- which(rowSums(is.na(as.matrix(z))) == 0)
  p <- out_train$P[id]
  z <- z[id,]
  indep_snps_inform <- out_train$indep_snp_inform[id]

  # apply sffdr
  knots <- c(0.005, 0.01, 0.025, 0.05, 0.1)
  out_sffdr_indp_i <- apply_sffdr(p, z,
                                  indep_snps_inform,
                                  epsilon = 1e-15,
                                  method =  "gam",
                                  lambda = seq(0.05, 0.9, 0.05),
                                  knots = knots)

  fp_i <- out_sffdr_indp_i$fpvalues
  ret <-  data.frame(id = 1:length(p),
                     rep = rep,
                     p = p,
                     fp  = fp_i,
                     fq  = out_sffdr_indp_i$fqvalues,
                     fpi0 = out_sffdr_indp_i$fpi0,
                     indep_snps = indep_snps_inform,
                     trait = trait)
  ret_df <- rbind(ret, ret_df)
}

saveRDS(ret_df, file = "./04-functional-pvalue-null-inform-uncorrelated.rds")
