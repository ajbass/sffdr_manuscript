set.seed(12345)

library(sffdr)
library(tidyverse)

source("../00-helper.R")
df <- readRDS("../ukbb_ldblocks_01.rds")
df <- as_tibble(df) %>% dplyr::select(rsid, group)
colnames(df) <- c("snpid", "LD")

# Juvenile arthritis from FinnGen
finn_r5 <- read_tsv("./data/finngen_R5_JUVEN_ARTHR.gz", col_select = c(1:7, 9:10, 11:13)) 
colnames(finn_r5) <- c("chr", "bpos", "a2", "a1", "snpid", "nearest_gene", "pval", "beta", "se", "freq", "af_alt_cases", "af_alt_controls")
finn_r5 <- finn_r5 %>% 
  mutate(maf = ifelse(freq > 0.5, 1 - freq, freq)) %>% 
  filter(maf > 0.01) %>% 
  filter(!(chr =="6" & bpos > 24e6  & bpos < 45e6),
         chr %in% 1:22)

# informative traits RA + hypothyroid
ra <- read_tsv("./RA_Okada_24390342_1-hg38.tsv.gz", col_select = c(1:7, 7, 11))
ra <- ra %>%
  rename(snpid = SNPID,
         chr = CHR38,
         bpos = BP38, 
         a1 = ALT,
         a2 = REF,
         ra = P)  %>% 
  filter(!(chr =="6" & bpos > 24e6  & bpos < 45e6), chr %in% 1:22) 

hypothyroid <- read_tsv("./hypothyroid-gcst90013893-37.tsv", col_select = c(2:6, 11,15))
hypothyroid <- hypothyroid %>% 
  mutate(maf_ht = ifelse(hm_effect_allele_frequency > 0.5, 1 - hm_effect_allele_frequency, hm_effect_allele_frequency)) %>% 
  rename(snpid = hm_rsid,
         chr = hm_chrom,
         bpos = hm_pos,
         a1 = hm_effect_allele,
         a2 = hm_other_allele,
         hypothyroid = p_value)  %>% 
  filter(!(chr =="6" & bpos > 24e6  & bpos < 45e6), chr %in% 1:22) 

# merge data
eco_merge <- hypothyroid %>%
  left_join(ra %>% mutate(chr = as.numeric(chr)),
            multiple = "first", by = c("snpid", "chr", "bpos", "a1", "a2"))  %>%
  filter(!is.na(ra))

tmp <- finn_r5 %>% 
  mutate(chr = as.numeric(chr)) %>%
  left_join(eco_merge %>% select(chr, bpos, a2, a1, snpid, hypothyroid, ra) %>%
              mutate(chr = as.numeric(chr))) %>%
  left_join(df %>% select(snpid, LD), multiple = "first") %>% 
  filter(!is.na(hypothyroid), !is.na(pval), !is.na(ra))

# select independent SNPs
out_tmp2 <-  tmp %>%
  filter(is.na(LD)) %>%
  mutate(indep_snp_rand = FALSE,
         indep_snp_inform = FALSE)

out_tmp <-  tmp %>% filter(!is.na(LD)) %>%
  group_by(LD) %>%
  mutate(indep_snp_rand = sample(c(rep(FALSE, length(LD) - 1), TRUE)),
         indep_snp_inform = create_bool(cbind(ra, hypothyroid)))
out_train <- rbind(out_tmp, out_tmp2)

z <- out_train[, 14:15]

p <-  out_train$pval 
indep_snps <- out_train$indep_snp_inform

# apply sffdr 
out_sffdr <- apply_sffdr(p,
                         z,
                         indep_snps,
                         epsilon = min(p), 
                         lambda = seq(0.05, 0.9, 0.05),
                         knots = c(0.005, 0.01, 0.025, 0.05, 0.1))

out_train$fpvalues <- out_sffdr$fpvalues
out_train$fqvalues <- out_sffdr$fqvalues
out_train$fpi0 <- out_sffdr$fpi0
out_train$flfdr <- out_sffdr$flfdr

saveRDS(out_train, "./04-juven-r5.rds")
