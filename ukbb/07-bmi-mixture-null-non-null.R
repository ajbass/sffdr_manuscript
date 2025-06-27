#####################################################
### Applying sfFDR to mixture null/non-null traits ##
#####################################################
source("../00-helper.R")

library(locfit)
library(qvalue)
library(splines)
library(sffdr)
library(tidyverse)

# load LD blocks
df <- readRDS("../ukbb_ldblocks_01.rds")
df$BP38 <- as.numeric(unlist(str_split(df$alt38, "[:_]", simplify = T))[,2])
df_LD <- as_tibble(df) %>% dplyr::select(rsid, BP38, group)
colnames(df_LD) <- c("snpid", "bpos", "LD")

# load summary statistics
pdf <- readRDS("./data/mtag_pcombined_pvals.rds")
fcdf <- readRDS("./data/mtag_cfcombined_pvals.rds")
cdf <- readRDS("./data/mtag_ccombined_pvals.rds")

triglycerides <- fcdf %>%
  filter(downsample ==  1, trait == "triglycerides") %>%
  mutate(z = beta / se)
bfp <- fcdf %>% 
  filter(downsample == 1, trait == "bfp") %>%
  mutate(z = beta / se)
cholesterol <- fcdf %>% 
  filter(downsample ==  1, trait == "cholesterol") %>%
  mutate(z = beta / se)


# Reorganize data
df <- cbind(bfp[, 3:7],
            data.frame(triglycerides = triglycerides$pval,
                       bfp = bfp$pval,
                       cholesterol = cholesterol$pval))
fcdf_cor_all <- df %>%
  filter(!is.na(bfp), !is.na(cholesterol), !is.na(triglycerides)) %>% distinct()

# Load null traits.
fcdf_null <- readRDS("./data/mtag_null_pvals.rds")

triglycerides_null <- fcdf_null %>%
  filter(downsample ==  1, trait == "triglycerides") %>%
  mutate(z = beta / se)
bfp_null <- fcdf_null %>% 
  filter(downsample == 1, trait == "bfp") %>%
  
  mutate(z = beta / se)
cholesterol_null <- fcdf_null %>% 
  filter(downsample ==  1, trait == "cholesterol") %>%
  mutate(z = beta / se)

# Reorganize data
df_null <- cbind(bfp_null[, c(3:7,13)],
            data.frame( triglycerides_null = triglycerides_null$pval,
                        bfp_null = bfp_null$pval,
                        cholesterol_null = cholesterol_null$pval))
fcdf_cor_null_all <- df_null %>%
  filter(!is.na(bfp_null), !is.na(cholesterol_null), !is.na(triglycerides_null)) %>% distinct()

# Go through each downsample proportion and apply sffdr
prop <- seq(0.1, 1, 0.1)

for (iii in 1:3) { # non-null traits
  fcdf_cor <- fcdf_cor_all[, c(1:5, 6:(5+iii))]
  print(colnames(fcdf_cor))
  for (jjj in 1:3) { # null trait
    fcdf_cor_null <- fcdf_cor_null_all %>% filter(replicate == 1)
    fcdf_cor_null <- fcdf_cor_null[, c(1:5, 7:(6+jjj))]
    print(colnames(fcdf_cor_null))
    lfile <- NULL
    for (i in 10:1) {
      set.seed(i)
      # randomly select representative SNP within each block
      tmp_pdf <- pdf %>%
        filter(trait == "bmi", downsample == prop[i]) 
      
      names <- c(colnames(fcdf_cor)[-c(1:5)], colnames(fcdf_cor_null)[-c(1:5)])
      fcdf_tmp <- tmp_pdf %>%
        left_join(fcdf_cor, by = c("bpos", "snpid", "chr", "a1", "a2")) %>%
        left_join(fcdf_cor_null, by = c("bpos", "snpid", "chr", "a1", "a2")) %>%
        filter(!is.na(triglycerides), !is.na(pval)) %>%
        left_join(df_LD %>% select(-bpos), by = c("snpid"), multiple = "first") %>%
        filter(!is.na(LD)) 
      
      fcdf_tmp <- fcdf_tmp %>%
        group_by(LD) %>%
        mutate(indep_snps = create_bool(eval(parse(text = paste0("cbind(",paste0(names, collapse=","), ")")))))

      tmp_cdf <- cdf %>%
        filter(trait == "bmi", downsample == prop[i])
      
      tmp_cdf <- fcdf_tmp %>%
        select(-pval)  %>%
        left_join(tmp_cdf %>% select(snpid, chr, a1, a2, bpos, pval)) 
      
      z <- fcdf_tmp[, 14:(13+iii+jjj)]
      p <- fcdf_tmp$pval
      
      indep_snps <- fcdf_tmp$indep_snps
      
      # apply sffdr
      t1 <- proc.time()[3]
      out <- apply_sffdr(p,
                         z,
                         indep_snps,
                         knots = c(0.005, 0.01, 0.025, 0.05, 0.1),
                         lambda = seq(0.05, 0.9, 0.05),
                         method = "gam")
      t2 <- proc.time()[3] - t1
      
      # q-values for raw + meta-analysis approach
      meta <- sffdr:::gwasQvalue(tmp_cdf$pval, indep_snps = indep_snps, pi0.method = "bootstrap")
      marg <- sffdr:::gwasQvalue(p, indep_snps = indep_snps, pi0.method = "bootstrap")
      
      df_cor <- data.frame(type = "correlated",
                           CHR = fcdf_tmp$chr,
                           POS = fcdf_tmp$bpos,
                           ID = fcdf_tmp$snpid,
                           replicate = 1,
                           num_non_null = iii,
                           num_null = jjj,
                           p = p,
                           p_meta = tmp_cdf$pval,
                           indep_snps = indep_snps,
                           downsample = prop[i],
                           q = marg$qvalues,
                           meta_q =  meta$qvalues,
                           fq.rand = out$fqvalues,
                           fp.rand = out$fpvalues, 
                           time.sfFDR_constrained = t2)
      
      df_cor <- df_cor %>% filter(p < 5e-7 | p_meta < 5e-7 | fp.rand < 5e-7)
      
      save(df_cor, file = paste0("./data/07-sffdr-mixture-", iii, "-", jjj, "-", i, ".rds"))
    }
  }
}
