#################################
### Applying sfFDR w/ pruning ###
#################################
source("../00-helper.R")

library(tidyverse)
library(locfit)
library(qvalue)
library(splines)
library(sffdr)

# load summary statistics
pdf <- readRDS("./data/mtag_pcombined_pvals.rds")
fcdf <- readRDS("./data/mtag_cfcombined_pvals.rds")
cdf <- readRDS("./data/mtag_ccombined_pvals.rds")

triglycerides <- fcdf %>%
  filter(downsample ==  1, trait == "triglycerides") %>%
  mutate(z = beta / se)
bfp <- fcdf %>% 
  filter(downsample ==  1, trait == "bfp") %>%
  mutate(z = beta / se)
cholesterol <- fcdf %>% 
  filter(downsample ==  1, trait == "cholesterol") %>%
  mutate(z = beta / se)

# Reorganize data
df <- cbind(bfp[, 3:7],
            data.frame(bfp = bfp$pval,
                       cholesterol = cholesterol$pval,
                       triglycerides = triglycerides$pval))
fcdf_cor <- df %>%
  filter(!is.na(bfp), !is.na(cholesterol), !is.na(triglycerides)) %>% distinct()

# Go through each downsample proportion and apply sffdr
prop <- seq(0.1, 1, 0.1)
for (prune in c("0.05", "0.1", "0.2", "0.3", "0.4", "all")) {
  # combine information on pruned SNPs from PLINK
  if (prune != "all") {
    pf <- list.files(paste0("./prune_genotype/", prune), pattern = "in", full.names = T)
    
    tmp  <- data_frame(filename = pf) %>% 
    mutate(file_contents = map(filename,       
                               ~ read_delim(file.path(.), col_names = F, delim = "\t"))) %>%
      rowwise() %>% 
      unnest(file_contents) %>%
      select(-filename) %>%
      rename(snpid = X1) %>% 
      mutate(LD = TRUE)
     pf <- list.files(paste0("./prune_genotype/", prune), pattern = "out", full.names = T)
  
    tmp1  <- data_frame(filename = pf) %>% 
      mutate(file_contents = map(filename,       
                                 ~ read_delim(file.path(.), col_names = F, delim = "\t"))) %>%
      rowwise() %>% 
      unnest(file_contents) %>%
      select(-filename) %>%
      rename(snpid = X1) %>% 
      mutate(LD = FALSE)
    
    df_LD <- rbind(tmp, tmp1)
  } 
  for (i in 10:1) {
    set.seed(i)
    tmp_pdf <- pdf %>%
      filter(trait == "bmi", downsample == prop[i]) 
    
    if (prune != "all") {
      fcdf_tmp <- tmp_pdf %>%
        left_join(fcdf_cor, by = c("bpos", "snpid", "chr", "a1", "a2")) %>%
        filter(!is.na(bfp), !is.na(pval)) %>%
        left_join(df_LD, by = c("snpid"), multiple = "first") %>%
        filter(!is.na(LD))
    } else {
      fcdf_tmp <- tmp_pdf %>%
        left_join(fcdf_cor, by = c("bpos", "snpid", "chr", "a1", "a2")) %>%
        filter(!is.na(bfp), !is.na(pval)) %>%
        mutate(LD = TRUE)
    }
    
    tmp_cdf <- cdf %>%
      filter(trait == "bmi", downsample == prop[i])
    tmp_cdf <- fcdf_tmp %>%
      select(-pval)  %>%
      left_join(tmp_cdf %>% select(snpid, chr, a1, a2, bpos, pval)) 
    
    z <- fcdf_tmp[, 14:16]
    p <- fcdf_tmp$pval
    
    indep_snps <- fcdf_tmp$LD
    
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
    meta <- sffdr:::gwasQvalue(tmp_cdf$pval, indep_snps = indep_snps_rand, pi0.method = "bootstrap")
    marg <- sffdr:::gwasQvalue(p, indep_snps = indep_snps_rand, pi0.method = "bootstrap")
   
    df_cor <- data.frame(type = "correlated",
                         prune = prune,
                         CHR = fcdf_tmp$chr,
                         POS = fcdf_tmp$bpos,
                         ID = fcdf_tmp$snpid,
                         p = p,
                         p_meta = tmp_cdf$pval,
                         indep_snps = indep_snps,
                         downsample = prop[i],
                         q = marg$qvalues,
                         meta_q =  meta$qvalues,
                         fq.rand = out$fqvalues,
                         fp.rand = out$fpvalues,
                         time.sfFDR_constrained = t2)
     
    save(df_cor, file = paste0("./data/11-prune", i, "-", prune, ".rds"))
  }
}
