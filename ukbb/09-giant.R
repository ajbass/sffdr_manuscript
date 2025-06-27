###############################
### Applying sfFDR to GIANT ###
###############################
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

giant <- read_delim("./SNP_gwas_mc_merge_nogc.tbl.uniq", "\t") %>%
 filter(!is.na(p))

colnames(giant) <- c("snpid",  "a1", "a2", "freq", "beta", "se", "giant", "n")
giant <- giant %>% mutate(z = beta/se)

# Go through each downsample proportion and apply sffdr
prop <- seq(0.1, 1, 0.1)
for (i in 10:1) {
  set.seed(i)
  # randomly select representative SNP within each block
  tmp_pdf <- pdf %>%
    filter(trait == "bmi", downsample == prop[i])
  
  fcdf_tmp <- tmp_pdf %>%
    left_join(giant %>% select(snpid, giant), by = c("snpid")) %>%
    filter(!is.na(giant), !is.na(pval)) %>%
    group_by(snpid) %>% 
    filter(length(snpid) == 1) %>%
    left_join(df_LD %>% select(-bpos), by = c("snpid"), multiple = "first") %>%
    filter(!is.na(LD)) %>%
    group_by(LD) %>%
    mutate(indep_snps = create_bool(cbind(giant)))
  
  tmp_cdf <- cdf %>%
    filter(trait == "bmi", downsample == prop[i])
  tmp_cdf <- fcdf_tmp %>%
    select(-pval)  %>%
    left_join(tmp_cdf %>% select(snpid, chr, a1, a2, bpos, pval)) 
  
  z <- fcdf_tmp[, 14,drop = F]
  p <- fcdf_tmp$pval
  
  indep_snps <- fcdf_tmp$indep_snps
  
  # apply sffdr
  t1 <- proc.time()[3]
  out <- apply_sffdr(p,
                     z,
                     indep_snps,
                     knots = c(0.01, 0.025, 0.05, 0.1),
                     lambda = seq(0.05, 0.9, 0.05),
                     method = "gam")
  t2 <- proc.time()[3] - t1
  
  # q-values for raw + meta-analysis approach
  meta <- sffdr:::gwasQvalue(tmp_cdf$pval, indep_snps = indep_snps, pi0.method = "bootstrap")
  marg <- sffdr:::gwasQvalue(p, indep_snps = indep_snps, pi0.method = "bootstrap")
  
  # # Run MTAG
  write_delim(fcdf_tmp %>%
                mutate(z = beta / se) %>%
                ungroup() %>%
                distinct() %>%
                select(snpid, chr, bpos, a1, a2, freq, n, beta, se, z, pval),
              file = "./data/bmi.txt", delim = " ")

  write_delim(giant %>% ungroup() %>%
                left_join(fcdf_tmp %>% ungroup() %>% select(snpid, chr, bpos, trait)) %>%
                ungroup() %>% 
                filter(!is.na(trait)) %>% 
                rename(pval = giant) %>%
                ungroup() %>%
                select(snpid, chr, bpos, a1, a2, freq, n, beta, se, z, pval),
              file = "./data/giant.txt", delim = " ")
  mtag_path <- "./mtag/mtag.py"
  arg_sumstats <- paste0("./data/bmi.txt,", "./data/giant.txt")
  output <- paste0("./data/output_MTAG/09-bmi-", i)
  cmd <- paste("python2.7", mtag_path,
               "--sumstats", arg_sumstats,
               "--n_min 0",
               "--out", output,
               "--stream_stdout")
  # # Run
  t1 <- proc.time()[3]
  system(cmd)
  t4 <- proc.time()[3] - t1
  
  # # Save subset of results to save space
  tmp <- read_tsv(paste0(output, "_trait_1.txt"))
  tmp <- tmp %>%
    select(SNP, CHR, BP, mtag_pval) %>%
    rename(ID = SNP, POS = BP)
  system(paste0("rm ", output, "*"))
  
  df_cor <- data.frame(type = "correlated",
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
                       time.sfFDR_constrained = t2,
                       time.MTAG = t4)
  
  df_cor <- df_cor %>% left_join(tmp)
  save(df_cor, file = paste0("./data/09-giant-", i, ".rds"))
}
