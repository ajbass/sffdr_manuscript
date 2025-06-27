#################################
### Applying sfFDR to FinnGen ###
#################################
source("../00-helper.R")

library(tidyverse)
library(locfit)
library(qvalue)
library(splines)
library(sffdr)

# load LD blocks
df <- readRDS("../ukbb_ldblocks_01.rds")
df$BP38 <- as.numeric(unlist(str_split(df$alt38, "[:_]", simplify = T))[,2])
df_LD <- as_tibble(df) %>% dplyr::select(rsid, BP38, group)
colnames(df_LD) <- c("snpid", "bpos", "LD")

# load summary statistics
pdf <- readRDS("./data/mtag_pcombined_pvals.rds")
fcdf <- readRDS("./data/mtag_cfcombined_pvals.rds")
cdf <- readRDS("./data/mtag_ccombined_pvals.rds")

finngen <- read_delim("./finngen_R12_BMI_IRN.gz", "\t") 
colnames(finngen) <- c("chr", "bpos", "a1", "a2", "snpid", "ngene", "finngen", "mlogp", "beta", "se", "freq")
finngen <- finngen %>% mutate(z = beta / se)
 
# Go through each downsample proportion and apply sffdr
prop <- seq(0.1, 1, 0.1)
for (i in 10:1) {
  set.seed(i)
  # randomly select representative SNP within each block
  tmp_pdf <- pdf %>%
    filter(trait == "bmi", downsample == prop[i]) 
  # 
  fcdf_tmp <- tmp_pdf %>%
    left_join(finngen %>% select(snpid, finngen), by = c( "snpid" )) %>%
    filter(!is.na(finngen), !is.na(pval)) %>%
    group_by(snpid) %>% 
    filter(length(snpid) == 1) %>%
    left_join(df_LD %>% select(-bpos), by = c("snpid"), multiple = "first") %>%
    filter(!is.na(LD)) %>% 
    group_by(LD) %>%
    mutate(indep_snps = create_bool(cbind(finngen)))
  
  tmp_cdf <- cdf %>%
    filter(trait == "bmi", downsample == prop[i])
  tmp_cdf <- fcdf_tmp %>%
    select(-pval)  %>%
    left_join(tmp_cdf %>% select(snpid, chr, a1, a2, bpos, pval)) 
  
  z <- fcdf_tmp[, 14, drop = F]
  p <- fcdf_tmp$pval
  
  indep_snps <- fcdf_tmp$indep_snps
 
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
  
  # # Run MTAG
  write_delim(fcdf_tmp %>%
                mutate(z = beta / se) %>%
                ungroup() %>%
                distinct() %>%
                select(snpid, chr, bpos, a1, a2, freq, n, beta, se, z, pval),
              file = "./data/bmi.txt", delim = " ")
  
  # # Run MTAG
  write_delim(finngen %>% ungroup() %>%
                select(snpid, a1, a2, finngen, freq, beta, se) %>% 
                left_join(fcdf_tmp %>% ungroup() %>% select(snpid, chr, bpos, trait)) %>%
                ungroup() %>% 
                mutate(n = 500348) %>% #FinnGenn sample size
                filter(!is.na(trait)) %>% 
                rename(pval = finngen) %>%
                mutate(z = beta / se) %>%
                ungroup() %>%
                select(snpid, chr, bpos, a1, a2, freq, n, beta, se, z, pval),
              file = "./data/finngen.txt", delim = " ")
  
  mtag_path <- "./mtag/mtag.py"
  arg_sumstats <- paste0("./data/bmi.txt,", "./data/finngen.txt")
  output <- paste0("./data/output_MTAG/08-bmi-", i)
  cmd <- paste("python2.7", mtag_path,
               "--sumstats", arg_sumstats,
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
  save(df_cor, file = paste0("./data/08-finngen-", i, ".rds"))
}
