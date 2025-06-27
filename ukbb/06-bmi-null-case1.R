#############################################
### Applying sfFDR and MTAG to null case 1 ##
#############################################
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
fcdf <- readRDS("./data/mtag_null_pvals.rds")
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
df <- cbind(bfp[, c(3:5, 13)],
            data.frame(bfp = bfp$pval,
                       cholesterol = cholesterol$pval,
                       triglycerides = triglycerides$pval))
fcdf_cor <- df %>%
  filter(!is.na(bfp), !is.na(cholesterol), !is.na(triglycerides)) %>% distinct()

for (pp in 1:10) {
  for (i in  1:10) {
    set.seed(i)
    prop <- unique(pdf$downsample)[pp]
    tmp_pdf <- pdf %>%
      filter(trait == "bmi", downsample == prop) 
    
    fcdf_tmp <- tmp_pdf %>%
      left_join(fcdf_cor %>%
                  filter(replicate == i),
                by = c("bpos", "snpid", "chr")) %>%
      filter(!is.na(bfp), !is.na(pval)) %>%
      left_join(df_LD %>% select(-bpos), by = c("snpid"), multiple = "first") %>%
      filter(!is.na(LD)) %>%
      group_by(LD) %>%
      mutate(indep_snps = create_bool(cbind(bfp, cholesterol, triglycerides)))
     # mutate(indep_snps = create_bool(cbind(bfp)))
    
    z <- fcdf_tmp[, 15:17]
    # BFP only
   # z <- fcdf_tmp[, 15]
    p <- fcdf_tmp$pval
    
    indep_snps <- fcdf_tmp$indep_snps
    
    # apply sffdr
    t1 <- proc.time()[3]
    out <- apply_sffdr(p,
                       z,
                       indep_snps,
                       maxk = 50000,  
                       epsilon = min(p),
                       knots = c(0.005, 0.01, 0.025, 0.05, 0.1),
                       lambda = seq(0.05, 0.9, 0.05),
                       method = "gam")
    t2 <- proc.time()[3] - t1

    # Run MTAG
    write_delim(triglycerides %>% filter(replicate == i) %>% select(-trait, -type, -downsample, -replicate), "./data/tri.txt", delim = " ")
    write_delim(bfp %>% filter(replicate == i) %>% select(-trait, -type, -downsample, -replicate), "./data/bfp.txt", delim = " ")
    write_delim(cholesterol %>% filter(replicate == i) %>% select(-trait, -type, -downsample, -replicate), "./data/cho.txt", delim = " ")
    
    write_delim(fcdf_tmp %>%
                  mutate(z = beta / se) %>%
                  ungroup() %>% 
                  select(snpid, chr, bpos, a1, a2, freq, beta, se, pval, n, z), 
                file = "./data/bmi.txt", delim = " ")
    mtag_path <- "./mtag/mtag.py"
       
    arg_sumstats <- paste0("./data/bmi.txt,", "./data/tri.txt,./data/bfp.txt,./data/cho.txt")
    # BFP only
    # arg_sumstats <- paste0("./data/bmi.txt,", "./data/bfp.txt")
    output <- paste0("./data/output_MTAG/TMP-", prop,"-", i)
    cmd <- paste("python2.7", mtag_path,
                 "--sumstats", arg_sumstats,
                 "--out", output, 
                 "--p_name", "pval",
                 "--force",
                 "--stream_stdout") 
    # Run
    t1 <- proc.time()[3]
    system(cmd)
    t4 <- proc.time()[3] - t1
    
    # Save subset of results to save space
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
                         indep_snps = indep_snps,
                         downsample = prop,
                         replicate = i,
                         fq.rand = out$fqvalues,
                         fp.rand = out$fpvalues,
                         time.sfFDR_constrained = as.numeric(t2),
                         time.MTAG = as.numeric(t4))
    
    df_cor <- df_cor %>% left_join(tmp)
    save(df_cor, file = paste0("./data/06-null-bmi-tri-bfp-cho-", i, "-", pp, ".rds"))
    # Uncomment to save BFP setting 
   # save(df_cor, file = paste0("./data/06-null-bmi-bfp-", i, "-", pp, ".rds"))  
  }
}
