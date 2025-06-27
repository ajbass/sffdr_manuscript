#############################
### Applying sfFDR to MVP ###
#############################
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

# Go through each downsample proportion and apply sffdr
for (pop in c("AFR", "AMR", "EAS", "EUR", "META")) {
  mvp <- read_delim(paste0("./MVP_BMI/MVP_R4.1000G_AGR.BMI_Mean_INT.", pop, ".GIA.dbGaP.txt.gz"), "\t",
                      col_select = c(1:6, 7:11)) %>%
    filter(!is.na(pval))
  mvp <- mvp %>% 
    mutate(maf = pmin(af, 1-af)) %>%
    filter(maf > 0.01) %>%
    rename(snpid = SNP_ID,
           chr = chrom,
           bpos = pos,
           a1 = alt,
           a2 = ref,
           n = num_samples,
           se = sebeta)
  mvp <- mvp %>% mutate(z = beta/se)
  # bound p-value (0,1)
  mvp$pval[mvp$pval==0] <- min(mvp$pval[mvp$pval!=0])
  mvp$pval[mvp$pval== 1] <- max(mvp$pval[mvp$pval!=1])
  prop <- seq(0.1, 1, 0.1)
  for (i in 10:1) {
    set.seed(i)
    # randomly select representative SNP within each block
    tmp_pdf <- pdf %>%
      filter(trait == "bmi", downsample == prop[i]) 
  
    fcdf_tmp <- tmp_pdf %>%
      left_join(mvp %>% mutate(mvp = pval) %>% select(snpid, mvp), by = c("snpid" )) %>%
      filter(!is.na(mvp), !is.na(pval)) %>%
      group_by(snpid) %>% 
      filter(length(snpid) == 1) %>%
      left_join(df_LD %>% select(-bpos), by = c("snpid"), multiple = "first") %>%
      filter(!is.na(LD)) %>% 
      group_by(LD) %>%
      mutate(indep_snps = create_bool(cbind(mvp)))
    
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
    
    write_delim(mvp %>% ungroup() %>% select(-chr, -bpos) %>%
                  left_join(fcdf_tmp %>% ungroup() %>% select(snpid, chr, bpos, trait), by = "snpid") %>%
                  ungroup() %>% 
                  filter(!is.na(trait)) %>% 
                  mutate(freq = af, beta =  -1 * beta, z = -1 * z) %>% # MVP inverse-rank transformed other direction
                  ungroup() %>%
                  select(snpid, chr, bpos, a1, a2, freq, n, beta, se, z, pval),
                file = "./data/mvp.txt", delim = " ")
    mtag_path <- "./mtag/mtag.py"
    arg_sumstats <- paste0("./data/bmi.txt,", "./data/mvp.txt")
    output <- paste0("./data/output_MTAG/10-bmi-", i)
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
    
    df_cor <- data.frame(pop=pop,
                         type = "correlated",
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
    save(df_cor, file = paste0("./data/12-mvp-", pop, "-", i, ".rds"))
  }
}
