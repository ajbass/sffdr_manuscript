##################################################
### Applying sfFDR to UK Biobank null data ###
##################################################
source("../00-helper.R")
set.seed(12345)

library(tidyverse)
library(locfit)
library(qvalue)
library(splines)
library(fFDR)
library(sffdr)
library(tidyverse)

# data preparation: p-values from PLINK on permuted phenotypes
library(tidyverse)
files <- list.files("./data/assoc_null/",full.names = T)
pdf <- fcdf <- NULL
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

# load LD blocks
ld <- readRDS("../ukbb_ldblocks_01.rds")

ld <- ld %>%
  select(rsid,  group)

colnames(ld) <- c("ID",  "LD")
cdf2 <- cdf %>%
  left_join(ld, multiple = 'first')

saveRDS(cdf2, "./null_data.rds")

# Load data
pdf <- readRDS("./pcombined_pvals.rds")
fcdf <- readRDS("./null_data.rds")
pdf <- pdf %>%
  select(trait, downsample, `#CHROM`, POS, ID, indep_snps, P, type) %>%
  left_join(ld, multiple = "first")
fcdf <- fcdf %>%
  select(trait, replicate, downsample, `#CHROM`, POS, ID,  P, LD)

bfp <- fcdf %>%
  filter(trait == "bfp")
cholesterol <- fcdf %>%
  filter(trait == "cholesterol")
triglycerides <- fcdf %>%
  filter(trait== "triglycerides")
df <- cbind(bfp[, 1:8],
            data.frame(bfp = bfp$P,
                       cholesterol = cholesterol$P,
                       triglycerides = triglycerides$P))
fcdf <- df %>%
  filter(!is.na(bfp), !is.na(cholesterol), !is.na(triglycerides)) %>% distinct()
fcdf_bmi <- fcdf %>% filter(downsample == 1)
prop = seq(0.1, 1, 0.1)
for (pp in 1:10) {
  for (i in  1:10) {
    prop <- unique(pdf$downsample)[pp]

    # Correlated traits
    tmp_pdf <- pdf %>% filter(trait == "bmi", downsample == prop)
    fcdf_bmitmp <- fcdf_bmi %>% filter(replicate == i)
    tmp_z <- tmp_pdf[, c(1, 3:5)] %>%
      left_join(fcdf_bmitmp, by = c("POS", "ID"))
    ff = !is.na(tmp_z$bfp)
    tmp_z <- tmp_z[ff,]
    tmp_pdf <- tmp_pdf[ff,]
    tmp_z2 <- tmp_z %>% filter(!is.na(LD)) %>%
      group_by(LD) %>%
      mutate(indep_snp_rand = sample(c(rep(FALSE, length(LD) - 1), TRUE)))

    indep_snps_rand <- tmp_z2 %>%
      ungroup() %>%
      filter(indep_snp_rand) %>%
      select(ID)
    indep_snps_inform <- tmp_z2 %>%
      ungroup() %>%
      filter(indep_snp_inform) %>%
      select(ID)

    z <- tmp_z[, 11:13]
    p <- tmp_pdf$P
    indep_snps_rand <- tmp_pdf$ID %in% indep_snps_rand$ID

    # apply sffdr
    out_r <- apply_sffdr(p,
                         z,
                         indep_snps_rand,
                         knots = c(0.005, 0.01, 0.025, 0.05, 0.1),
                         lambda = seq(0.05, 0.9, 0.05),
                         epsilon = min(p),
                         method = "gam")


    df_cor <- data.frame(type = "correlated",
                         ID = tmp_pdf$ID,
                         prop = prop,
                         rep = i,
                         p = p,
                         fq_r = out_r$fqvalues,
                         flfdr_r = out_r$flfdr,
                         fp_r = out_r$fpvalues,
                         indep_snps = indep_snps_rand)

    save(df_cor, file = paste0("./", i, "-", prop, "null.rds"))
    print(i)
  }
}
