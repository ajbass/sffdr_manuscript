# Data preparation
files <- list.files("../ukbb/data/assoc/",full.names = T)
pdf <- cdf <- fcdf <- NULL
prop <- seq(0.1, 1, 0.1)
i <- j <- 1
for (f in files) {
  print(j / length(files))
  j <- j + 1
  out <- read_tsv(f)
  str <-  str_split(str_split(f, pattern = "subset")[[1]][2], "\\.")[[1]]
  i <- str[1]
  trait <- str[2]
  maf <- str_split(str_split(f, pattern = "\\.")[[1]][3], pattern = "/assoc/")[[1]][2]
  out2 <- read_tsv(paste0("../../../cam_postdoc/ukbb/scripts/maf/", maf, ".afreq"))
  type <- str_split(str_split(f, pattern = "subset")[[1]][1], "_")[[1]][2]
  out$type = type
  out$trait <- trait
  out <- out %>% left_join(out2 %>% select("#CHROM", "ID", "ALT_FREQS"), by = c("#CHROM", "ID"))
  colnames(out) <- c("chr", "bpos", "snpid", "a2", "a1", "A1", "TEST", "n", "beta", "se", "T_STAT", "pval", "type", "trait", "freq")
  out <- out %>% select(trait, type, snpid, chr, bpos, a1, a2, freq, beta, se, pval, n)
  if (type == "p") {
    if (i != 0) {
      out$downsample <- prop[as.numeric(i)]
      pdf <- dplyr::bind_rows(out, pdf)
    }
  } else {
    if (i != 0) {
      out$downsample <- prop[as.numeric(i)]
      cdf <- dplyr::bind_rows(out, cdf)
    } else {
      out$downsample <- 1
      fcdf <- dplyr::bind_rows(out, fcdf)
    }
  }
}

pdf <- pdf %>%
  filter(!(chr == 6 & bpos < 45e6 & bpos > 24e6),
         !(nchar(a1) > 1 | nchar(a2) > 1), 
         !((a1 == "A" & a2 == "T")  | (a1 == "T" & a2 == "A")),
         !((a1 == "C" & a2 == "G")  | (a1 == "G" & a2 == "C"))) %>%
  group_by(trait, type, downsample, snpid)  %>%
  filter(length(snpid) == 1) %>% 
  ungroup() %>% 
  group_by(trait, type, downsample, bpos, chr) %>%
  filter(length(snpid) == 1) %>%
  ungroup()

cdf <- cdf %>%
    filter(!(chr == 6 & bpos < 45e6 & bpos > 24e6),
           !(nchar(a1) > 1 | nchar(a2) > 1), 
           !((a1 == "A" & a2 == "T")  | (a1 == "T" & a2 == "A")),
           !((a1 == "C" & a2 == "G")  | (a1 == "G" & a2 == "C"))) %>%
    group_by(trait, type, downsample, snpid)  %>%
    filter(length(snpid) == 1) %>% 
    ungroup() %>% 
    group_by(trait, type, downsample, bpos, chr) %>%
    filter(length(snpid) == 1) %>%
    ungroup()
fcdf <- fcdf %>%
  filter(!(chr == 6 & bpos < 45e6 & bpos > 24e6),
         !(nchar(a1) > 1 | nchar(a2) > 1), 
         !((a1 == "A" & a2 == "T")  | (a1 == "T" & a2 == "A")),
         !((a1 == "C" & a2 == "G")  | (a1 == "G" & a2 == "C"))) %>%
  group_by(trait, type, downsample, snpid)  %>%
  filter(length(snpid) == 1) %>% 
  ungroup() %>% 
  group_by(trait, type, downsample, bpos, chr) %>%
  filter(length(snpid) == 1) %>%
  ungroup()

saveRDS(pdf, "./data/mtag_pcombined_pvals.rds")
saveRDS(cdf, "./data/mtag_ccombined_pvals.rds")
saveRDS(fcdf, "./data/mtag_cfcombined_pvals.rds")
rm(pdf)
rm(cdf)
rm(fcdf)

# Merge null summary statistics
library(tidyverse)
files <- list.files("../ukbb/data/assoc_null/",full.names = T)
null_cdf <- NULL
i <- 1
for (f in files) {
  out <- read_tsv(f)
  str <-  str_split(str_split(f, pattern = "subset")[[1]][2], "\\.")[[1]]
  i <- str[1]
  trait <- str[2]
  maf <- str_split(str_split(f, pattern = "\\.")[[1]][3], pattern = "/assoc_null/")[[1]][2]
  maf <- paste0(str_split(maf, "_")[[1]][-2], collapse="_")
  out2 <- read_tsv(paste0("../../../cam_postdoc/ukbb/scripts/maf/", maf, ".afreq"))
  type <- str_split(str_split(f, pattern = "subset")[[1]][1], "_")[[1]][2]
  out$type = type 
  out$trait <- trait
  out <- out %>% left_join(out2 %>% select("#CHROM", "ID", "ALT_FREQS"), by = c("#CHROM", "ID"))
  colnames(out) <- c("chr", "bpos", "snpid", "a2", "a1", "A1", "TEST", "n", "beta", "se", "T_STAT", "pval", "type", "trait", "freq")
  out <- out %>% select(trait, type, snpid, chr, bpos, a1, a2, freq, beta, se, pval, n)
  
  out$replicate <- as.numeric(i)
  out$downsample <- 1
  
  null_cdf <- dplyr::bind_rows(out, null_cdf)
}

null_cdf <- null_cdf %>%
  filter(!(chr == 6 & bpos < 45e6 & bpos > 24e6),
         !(nchar(a1) > 1 | nchar(a2) > 1), 
         !((a1 == "A" & a2 == "T")  | (a1 == "T" & a2 == "A")),
         !((a1 == "C" & a2 == "G")  | (a1 == "G" & a2 == "C"))) %>%
  group_by(trait, type, replicate, downsample, snpid)  %>%
  filter(length(snpid) == 1) %>% 
  ungroup() %>% 
  group_by(trait, type, replicate, downsample, bpos, chr) %>%
  filter(length(snpid) == 1) %>%
  ungroup()

saveRDS(null_cdf, "./data/mtag_null_pvals.rds")
