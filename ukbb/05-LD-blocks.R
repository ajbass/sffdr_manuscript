##################################################
### Create independent LD blocks ###
##################################################
library(tidyverse)
library(magrittr)
library(data.table)

DIR <- "/home/ab3105/rds/rds-basis-YRWZsjDGyaU/People/CHRIS/coloc-all/" # root of all output
files <- list.files("/home/ab3105/rds/rds-basis-YRWZsjDGyaU/People/CHRIS/coloc-all/reference/byblock/")
chr  <- sapply(files, FUN = function(x) str_split(x, "_")[[1]][1])
block <- sapply(files, FUN = function(x) str_split(x, "_")[[1]][2])
df <- data.frame(chr = chr, block = block) %>% distinct()
tot <- tot2 <- 0
tmp <- "./del"
for (i in 1:nrow(df)) {
  tmp2 <- df[i,]
  str = paste0(tmp2$chr, "_", tmp2$block)
  args <- list(block = str)

  snps37 <- fread(file.path(DIR, "reference", "byblock", paste0(args$block, "_ukbb_hg37.bim")))
  snps38 <- fread(file.path(DIR, "reference", "byblock", paste0(args$block, "_ukbb_hg38.bim")))
  setnames(snps37, c("chr", "rsid", "kk", "bp37", "a1", "a2"))
  setnames(snps38, c("chr", "rsid", "kk", "bp38", "a1", "a2"))
  snps <- cbind(snps37[, .(rsid, chr, bp37, a1, a2)], snps38[, .(bp38)])
  snps[, alt37 := paste0(chr, ":", bp37, "_", a1, "_", a2)] # these are all the block snps in build 37
  snps[, alt38 := paste0(chr, ":", bp38, "_", a1, "_", a2)] # these are all the block snps in build 38

  ## LD and MAF
  stub <- file.path(DIR,"reference","byblock", paste0(args$block,"_ukbb40k_hg37"))
  paste0("plink --r square --keep-allele-order --bfile ", stub," --out ", tmp) %>%
    system()
  LD <- paste0(tmp,".ld")  %>% scan(., what=0)
  n <- length(LD)  %>% sqrt()
  LD  %<>%  matrix(., n, n, byrow=TRUE)
  LDsnps <- paste0(stub,".bim")  %>% fread()
  LDsnps[, alt37:=paste0(V1,":",V4,"_",V5,"_",V6)] # these are the colnames of LD, in build 37, you can match to snps above to get build 38 identifiers
  colnames(LDsnps) <- c("chr", "rsid", "V3", "bp37", "a1", "a2", "alt37")
  LDsnps <- LDsnps %>% left_join(snps) %>% select(chr, rsid, alt38)

  D <- as.dist(1 - LD ^ 2)
  hc <- hclust(D)
  cuts <- cutree(hc, h = 1 - 0.01 ^ 2)
  LDsnps$group <- paste0(str, "_", cuts)
  print(length(unique(cuts)))
  saveRDS(LDsnps, file = paste0("~/rds/hpc-work/LDblocks/ukbiobank/", str, ".rds") )
}
