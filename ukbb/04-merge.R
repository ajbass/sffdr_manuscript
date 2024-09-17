##################################################
### Combines UKBB significance results ###
##################################################
library(tidyverse)

# combine PLINK significance results
files <- list.files("./data/assoc", full.names = T)
pdf <- cdf <- fcdf <- NULL
prop <- seq(0.1, 1, 0.1)
i <- 1
for (f in files) {
  out <- read_tsv(f)
  str <-  str_split(str_split(f, pattern = "subset")[[1]][2], "\\.")[[1]]
  i <- str[1]
  trait <- str[2]
  type <- str_split(str_split(f, pattern = "subset")[[1]][1], "_")[[1]][2]
  out$type <- type
  out$trait <- trait
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

saveRDS(pdf, "./pcombined_pvals.rds")
saveRDS(cdf, "./ccombined_pvals.rds")
saveRDS(fcdf, "./cfcombined_pvals.rds")
