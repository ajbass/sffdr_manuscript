##################################################
### File downsamples the UK Biobank data ###
##################################################
# Load individuals
library(tidyverse)
data <- read_tsv("obesity_individuals.tab")

# Split into two data sets: (i) primary GWAS (ii) conditional GWAS
set.seed(12345)
n <- nrow(data) / 2
id <- sample(1:nrow(data), n, replace = F)
pdata <- data[id,]
cdata <- data[-id,]

write_tsv(pdata,
          col_names = T,
          file = paste0("pgwas.tab"))
write_tsv(cdata,
          col_names = T,
          file = paste0("cgwas.tab"))

# STEP 1 + 2: downsampling / upsampling the primary GWAS
prop <- seq(0.1, 1, 0.1)
for (i in 1:length(prop)) {
  id <- sample(1:nrow(pdata),
               replace = F,
               size = round(prop[i] * nrow(pdata)))
   sdata <- pdata[id, ]
   write_tsv(sdata, col_names = T, file = paste0("pgwas", i, ".tab"))

   # STEP 2: upsampling the conditional GWAS with chosen primary
   sdata <- rbind(cdata, pdata[id, ])
   write_tsv(sdata, col_names = T, file = paste0("cgwas", i, ".tab"))
}

# Generate null informative variables
library(tidyverse)
data <- read_tsv("obesity_phenotypes.tab")
indv <- read_tsv("./subset/cgwas.tab")
indv <- indv %>% left_join(data)
for (ii in 1:10) {
  set.seed(ii)
  null_data <- indv
  id <- sample(1:nrow(indv), replace = FALSE)
  null_data[,-c(1:2)] <- indv[id, -c(1:2)]
  write_tsv(null_data, col_names = T, file = paste0("cgwas_pheno_null", ii,".tab"))
}
