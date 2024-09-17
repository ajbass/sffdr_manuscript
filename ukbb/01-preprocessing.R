##################################################
### File preprocesses UK Biobank data  ###
## Interested in obesity-related traits for
## individuals (unrelated) with British ancestry
##################################################
library(tidyverse)

path <- "/home/ab3105/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank"
path_to_phenotype <- paste0(path, "/phenotypes/ukb672224.tab")

# Grab phenotype names
df <- read_delim(path_to_phenotype, delim = "\t", n_max = 0)
names <- colnames(df)

PCs <- which(names %in% paste0("f.22009.0.", 1:20))
used_in_PCA <- which(names %in% "f.22020.0.0")
sex <- which(names %in% "f.22001.0.0")
ID <- which(names %in% "f.eid")
BMI <- which(names %in% "f.21001.0.0")
age <- which(names %in% "f.21003.0.0")
BFP <- which(names %in% "f.23099.0.0")
cholesterol <-  which(names %in% "f.30690.0.0")
tri <- which(names %in% "f.30870.0.0")

id <- c(ID, age, sex, BMI, BFP, cholesterol, tri, used_in_PCA, PCs)

## Load age + BMI
df <- read_delim(path_to_phenotype, delim = "\t",
                col_select = id)
colnames(df) <- c("IID", "age", "sex", "bmi", "bfp", "cholesterol", "triglycerides",  "used.in.pca", paste0("PC", 1:20))

id <- !is.na(df$used.in.pca)
df0 <- df[id,]
df0 <- df0 %>%
  filter(!is.na(bmi), !is.na(bfp), !is.na(cholesterol), !is.na(triglycerides))

df0 <- df0 %>%
  mutate(bmi = RNOmni::RankNorm(bmi),
         bfp = RNOmni::RankNorm(bfp),
         cholesterol = RNOmni::RankNorm(cholesterol),
         triglycerides = RNOmni::RankNorm(triglycerides))

# Organize for PLINK
out <- data.frame(IID = df0$IID,
                  FID = df0$IID,
                  sex = df0$sex,
                  age = df0$age,
                  bmi = df0$bmi,
                  bfp = df0$bfp,
                  cholesterol = df0$cholesterol,
                  triglycerides = df0$triglycerides,
                  df0[, 9:28])

# Withdrawals
withdrawals <- read_tsv("/home/ab3105/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/phenotypes/withdraw98032_19.txt", col_names = FALSE)
withdrawals <- withdrawals %>% rename(IID = X1)
out <- out %>% filter(!(IID %in% withdrawals$IID))

# save sample info for R
saveRDS(out, "~/rds/hpc-work/ukbiobank/sample_info.rds")

## Save files for plink
write_tsv(out[, c(2,1)], file = "~/rds/hpc-work/ukbiobank/obesity_individuals.tab")
write_tsv(out[,-c(5:8)][,c(2,1, 3:24)], file = "~/rds/hpc-work/ukbiobank/obesity_covariates.tab")
write_tsv(out[, c("FID", "IID", "bmi", "bfp", "cholesterol", "triglycerides")],
          file = "~/rds/hpc-work/ukbiobank/obesity_phenotypes.tab")
