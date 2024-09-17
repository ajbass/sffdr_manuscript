##################################################
### File preprocesses EGPA analysis
##################################################
set.seed(12345)

library(sffdr)
library(tidyverse)

# Load primary and informative files
pfiles <- list.files("./primary/")
ifiles <- list.files("./conditional/")[-4]
primary <- info <- NULL
for (f in pfiles) {
  print(f)
  col_select <- c("SNPID", "INFO", "CHR38", "BP38", "ALT", "REF", "ALTFREQ",   "BETA", "SE")
  tmp <- read_tsv(paste0("./primary/", f),  col_select = col_select)
  tmp$CHR38 <- as.character(tmp$CHR38)
  tmp$BP38 <- as.numeric(tmp$BP38)
  tmp <- tmp %>% filter(CHR38 %in% as.character(1:22),
                        !is.na(BP38),
                        !is.na(BETA),!is.na(SE),
                        !(CHR38 =="6" & BP38 > 24e6  & BP38 < 45e6)) %>%
    group_by(BP38, CHR38) %>%
    filter(!(length(SNPID) > 1) & !(ALT %in% c("I", "D", "NA")))

  tmp <- tmp %>% filter(INFO > 0.8)
  tmp <- tmp %>% ungroup() %>%
    mutate(P = 2 * pnorm(abs(BETA / SE), lower.tail = F)) %>%
    select(-BETA, -SE)

  for (cf in ifiles) {
    col_select <- c("SNPID",  "CHR38", "BP38", "ALT",   "BETA", "SE", "P")
    tmp2 <- read_tsv(paste0("./conditional/", cf), col_select = col_select)
    if (cf != "EOSC_Astle_27863252_1-hg38.tsv.gz") {
      tmp2 <- tmp2 %>%
        dplyr::select(SNPID, CHR38, BP38,ALT, BETA, SE, P) %>%
        dplyr::filter(CHR38 %in% as.character(1:22),
                      !is.na(BP38),
                      !is.na(BETA),!is.na(SE),
                      !(CHR38 =="6" & BP38 > 24e6  & BP38 < 45e6)) %>%
        dplyr::group_by(BP38, CHR38) %>%
        dplyr::filter(!(length(SNPID) > 1) & !(ALT %in% c("I", "D", "NA"))) %>%
        dplyr::select(-SNPID, -ALT) %>%
        mutate(P = 2 * pnorm(abs(BETA / SE), lower.tail = F)) %>%
        dplyr::select(-BETA, -SE)
    } else {
      tmp2 <- tmp2 %>%
        dplyr::select(SNPID, CHR38, BP38,ALT, P) %>%
        dplyr::filter(CHR38 %in% as.character(1:22),
                      !is.na(BP38),
                      !is.na(P),
                      !(CHR38 =="6" & BP38 > 24e6  & BP38 < 45e6)) %>%
        dplyr::group_by(BP38, CHR38) %>%
        dplyr::filter(!(length(SNPID) > 1) & !(ALT %in% c("I", "D", "NA"))) %>%
        dplyr::select(-SNPID, -ALT)
    }

    tmp2$CHR38 <- as.character(tmp2$CHR38)
    tmp2$BP38 <- as.numeric(tmp2$BP38)
    colnames(tmp2)[colnames(tmp2) %in% "P"] <- str_split(cf, "_")[[1]][1]
    tmp <- tmp %>%
      left_join(tmp2)
  }

  trait <- str_split(f, "_")[[1]][1]
  write_tsv(tmp, file = paste0("processed/", paste0(trait, ".tsv")))
  rm(tmp)
}
