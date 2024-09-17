#######################################
### Code to analyze simulation study ##
### UK Biobank study and EGPA study ###
#######################################

library(tidyverse)
library(sffdr)
library(splines)
library(latex2exp)
library(xtable)
library(forcats)
library(scales)
library(patchwork)
library(coloc)

############ UK Biobank ###############

# load data
dfc <- dfc_cor <- NULL
prop <- seq(0.1, 1, 0.1)
for (i in 1:10) {
  print(i)
  load(paste0("./ukbb/data/subsample/ssfdr-results-fpvalues-", i, ".rds"))
  dfc <- rbind(df, dfc)
  dfc_cor <- rbind(dfc_cor, df_cor)
}

# Get prediction
oracle <- dfc_cor %>%
  dplyr::group_by(downsample) %>%
  dplyr::summarise(tot = sum(p_meta < 5e-8, na.rm = T))
oracle$downsample = oracle$downsample + 1
raw0 <- dfc_cor %>%
  dplyr::group_by(downsample) %>%
  dplyr::summarise(tot = sum(p < 5e-8, na.rm = T))
df_pred <- rbind(oracle,raw0)

# Prediction model
x <- df_pred$tot
lmout <- lm(df_pred$downsample ~ ns(x, df = 3))

df_correlated <- dfc_cor %>%
  dplyr::group_by(downsample) %>%
  dplyr::filter(!is.na(p)) %>%
  dplyr::summarise(tot = sum(fp.rand < 5e-8, na.rm = T))
pred.cor <- predict(lmout, newdata = data.frame(x=df_correlated$tot))
pred.c <- pred.cor

N <- 190300 # sample size at downsample proportion = 1
plot.df <- data.frame(downsample = df_correlated$downsample,
           predicted = N *  pred.c - N * df_correlated$downsample,
           Method = rep("Functional P", each = 10))

# average predicted increase across all downsampling proportions
mean(plot.df$predicted / (N * rep(df_correlated$downsample, 1)))

p0 <- plot.df %>%
  ggplot(aes(x = downsample, y = predicted / 1000, color = Method)) +
  geom_point(size = 2.5) +
  geom_line(size = 1.3) +
  theme_bw(base_size = 15) +
  xlim(0,1) + ylim(0,115) +
  scale_color_manual("Method",
                     breaks = c("Oracle", "Functional P-value"),
                     labels = c("Oracle",  "Functional P-value"),
                     values = c("black", "black")) +
  xlab("Downsampling proportion") +theme(axis.text.x=element_text(size=13.5))+
  ylab("Additional samples (x1000)") + theme(legend.position = "none")

p1 <- dfc_cor %>%
  dplyr::group_by(downsample) %>%
  dplyr::summarise(tot0 = sum(p<5e-8, na.rm = T), tot  = sum(fp.rand < 5e-8, na.rm= T)) %>%
  ggplot(aes(x = downsample, y = tot)) +
  geom_point(aes(x=downsample, y= tot0, color = "gray"), size = 2.5) +
  geom_line(aes(x=downsample, y= tot0, color = "gray"),size = 1.3) +
  geom_point(data = df_correlated, aes(x=downsample, y= tot, color = "darkgreen"), size = 2.5) +
  geom_line(data = df_correlated, aes(x=downsample, y= tot, color = "darkgreen"), size = 1.3) +
  scale_color_manual("Method",
                    breaks = c("darkred", "darkgreen", "gray"),
               labels = c("Oracle",  "Functional p-value", "Standard p-value"),
                     values = c("black", "black", "#999999")) +
  theme_bw(base_size = 15) +theme(axis.text.x=element_text(size=13.5))+
  xlim(0.0,1) +
  ylab("Discoveries") +
  xlab("Downsampling proportion")+ theme(legend.position = c(0.26, 0.75),
                                         legend.background=element_rect(fill = alpha("white", 0)),
                                         legend.key=element_rect(fill = alpha("white", 0)))
library(patchwork)
p= p1 + p0 + plot_annotation(tag_levels = c('a'), tag_prefix = '(',
                       tag_suffix = ')')
ggsave(p,
       file = "./figures/09-applied-combined.png",
       dpi = 1000,
       width = 10,
       height = 4)

# Overlap with meta-analysis approach
status <- dfc_cor %>%
  dplyr::mutate(status = meta_fp < 5e-8) %>%
  dplyr::select(CHR, POS, ID, downsample, status)
colnames(status) <- c("CHR", "BP", "SNP", "downsample", "status")

t1 <- dfc_cor %>%
  dplyr::mutate(status = meta_fp < 5e-8) %>%
  dplyr::mutate(method = "Functional P") %>%
  dplyr::group_by(downsample, method) %>%
  dplyr::summarise(Overlap = sum(fp.rand[status] < 5e-8, na.rm = T),
            Distinct = sum(fp.rand[!status] < 5e-8, na.rm = T),
            Prop = Overlap / (Overlap + Distinct))

t2 <- dfc_cor %>%
  dplyr::mutate(status = meta_fp < 5e-8) %>%
  dplyr::mutate(method = "Raw P") %>%
  dplyr::group_by(downsample, method) %>%
  dplyr::summarise(Overlap = sum(p[status] < 5e-8, na.rm = T),
            Distinct = sum(p[!status] < 5e-8, na.rm = T),
            Prop = Overlap / (Overlap + Distinct))

df <- rbind(t1, t2)
p2 <- df %>%
  ggplot(aes(x = downsample, y = Prop, color = method)) +
  geom_point() +
  geom_line() +
  scale_color_manual("Method",
                     breaks = c("Functional P", "Raw P"),
                     labels = c("sfFDR\nfunctional p-value", "Standard p-value"),
                     values = c("black",  "grey")) +
  theme_bw() +
  xlim(0.1,1) +
  ylim(0.0,1 ) +
  ylab("Overlap") +
  xlab("Downsampling proportion")

ggsave(p2,
       file = "./figures/09-functional-pvalues-overlap_new.png",
       dpi = 1000, width = 5, height = 2.75)

# UK Biobank: sffdr w/ permuted null traits
rm(df)
rm(df_cor)
df <- NULL
prop <- seq(0.1, 1, 0.1)
for (i in prop) {
  for (ii in 1:10) {
    print(ii)
    load(paste0("./ukbb/data/sffdr-results-null-LD-", ii, "-", i, ".rds"))
    df_cor$prop <- i
    df_cor$rep <- ii
    df <- rbind(df, df_cor)
  }
}

prop <- paste0("Downsample: ", seq(0.1, 1.0, 0.1))
names(prop) <- seq(0.1, 1.0, 0.1)

dfs <- df %>%
  dplyr::select(ID, prop, rep, p, fp_r)

p3 <- dfs %>%
  dplyr::group_by(ID, prop) %>%
  dplyr::summarise(p = mean(p), fp = mean(fp_r)) %>%
  ggplot(aes(x = -log10(p), y = -log10(fp))) +
  geom_point(alpha = 0.35,
             size = 1.5) +
  theme_bw() +
  facet_wrap(~prop,
             scales = "free",
             labeller = labeller(prop = prop)) +
  xlab("p-values") +
  ylab("functional p-values") +
  geom_abline(slope = 1)

ggsave(p3,
       file = "./figures/09-functional-pvalues-null-oracle.png",
       dpi = 1000, width = 6.5, height = 6.5)

############ Simulation study: independent tests ############
fnames <- list.files("./simulation_study/data-fdr/", recursive = T)

# combine simulation files
out <- list()
for (file in fnames) {
  print(file)
  df <- readRDS(paste0("./simulation_study/data-fdr/", file))
  out <- rbind(out, df)
}

out$method <- factor(out$method, labels = c("AdaPT", "Boca-Leek", "CAMT", "Oracle", "Q-value", "sfFDR" ))
out$method <- fct_relevel(out$method,
            "Q-value",
            "Boca-Leek",
            "AdaPT",
            "CAMT",
            "sfFDR",
            "Oracle")

# sfFDR example fq-values in manuscript
out %>% filter(method == "sfFDR",
               signal.density == "Medium",
               signal.density.z == "High",
               quantity == "fqvalues",
               prior.strength == "Large",
               threshold == 0.01) %>%
  ungroup() %>%
  dplyr::summarise(mean(emptdr))

out %>% filter(method == "sfFDR",
               signal.density == "Medium",
               signal.density.z == "Low",
               quantity == "fqvalues",
               prior.strength == "Large",
               threshold == 0.01) %>%
  ungroup() %>%
  dplyr::summarise(mean(emptdr))

# sfFDR example fp-values in manuscript
out %>% filter(method == "sfFDR",
               signal.density == "Low",
               signal.density.z == "Low",
               quantity == "fpvalues",
               prior.strength == "Large",
               threshold == 5e-8) %>%
  ungroup() %>%
  dplyr::summarise(mean(discoveries, na.rm = T))

out %>% filter(method == "sfFDR",
               signal.density == "Low",
               signal.density.z == "High",
               quantity == "fpvalues",
               prior.strength == "Large",
               threshold == 5e-8) %>%
  ungroup() %>%
  dplyr::summarise(mean(discoveries))

# Figures for simulation results
cbbPalette <- c("High" = "black", "Medium"= "darkgrey", "Low" = "lightgrey")
sd.labs <- c("High power: P", "Medium power: P", "Low power: P")
names(sd.labs) <- c("High", "Medium", "Low")
z.labs <- c("High power: Z", "Medium power: Z", "Low power: Z")
names(z.labs) <- c("High", "Medium", "Low")

# Simulations evaluating fq-values
cbbPalette <- c("High" = "black", "Medium"= "darkgrey", "Low" = "lightgrey")

p1 <- out %>%
  dplyr::filter(quantity == "fqvalues",
                prior.strength != "None",
                prior.coverage == 0.025,
                !(method %in% c("sffdr_mono" ))) %>%
  dplyr::group_by(prior.strength, threshold, method, signal.density, signal.density.z, prior.coverage) %>%
  dplyr::summarise(empfdr = mean(empfdr),
            emptdr = mean(emptdr)) %>%
  ggplot(aes(x = empfdr, y = emptdr, color = method, shape = factor(prior.strength)), group = signal.density.z) +
  geom_point() +
  geom_line() +
  facet_grid(signal.density ~ signal.density.z,
             scales = "free", labeller = labeller(signal.density = sd.labs,signal.density.z = z.labs)) +
  theme_bw()  +
  scale_color_manual("Method",
                     breaks = c("Oracle", "sfFDR", "CAMT", "Boca-Leek", "AdaPT", "Q-value"),
                     labels = c("Oracle", "sfFDR", "CAMT", "Boca-Leek", "AdaPT", "Q-value"),
                     values = c("black", "#009E73", "lightblue",  "#56B4E9", "#E69F00", "#D55E00")) +
  xlab("Empirical false discovery rate") +
  ylab("Empirical true discovery rate") +
  scale_shape_manual("Prior strength", values = c(15,16,17))+ theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1))
ggsave(p1,
       file = "./figures/04-simulations-TDR-FDR.png",
       dpi = 1000,
       width = 7,
       height = 4.5)

p1 <- out %>%
  dplyr::filter(quantity == "fqvalues",
                prior.strength != "None", method %in% c("sfFDR", "Oracle", "Q-value"),
                prior.coverage == 0.025,
                !(method %in% c("sffdr_mono" ))) %>%
  dplyr::group_by(prior.strength, threshold, method, signal.density, signal.density.z, prior.coverage) %>%
  dplyr::summarise(threshold = mean(threshold),
                   discoveries = mean(discoveries)) %>%
  ggplot(aes(x = threshold, y = discoveries, color = method, shape = factor(prior.strength)), group = signal.density.z) +
  geom_point() +
  geom_line() +
  facet_grid(signal.density ~ signal.density.z,
             scales = "free", labeller = labeller(signal.density = sd.labs,signal.density.z = z.labs)) +
  theme_bw()  +
  scale_color_manual("Method",
                     breaks = c( "Q-value", "sfFDR","Oracle" ),
                     labels = c( "Standard q-value", "sfFDR\nfunctional q-value", "Oracle\nfunctional q-value"),
                     values = c("#D55E00",  "#009E73", "black")) +
  xlab("Target FDR") +guides(color = guide_legend(order = 1)) + theme(legend.key.spacing.y =  unit(2, "pt")) +
  ylab("Discoveries") +
  scale_shape_manual("Prior strength", values = c(15,16,17))
ggsave(p1,
       file = "./figures/04-simulations-fdr-discoveries.png",
       dpi = 1000,
       width = 8,
       height = 4.5)

out$method <- factor(out$method, levels = c( "Oracle","Q-value", "sfFDR", "AdaPT", "Boca-Leek", "CAMT"))
p1 <- out %>%
  dplyr::filter(quantity == "fqvalues",
                prior.coverage == 0.025,
                threshold %in% c( 0.01),
                !(method %in% c("sffdr_mono" ))) %>%
  dplyr::group_by(prior.strength, threshold, method, signal.density, signal.density.z, prior.coverage) %>%
  ggplot(aes(x = as.factor(method), y = empfdr), group = signal.density.z) +
  geom_boxplot(outlier.size = 0.5) +
  facet_grid(signal.density ~ signal.density.z, scales = "free",
             labeller = labeller(signal.density = sd.labs,signal.density.z = z.labs)) +
  theme_bw()  +
  ylab("Empirical FDR") +
  xlab("Target FDR") +
  scale_shape_manual("Prior strength", values = c(15,16,17)) + geom_hline(yintercept = 0.01, size = 0.7, color = "darkred", linetype = "dashed")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")
ggsave(p1,
       file = "./figures/04-simulations-FDR-type1-error.png",
       dpi = 1000,
       width = 7,
       height = 5)

# Evaluating pi0
out$method <- factor(out$method, levels = c("Oracle", "Q-value", "sfFDR", "AdaPT", "Boca-Leek", "CAMT"))
out_tmp <- out %>%
  dplyr::filter(method == "Oracle")
out_tmp$Oracle <- out_tmp$pi0
out_tmp <- out_tmp %>%
  dplyr::select(seed, Oracle) %>% distinct()
out2 <- out %>%
  dplyr::filter(method != "Oracle") %>%
  dplyr::left_join(out_tmp)

p1 <- out2 %>%
  dplyr::filter(quantity == "fqvalues",
                prior.coverage == 0.025,
                !(method %in% c("sffdr_mono"))) %>%
  dplyr::group_by(prior.strength, Oracle, threshold, method, signal.density, signal.density.z, prior.coverage) %>%
  ggplot(aes(x = as.factor(method), y = pi0-Oracle, color = factor(prior.strength)), group = signal.density.z) +
  geom_boxplot() +
  facet_grid(signal.density ~ signal.density.z,
             scales = "free", labeller = labeller(signal.density = sd.labs,signal.density.z = z.labs)) +
  theme_bw()  +
  scale_color_manual("Prior strength",
                     breaks = c("Large", "Moderate", "None"),
                     values = c("black", "grey40", "grey")) +
  ylab("Proportion of null tests") +
  geom_hline(yintercept = 0, color = "darkred", size = 0.8, linetype = "dashed")+
  xlab("") +
  ylab( TeX("$\\hat{\\pi}_{0} - \\pi_{0}$"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(p1,
       file = "./figures/04-simulations-pi0.png",
       width = 7.5,
       height = 5.25)

fnames <- list.files("./simulation_study/data-fdr/", recursive =  T)

# Simulation results for fp-values
outQ <- out
out <- list()
for (file in fnames) {
  print(file)
  df <- readRDS(paste0("./simulation_study/data-fdr/", file))
  out <- rbind(out, df)
}

out$method <- factor(out$method, labels = c("AdaPT", "Boca-Leek", "CAMT", "Oracle", "Raw P", "Functional P" ))
out$method <- fct_relevel(out$method,
                          "Raw P",
                          "Boca-Leek",
                          "AdaPT",
                          "CAMT",
                          "Functional P",
                          "Oracle")
out <- out %>%
  dplyr::filter(method %in% c("Oracle", "Raw P", "Functional P"))

# Figures for 3 traits
cbbPalette <- c("High" = "black", "Medium"= "darkgrey", "Low" = "lightgrey")
# New facet label names for dose variable
sd.labs <- c("High power: P", "Medium power: P", "Low power: P")
names(sd.labs) <- c("High", "Medium", "Low")

fnames <- c(TeX("Standard $p$-value"), TeX("sfFDR $p_f$-value"), TeX("Oracle $p_f$-value"))
names(fnames) <- c("Raw P", "Functional P", "Oracle")


fnames <- c(
  'Raw P' = "Standard p-value" , #@TeX("Standard $p$-value"),
  'Functional P' = "sfFDR functional p-value",#TeX("sfFDR $p_f$-value"),
  'Oracle' = "Oracle functional p-value"#TeX("Oracle $p_f$-value")
)
global_labeller <- labeller(
  signal.density = sd.labs,
  method = fnames
)

p1 <- out %>%
  dplyr::filter(quantity == "fpvalues",
                prior.coverage == 0.025,
                threshold %in% c(0.0001),
                method != "sffdr_mono") %>%
  dplyr::group_by(prior.strength, threshold, method, signal.density, signal.density.z, prior.coverage) %>%
  ggplot(aes(x = factor(prior.strength), y = empfdr,color = signal.density.z), group = signal.density.z) + geom_boxplot() +
  facet_grid(signal.density~method,
             scales = "free",labeller = global_labeller) +
  geom_abline(slope = 1) + theme_bw()  +
  xlab("Prior strength") +
  ylab("False positive rate") +#theme(legend.position = "top", legend.justification = "left") +
  scale_color_manual("Z density", values = cbbPalette)  +
  geom_hline(data = data.frame(empfdr = rep(c( 0.0001),each = 1), threshold = rep(c( 0.0001),each = 1)), aes(yintercept = empfdr), linetype = "dashed")
ggsave(p1,
       file = "./figures/04-simulations-type1error.png",
       dpi = 1000,
       width = 7,
       height = 4.5)


cbbPalette <- c("Oracle" = "black", "Functional P"= "#0072B2", "Raw P" = "#999999")
sd.labs <- c("High power: P", "Medium power: P", "Low power: P")
names(sd.labs) <- c("High", "Medium", "Low")
z.labs <- c("High power: Z", "Medium power: Z", "Low power: Z")
names(z.labs) <- c("High", "Medium", "Low")

p1 <- out %>%
  dplyr::filter(quantity == "fpvalues",
                prior.coverage == 0.025,
                threshold %in% c(5e-8),
                method != "sffdr_mono") %>%
  dplyr::group_by(prior.strength, threshold, method, signal.density, signal.density.z, prior.coverage) %>%
  ggplot(aes(x = factor(prior.strength),
             y = discoveries,
             color = method), group = signal.density.z) +
  geom_boxplot() +
  facet_grid(signal.density ~ signal.density.z,
             scales = "free", labeller = labeller(signal.density = sd.labs,signal.density.z = z.labs)) +
  theme_bw()  +
  xlab("") +
  ylab("Discoveries") +guides(color = guide_legend(order = 1)) +
  theme(legend.key.spacing.y =  unit(2, "pt")) +
  scale_color_manual("Method", values = cbbPalette,
                     labels = c("Standard p-value", "sfFDR\nfunctional p-value", "Oracle\nfunctional p-value"))

ggsave(p1,
       file = "./figures/04-simulations-discoveries.png",
       dpi = 1000,
       width = 8,
       height = 4.5)

# Main Figure 2 in manuscript
p0 <- out %>%
  dplyr::filter(quantity == "fpvalues",
                prior.coverage == 0.025,
                threshold %in% c(5e-8), signal.density == "Medium",
                method != "sffdr_mono") %>%
  dplyr::group_by(prior.strength, threshold, method, signal.density, signal.density.z, prior.coverage) %>%
  ggplot(aes(x = factor(prior.strength),
             y = discoveries,
             color = method), group = signal.density.z) +
  geom_boxplot() +
  facet_grid(~ signal.density.z,
             scales = "free", labeller = labeller( signal.density.z = z.labs)) +
  theme_bw()  +
  xlab("") +
  ylab("Discoveries") +
  scale_color_manual("Method", values = cbbPalette,
                     labels = c("Standard p-value", "sfFDR\nfunctional p-value", "Oracle\nfunctional p-value"))

p2 <- outQ %>%
  dplyr::filter(quantity == "fqvalues",
                prior.strength != "None", method %in% c("sfFDR", "Oracle", "Q-value"),
                prior.coverage == 0.025,
                signal.density == "Medium",
                !(method %in% c("sffdr_mono" ))) %>%
  dplyr::group_by(prior.strength, threshold, method, signal.density, signal.density.z, prior.coverage) %>%
  dplyr::summarise(threshold = mean(threshold),
                   discoveries = mean(discoveries)) %>%
  ggplot(aes(x = threshold, y = discoveries, color = method, shape = factor(prior.strength)), group = signal.density.z) +
  geom_point() +
  geom_line() + #guides(colour = guide_legend(reverse=T)) +
  facet_grid(  ~ signal.density.z,
             scales = "free", labeller = labeller( signal.density.z = z.labs)) +
  theme_bw()  +
  #theme(legend.box.spacing = unit(0, "pt")) +
  scale_color_manual("Method",
                     breaks = c("Q-value", "sfFDR", "Oracle" ),
                     labels = c( "Standard q-value", "sfFDR\nfunctional q-value", "Oracle\nfunctional q-value"),
                     values = c("#D55E00",  "#009E73", "black")) +
  xlab("Target FDR") +guides(color = guide_legend(order = 1)) +
  ylab("Discoveries") +
  scale_shape_manual("Prior strength", values = c(15,16,17))

p <- (p2) / (p0) +
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") &  theme(legend.justification = "left", legend.box.spacing = unit(9, "pt"),  legend.key.spacing.y =  unit(2, "pt"))
ggsave(p,
       file = "./figures/04-simulations-FIG2.png",
       dpi = 1000,
       width = 8,
       height = 5)

############ Simulation study: linkage disequilibrium ############
# Simulation results with LD
fnames <- list.files("./simulation_study/data-LD/")

out <- list()
for (file in fnames) {
  print(file)
  df <- readRDS(paste0("./simulation_study/data-LD/", file))
  out <- rbind(out, df)
}

out$method <- factor(out$method, labels = c( "Oracle", "Raw P", "Functional P" ))
out$method <- fct_relevel(out$method,
                          "Raw P",
                          "Functional P",
                          "Oracle")

cbbPalette <- c("High" = "black", "Medium"= "darkgrey", "Low" = "lightgrey")

sd.labs <- c("High power: P", "Medium power: P", "Low power: P")
names(sd.labs) <- c("High", "Medium", "Low")

sd.labs <- c("High power: P", "Medium power: P", "Low power: P")
names(sd.labs) <- c("High", "Medium", "Low")

fnames <- c(TeX("Standard $p$-value"), TeX("sfFDR $p_f$-value"), TeX("Oracle $p_f$-value"))
names(fnames) <- c("Raw P", "Functional P", "Oracle")


fnames <- c(
  'Raw P' = "Standard p-value" , #@TeX("Standard $p$-value"),
  'Functional P' = "sfFDR functional p-value",#TeX("sfFDR $p_f$-value"),
  'Oracle' = "Oracle functional p-value"#TeX("Oracle $p_f$-value")
)
global_labeller <- labeller(
  signal.density = sd.labs,
  method = fnames
)

p1 <- out %>%
  dplyr::filter(quantity == "fpvalues",
         prior.coverage == 0.025,
         threshold %in% c(0.0001),
         method != "sffdr_mono") %>%
  dplyr::group_by(prior.strength, threshold, method, signal.density, signal.density.z, prior.coverage) %>%
  ggplot(aes(x = factor(prior.strength), y = empfdr,color = signal.density.z), group = signal.density.z) + geom_boxplot() +
  facet_grid(signal.density~method,
             scales = "free",labeller =  global_labeller) +
  geom_abline(slope = 1) + theme_bw()  +
  xlab("Prior strength") +
  ylab("False positive rate") +
  scale_color_manual("Z density", values = cbbPalette)  +
  geom_hline(data = data.frame(empfdr = rep(c( 0.0001),each = 1), threshold = rep(c( 0.0001),each = 1)), aes(yintercept = empfdr), linetype = "dashed")
ggsave(p1,
       file = "./figures/04-simulations-type1error_LD.png",
       dpi = 1000,
       width = 7,
       height = 4.5)

cbbPalette <- c("Oracle" = "black", "Functional P"= "#0072B2", "Raw P" = "#999999")
sd.labs <- c("High power: P", "Medium power: P", "Low power: P")
names(sd.labs) <- c("High", "Medium", "Low")
z.labs <- c("High power: Z", "Medium power: Z", "Low power: Z")
names(z.labs) <- c("High", "Medium", "Low")

p1 <- out %>%
  dplyr::filter(quantity == "fpvalues",
         prior.coverage == 0.025,
         threshold %in% c(5e-8),
         method != "sffdr_mono") %>%
  dplyr::group_by(prior.strength, threshold, method, signal.density, signal.density.z, prior.coverage) %>%
  ggplot(aes(x = factor(prior.strength),
             y = discoveries,
             color = method), group = signal.density.z) +
  geom_boxplot() +
  facet_grid(signal.density ~ signal.density.z,
             scales = "free", labeller = labeller(signal.density = sd.labs,signal.density.z = z.labs)) +
  theme_bw()  +
  xlab("") +
  ylab("Discoveries") +
   guides(color = guide_legend(order = 1)) +
  theme(legend.key.spacing.y =  unit(2, "pt")) +
  scale_color_manual("Method", values = cbbPalette,
                     labels = c("Standard p-value", "sfFDR\nfunctional p-value", "Oracle\nfunctional p-value"))

ggsave(p1,
       file = "./figures/04-simulations-discoveries-LD.png",
       dpi = 1000,
       width = 8,
       height = 4.5)

sd.labs <- c("High power: P", "Medium power: P", "Low power: P")
names(sd.labs) <- c("High", "Medium", "Low")
z.labs <- c("High power: Z", "Medium power: Z", "Low power: Z")
names(z.labs) <- c("High", "Medium", "Low")
cbbPalette <- c("High" = "black", "Medium"= "darkgrey", "Low" = "lightgrey")

p1 <- out %>%
  dplyr::filter(quantity == "fqvalues",
                prior.coverage == 0.025) %>%
  dplyr::group_by(prior.strength, threshold, method, signal.density, signal.density.z, prior.coverage) %>%
  dplyr::summarise(empfdr = mean(empfdr),
                   emptdr = mean(emptdr)) %>%
  ggplot(aes(x = empfdr, y = emptdr, color = method, shape = factor(prior.strength)), group = signal.density.z) +
  geom_point() +
  geom_line() +
  facet_grid(signal.density ~ signal.density.z,
             scales = "free", labeller = labeller(signal.density = sd.labs,signal.density.z = z.labs)) +
  theme_bw()  +
  scale_color_manual("Method",
                     label = c("Q-value", "sfFDR", "Oracle"),
                     breaks = c( "Raw P", "Functional P", "Oracle"),
                     values = c("#D55E00", "#009E73","black")) +
  xlab("Empirical false discovery rate") +
  ylab("Empirical true discovery rate")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_shape_manual("Prior strength", values = c(15,16,17))
ggsave(p1,
       file = "./figures/04-simulations-TDR-FDR-LD.png",
       dpi = 1000,
       width = 7,
       height = 4.5)

out$method <- factor(out$method, levels = c("oracle", "raw", "sffdr"), labels = c("Oracle", "Q-value",  "sfFDR"))
out$method <- fct_relevel(out$method,
                          "Q-value",
                          "Functional Q-value",
                          "Oracle")
p1 <- out %>%
  dplyr::filter(quantity == "fqvalues",
                prior.coverage == 0.025,
                threshold %in% c( 0.01),
                !(method %in% c("sffdr_mono" ))) %>%
  dplyr::group_by(prior.strength, threshold, method, signal.density, signal.density.z, prior.coverage) %>%
  ggplot(aes(x = as.factor(method), y = empfdr), group = signal.density.z) +
  geom_boxplot(outlier.size = 0.5) +
  facet_grid(signal.density ~ signal.density.z, scales = "free",
             labeller = labeller(signal.density = sd.labs,signal.density.z = z.labs)) +
  theme_bw()  +
  ylab("Empirical FDR") +
  xlab("FDR threshold") +# scale_y_log10() +
  scale_shape_manual("Prior strength", values = c(15,16,17)) +
  geom_hline(yintercept = 0.01, size = 0.7, color = "darkred", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")
ggsave(p1,
       file = "./figures/04-simulations-FDR-type1-error-LD.png",
       dpi = 1000,
       width = 7,
       height = 5)

############ EGPA ############
plot.sffdr <- function (x, rng = c(0, 1e-4), ...) {
  rm_na <- !is.na(x$pvalues)
  pvalues <- x$pvalues[rm_na]
  qvalues <- x$fqvalues[rm_na]
  fpvalues <- x$fpvalues[rm_na]

  fp.ord <- fpvalues[order(fpvalues)]
  if (min(fp.ord) > rng[2]) {
    rng <- c(min(fp.ord), quantile(fp.ord, 1e-2))
  }

  p.ord <- pvalues[fp.ord]

  pi0 <- x$fpi0

  pi0.df <- data.frame(z = rank(pi0, ties = "random")[which(pvalues >= rng[1] & pvalues <= rng[2])] / length(pi0),
                       pi0 = pi0[which(pvalues >= rng[1] & pvalues <= rng[2])])

  # surrogate versus pi0
  p1 <- pi0.df %>%
    ggplot(aes(x = z, y = pi0)) +
    #geom_point() +
    geom_line() +
    theme_bw() +
    scale_x_log10() +
    xlab("surrogate variable") +
    ylab("prior probability of null") +
    scale_color_brewer(palette = "Set1")

  # functional p versus raw p
  p2 <- ggplot(data.frame(pvalue = pvalues[which(pvalues >= rng[1] & pvalues <= rng[2])],
                          fpvalue = fpvalues[which(pvalues >= rng[1] & pvalues <= rng[2])]),
               aes(x =  pvalue, y = fpvalue)) +
    ylab("functional p-value") + scale_x_log10() + scale_y_log10() +
    xlab("p-value") + geom_point(size =1.2,alpha = 0.5) + theme_bw()
  rng = c(0, 5e-8)
  # Functional p versus discoveries
  df1 <- data.frame(threshold = fp.ord[fp.ord >= rng[1] & fp.ord <= rng[2]],
                    sig = (1 + sum(fp.ord < rng[1])):sum(fp.ord <=  rng[2]))
  df2 <- data.frame(threshold = sort(pvalues[pvalues >= rng[1] & pvalues <= rng[2]]),
                    sig = (1 + sum(pvalues < rng[1])):sum(pvalues <=  rng[2]))
  if (max(df2$threshold) > max(df1$threshold)) {
    df1 <- rbind(df1, data.frame(threshold = max(df2$threshold), sig = max(df1$sig)))
  } else {
    df2 <- rbind(df2, data.frame(threshold = max(df1$threshold), sig = max(df2$sig)))
  }
  if (min(df2$threshold) < min(df1$threshold)) {
    df1 <- rbind(df1, data.frame(threshold = min(df2$threshold), sig = min(df1$sig)))
  } else {
    df2 <- rbind(df2, data.frame(threshold = min(df1$threshold), sig = min(df2$sig)))
  }

  p3 <- ggplot(df1,
               aes(x = threshold, y = sig, linetype = "functional p-value")) +
    xlab("threshold") +
    ylab("significant tests") + geom_line( ) +
    geom_line(data = df2,
              aes(x=threshold, y= sig, linetype = "standard p-value") ) +  theme_bw() +
    scale_linetype_manual("", values = c("solid", "dashed")) +
    theme(strip.text.x = element_blank(),
          strip.background = element_rect(  fill=NA),
          legend.background=element_blank(),
        #  legend.key =
          legend.key = element_rect(colour = NA, fill = NA),
          legend.position =c(.245,.92),
          legend.key.spacing.y =  unit(-6, "pt"),
         legend.box.background = element_blank(),
           legend.text = element_text(size=7)
    )

  # functional p versus q
  p4 <-  ggplot(data.frame(fpvalue = fpvalues[which(fpvalues >= rng[1] & fpvalues <= rng[2])],
                           fqvalue =  qvalues[which(fpvalues >= rng[1] & fpvalues <= rng[2])]),
                aes(x = fpvalue, y =  fqvalue )) +
    xlab("functional p-value") + scale_x_log10() + scale_y_log10() +
    ylab("functional q-value") + #geom_point() +
    theme_bw() +
    geom_line() + theme_bw() + scale_x_log10()

  (p1+p2) / (p3+p4) + plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")")
}

flist <- list.files("./egpa_analysis/sffdr/")

df.sum <- NULL
for (file in flist) {
  df <- readRDS(paste0("./egpa_analysis/sffdr/", file))
  x <- data.frame(pvalues = df$P, fpvalues = df$fPi,  fqvalues = df$fQi, fpi0 = df$fPIi)
  class(x) <- "sffdr"
  p1 <- plot.sffdr(x)
  ggsave(plot = p1,
         filename = paste0("./figures/", file, ".png"),
         width = 6.5,
         dpi = 1000,
         height = 4.5)
  # P-values
  png(filename = paste0("./figures/manhattan-", file, "-rawP.png"),
      width     =  12, #9 ,#12
      height    = 4, #5.25,#4
      units     = "in",
      res       = 1500,)
  qqman::manhattan(df[df$fPi < 0.05, ], cex.axis=1.2, cex = 1.34, cex.lab = 1.3, suggestiveline =  FALSE, chr="CHR38", bp="BP38", p = "fPi", snp = "SNPID")
  dev.off()
}

# EGPA only
df <- readRDS("./egpa_analysis/sffdr/EGPA")
significanceset <- df %>% filter(fPi < 5e-8)
df.sum <- rbind(df.sum, data.frame(discoveries = c(sum(df$fPi<5e-8), sum(df$P<5e-8)),
                                   method = c("Functional P", "Raw P"),
                                   data = "EGPA"))
saveRDS(df.sum, file = "./discoveries.rds")
p1 <- data.frame(threshold = seq(0.00001, 0.01, 0.00001)) %>%
  group_by(threshold) %>%
  summarise(Discoveries = sum(df$fQi < threshold)) %>%
  ggplot(aes(x = threshold, y = Discoveries)) + geom_line() + theme_bw() + xlab("Threshold") + ylab("Discoveries")

cbbPalette <- c("Functional P" = "black", "Raw P"= "darkgrey", "Low" = "lightgrey")

p1 <- df.sum %>%
  ggplot(aes(x=(data) ,y = discoveries, fill = method, group = method)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_fill_manual("Method", values = cbbPalette) +
  xlab("") + ylab("Discoveries")

ggsave(p1,
       file = "./figures/01-egpa-discoveries.png",
       width = 6,
       height = 3.25)

# Results for locus zoom
dfzoom <- df %>%
  dplyr::ungroup() %>%
  dplyr::select(SNPID, fPi, ALT,BP38, CHR38) %>%
  dplyr::filter(fPi < 0.01)

dfzoom <- dfzoom %>%
  dplyr::group_by(CHR38) %>%
  dplyr::arrange(BP38, .by_group = TRUE)
write_tsv(dfzoom, file = "./locuszoom.tsv")

# Null applied
ret_df <- readRDS(file = "./egpa_analysis/functional-pvalue-null-inform.rds")
df <- ret_df %>%
  dplyr::group_by(id, trait) %>%
  dplyr::mutate(p = mean(p),
         fp = mean(fp),
         fq = mean(fq),
         fpi0 = mean(fpi0))

# Create null figure p versus functional p
p1 <- df %>%
  ggplot(aes(x = -log10(p), y = -log10(fp))) +
  geom_point(alpha = 0.05,
             size = 1) +
  theme_bw() +
  xlab("p-values") +
  ylab("Functional p-values") +
  facet_wrap(~trait ) +
  geom_abline(slope = 1) +
  coord_fixed()
ggsave(plot = p1,
       filename = paste0("./figures/pversusfp-applied.png"),
       width = 7,
       dpi = 1000,
       height = 2.5)

# Create null figure p versus functional p
p1 <- df %>%
  dplyr::filter(trait == "EGPA") %>%
  ggplot(aes(x = -log10(p), y = -log10(fp))) +
  geom_point(alpha = 0.05,
             size = 1) +
  theme_bw() +
  xlab("p-values") +
  ylab("Functional p-values") +
  geom_abline(slope = 1) +
  xlim(0, 9.5) + ylim(0,9.5) +
  coord_fixed()

ggsave(plot = p1,
       filename = paste0("./figures/pversusfp-applied-EGPA.png"),
       width = 3,
       dpi = 1000,
       height = 3)

ret_df <- readRDS(file = "./egpa_analysis/04-functional-pvalue-null-inform-uncorrelated.rds")

# Create null figure p versus functional p
p1 <- ret_df %>%
  dplyr::filter(trait == "EGPA") %>%
  ggplot(aes(x = -log10(p), y = -log10(fp))) +
  geom_point(alpha = 0.25,
             size = 1) +
  theme_bw() +
  xlab("p-values") +
  ylab("Functional p-values") +
  geom_abline(slope = 1) +
  xlim(0, 9.5) + ylim(0,9.5) +
  coord_fixed()

ggsave(plot = p1,
       filename = paste0("./figures/pversusfp-applied-EGPA-null-uncorrelated.png"),
       width = 3,
       dpi = 1000,
       height = 3)

########### FINE MAPPING w/ coloc package ################
# List of lead SNPs
df <- readRDS(paste0("./egpa_analysis/sffdr/EGPA"))
chr <- c(5, 2, 5, 10, 6, 21, 3, 17, 12, 11)
bp <- c(111066174, 111150990, 132454924, 9010073,
        90116186, 35365944,188730411, 49370984,
        56007301, 76590331)
gene <- c("TSLP", "BCL2L11", "IRF1", "GATA3", "BACH2",
          "RUNX1",  "TPRG1",  "ZNF652", "IKZF4", "LRRC32")
sdf <- data.frame(CHR38 = chr,
                  BP38 = bp, gene = gene)

sdf <- data.frame(CHR38 = chr,
                  BP38 = bp)
out <- NULL
for (i in 1:length(chr)) {
  dat <- df %>%
    dplyr::filter(CHR38 == chr[i],
                  BP38 >= bp[i]-500000,
                  BP38 <= bp[i]+500000)
  dat$gene = gene[i]
  out <- rbind(out, dat)
}
theme_update(plot.title = element_text(hjust = 0.5))
rm(df)

col_select <- c("SNPID", "INFO", "CHR38", "BP38", "ALT", "REF", "ALTFREQ",   "BETA", "SE")
tmp <- read_tsv(paste0("./egpa_analysis/primary/EGPA_Lyons_31719529_1-hg38.tsv.gz"),  col_select = col_select)
outc <- out %>% left_join(tmp, by = c("CHR38", "BP38"))

ss <- outc[, c("gene", "BETA", "SE", "BP38", "P", "ALTFREQ", "SNPID.x")]
colnames(ss) <- c("gene", "beta", "varbeta", "BP38", "pvalues", "MAF", "snp")

genes <- ss %>% group_by(gene) %>% filter(sum(pvalues < 5e-8) > 0)

#  EGPA study: CI for two genes
ss$type = "cc"
ss$s <- 676
tmp <- ss %>% filter(gene == "TSLP") #=="BCL2L11"
tmp <- list(
  pvalues = tmp$pvalues,
  MAF= tmp$MAF,
  snp = tmp$snp,
  type = "cc",
  N = (6809 + 676),
  s = 676 / (6809 + 676))

tt <- coloc::finemap.abf(tmp)
PP <- tt$SNP.PP[-length(tt$lABF)]
PP <- PP / sum(PP)
o <- order(PP,decreasing=TRUE)
cs <- cumsum(PP[o])
w <- which(cs > 0.95)[1]

# Get list of genes with credible sets
col_select <- c("CHR38", "BP38", "ALT_FREQ", "P", "SE")
tmp <- read_tsv(paste0("./egpa_analysis/conditional/ASTAO_Ferreira_30929738_1-hg38.tsv.gz"),  col_select = col_select)
tmp$CHR38 = as.integer(tmp$CHR38)
outc <- out %>% left_join(tmp, by = c("CHR38", "BP38"))
rm(tmp)
ss <- outc[, c("gene",  "SNPID", "SE", "BP38", "CHR38", "P.y", "ALT_FREQ" )]
colnames(ss) <- c("gene", "snp",  "varbeta", "BP38","CHR38", "pvalues", "MAF"  )

genes <- ss %>% group_by(gene) %>% filter(sum(pvalues < 5e-8) > 0)

glist <- c("TSLP", "BCL2L11", "IRF1", "GATA3", "BACH2",
           "RUNX1",  "TPRG1",  "ZNF652", "IKZF4", "LRRC32")
library(coloc)
df <- NULL
for (g in glist) {
  ss$type = "cc"
  ss$s <- 26582
  tmp0 <- ss %>% filter(gene == g)
  tmp <- list(
    pvalues = tmp0$pvalues,
    MAF= tmp0$MAF,
    snp = tmp0$snp,
    type = "cc",
    N = (26582 + 300671),
    s = 26582 / (26582 + 300671))

  tt <- coloc::finemap.abf(tmp)
  tt$CHR38 <- c(tmp0$CHR38, NA)
  tt$BP38 <- c(tmp0$BP38, NA)
  PP <- tt$SNP.PP[-length(tt$lABF)]
  PP <- PP / sum(PP)
  o <- order(PP, decreasing = TRUE)
  ro <- order(o)
  cs <- cumsum(PP[o])
  w <- which(cs > 0.95)[1]
  cs <- cs[ro]
  tts <- tt[-nrow(tt),]
  tts$PP <- cs
  subset <- tts %>% top_n(n=w, wt = -PP)
  subset$gene <- g
  subset$data <- "ASTAO"
  df <- rbind(subset, df)
}

# ASTCO
col_select <- c(  "CHR38", "BP38",   "ALT_FREQ", "P",    "SE")
tmp <- read_tsv(paste0("./egpa_analysis/conditional/ASTCO_Ferreira_30929738_1-hg38.tsv.gz"),  col_select = col_select)
tmp$CHR38 = as.integer(tmp$CHR38)
outc <- out %>% left_join(tmp, by = c("CHR38", "BP38"))
rm(tmp)
ss <- outc[, c("gene", "SE", "BP38","CHR38", "P.y", "ALT_FREQ", "SNPID")]
colnames(ss) <- c("gene", "varbeta", "BP38", "CHR38", "pvalues", "MAF", "snp")
for (g in glist) {
  ss$type = "cc"
  ss$s <- 13962
  tmp0 <- ss %>% filter(gene == g)
  tmp <- list(
    pvalues = tmp0$pvalues,
    MAF= tmp0$MAF,
    snp = tmp0$snp,
    type = "cc",
    N = (13962 + 300671),
    s = 13962 / (13962 + 300671))

  tt <- coloc::finemap.abf(tmp)
  tt$CHR38 <- c(tmp0$CHR38, NA)
  tt$BP38 <- c(tmp0$BP38, NA)
  PP <- tt$SNP.PP[-length(tt$lABF)]
  PP <- PP / sum(PP)
  o <- order(PP, decreasing = TRUE)
  ro <- order(o)
  cs <- cumsum(PP[o])
  w <- which(cs > 0.95)[1]
  cs <- cs[ro]
  tts <- tt[-nrow(tt),]
  tts$PP <- cs
  subset <- tts %>% top_n(n=w, wt = -PP)
  subset$gene <- g
  subset$data <- "ASTCO"
  df <- rbind(subset, df)
}

# EOSC
col_select <- c(  "CHR38", "BP38", "ma_freq",  "P", "SE")
tmp <- read_tsv(paste0("./egpa_analysis/conditional/EOSC_Astle_27863252_1-hg38.tsv.gz"),  col_select = col_select)
tmp$CHR38 = as.integer(tmp$CHR38)
outc <- out %>% left_join(tmp, by = c("CHR38", "BP38"))

ss <- outc[, c("gene",   "SE", "BP38", "CHR38", "P.y", "ma_freq", "SNPID")]
colnames(ss) <- c("gene",  "varbeta", "BP38", "CHR38", "pvalues", "MAF", "snp")
for (g in glist) {
  ss$type = "quant"
  tmp0 <- ss %>% filter(gene == g)
  tmp <- list(
    pvalues = tmp0$pvalues,
    MAF = tmp0$MAF,
    snp = tmp0$snp,
    varbeta = tmp0$varbeta,
    type = "quant",
    N = 172275)

  tt <- coloc::finemap.abf(tmp)
  tt$CHR38 <- c(tmp0$CHR38, NA)
  tt$BP38 <- c(tmp0$BP38, NA)
  PP <- tt$SNP.PP[-length(tt$lABF)]
  PP <- PP / sum(PP)
  o <- order(PP, decreasing = TRUE)
  ro <- order(o)
  cs <- cumsum(PP[o])
  w <- which(cs > 0.95)[1]
  cs <- cs[ro]
  tts <- tt[-nrow(tt),]
  tts$PP <- cs
  subset <- tts %>% top_n(n=w, wt = -PP)
  subset$gene <- g
  subset$data <- "EOSC"
  df <- rbind(subset, df)
}

saveRDS(df, file = "./overlap.rds")

########### FINE MAPPING SFFDR ################
# List of lead SNPs
df <- readRDS(paste0("./egpa_analysis/sffdr/EGPA"))
chr <- c(5, 2, 5, 10, 6, 21, 3, 17, 12, 11)
bp <- c(111066174, 111150990, 132454924, 9010073,
        90116186, 35365944,188730411, 49370984,
        56007301, 76590331)
gene <- c("TSLP", "BCL2L11", "IRF1", "GATA3", "BACH2",
          "RUNX1",  "TPRG1",  "ZNF652", "IKZF4", "LRRC32")
sdf <- data.frame(CHR38 = chr,
                  BP38 = bp)
out <- NULL
for (i in 1:length(chr)) {
  dat <- df %>%
    dplyr::filter(CHR38 == chr[i],
                  BP38 >= bp[i]-500000,
                  BP38 <= bp[i]+500000)
  dat$gene = gene[i]
  out <- rbind(out, dat)
}
theme_update(plot.title = element_text(hjust = 0.5))

# Fine mapping results w/ coloc
cdata <- readRDS("overlap.rds")
overlap <- cdata %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHR38,BP38) %>%
  dplyr::mutate(shape = ifelse(data== "ASTAO", 15, ifelse(data == "ASTCO", 17, 18))) %>%
  dplyr::select(CHR38, BP38, gene, shape)

p0 <- out %>%
  dplyr::ungroup() %>%
  dplyr::mutate(BF = (mean(fPIi) / (1-mean(fPIi))) * ((1-flfdri) / flfdri)) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(PIP = BF / sum(BF)) %>%
  dplyr::arrange(-1*PIP) %>%
  dplyr::mutate(CD = cumsum(PIP)) %>%
  dplyr::left_join(overlap) %>%
  dplyr::group_by(BP38, CHR38, gene) %>%
  dplyr::mutate(shape = ifelse(is.na(shape), 16, shape)) %>%
  dplyr::ungroup() %>%
  group_by(gene) %>%
  dplyr::mutate(color = c(rep(TRUE, length(CD[CD < 0.95]) + 1), rep(FALSE, length(CD[CD > 0.95])-1))) %>%
  dplyr::group_by(CHR38, BP38) %>%
  dplyr::mutate(color= ifelse(gene == "TSLP", ifelse(PIP >0.9, TRUE, FALSE), color)) %>%
  dplyr::ungroup() %>%
  ggplot(aes(x = BP38, y = PIP, color = color, shape  =  as.factor(shape))) +
  geom_point(size =1.5, alpha = 0.7) +
  xlab("Position") +
  facet_wrap(~gene,nrow = 1, scales = "free") +
  ylab("PP") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_shape_manual("CS Overlap", labels =  c("None", "ASTAO", "ASTCO", "EOSC"),  values = c(16,15,17,18), breaks = c(16,15,17,18)) +
  scale_colour_manual(guide = 'none', values=c("#000000", "darkred")) +theme(plot.title = element_text(hjust = 0.5)) #+ ylim(0,1)

labs <-  c("-log10(func. p)","-log10(p)", "PP")
names(labs) <- c( "fPi","P", "PIP")
c("TSLP", "BCL2L11", "IRF1", "GATA3", "BACH2",
  "RUNX1",  "TPRG1",  "ZNF652", "IKZF4", "LRRC32")
p1<-out %>% filter(gene %in% c("BACH2", "BCL2L11", "GATA3", "IRF1", "IKZF4")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(BF = (mean(fPIi) / (1-mean(fPIi))) * ((1-flfdri) / flfdri)) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(PIP = BF / sum(BF)) %>%
  dplyr::arrange(-1*PIP) %>%
  dplyr::mutate(CD = cumsum(PIP)) %>%
  dplyr::left_join(overlap) %>%
  dplyr::group_by(BP38, CHR38, gene) %>%
  dplyr::mutate(shape = ifelse(is.na(shape), 16, shape)) %>%
  dplyr::ungroup() %>%
  group_by(gene) %>%
  dplyr::mutate(color = c(rep(TRUE, length(CD[CD < 0.95]) + 1), rep(FALSE, length(CD[CD > 0.95])-1))) %>%

  dplyr::group_by(CHR38, BP38) %>%
  dplyr::mutate(color= ifelse(gene == "TSLP", ifelse(PIP >0.9, TRUE, FALSE), color)) %>%
  dplyr::ungroup() %>%
  mutate(fPi = -log10(fPi)) %>%
  mutate(P = -log10(P)) %>%
  pivot_longer(cols = c("fPi","P", "PIP"), names_to = "stat", values_to = "value") %>%
  mutate(stat = factor(stat, levels = c("P", "fPi", "PIP"))) %>%
  ggplot(aes(x = BP38, y = value, color = color, shape  =  as.factor(shape))) +
  geom_point(size =3.5, alpha = 0.7) +
  xlab("Position") +
  ggh4x::facet_grid2(stat~gene,independent = "y",scales = "free",  labeller = labeller(stat = labs)) +
  ylab("Value") +
  theme_bw(base_size = 23) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_shape_manual("CS Overlap", labels =  c("None", "ASTAO", "ASTCO", "EOSC"),  values = c(16,15,17,18), breaks = c(16,15,17,18)) +
  scale_colour_manual(guide = 'none', values=c("#000000", "darkred")) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.text = element_text(margin = margin(t=7,b=7, l=7, unit = "pt")))

p2 <- out %>% filter(!(gene %in%  c("BACH2", "BCL2L11", "GATA3", "IRF1", "IKZF4"))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(BF = (mean(fPIi) / (1-mean(fPIi))) * ((1-flfdri) / flfdri)) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(PIP = BF / sum(BF)) %>%
  dplyr::arrange(-1*PIP) %>%
  dplyr::mutate(CD = cumsum(PIP)) %>%
  dplyr::left_join(overlap) %>%
  dplyr::group_by(BP38, CHR38, gene) %>%
  dplyr::mutate(shape = ifelse(is.na(shape), 16, shape)) %>%
  dplyr::ungroup() %>%
  group_by(gene) %>%
  dplyr::mutate(color = c(rep(TRUE, length(CD[CD < 0.95]) + 1), rep(FALSE, length(CD[CD > 0.95])-1))) %>%
  dplyr::group_by(CHR38, BP38) %>%
  dplyr::mutate(color= ifelse(gene == "TSLP", ifelse(PIP >0.9, TRUE, FALSE), color)) %>%
  dplyr::ungroup() %>%
  mutate(fPi = -log10(fPi)) %>%
  mutate(P = -log10(P)) %>%
  pivot_longer(cols = c("fPi","P", "PIP"), names_to = "stat", values_to = "value") %>%
  mutate(stat = factor(stat, levels = c("P", "fPi", "PIP"))) %>%
  ggplot(aes(x = BP38, y = value, color = color, shape  =  as.factor(shape))) +
  geom_point(size =3.5, alpha = 0.7) +
  xlab("Position") +
  ggh4x::facet_grid2(stat~gene,independent = "y",scales = "free",  labeller = labeller(stat = labs)) +
  ylab("Value") +
  theme_bw(base_size = 23) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_shape_manual("CS Overlap", labels =  c("None", "ASTAO", "ASTCO", "EOSC"),  values = c(16,15,17,18), breaks = c(16,15,17,18)) +
  scale_colour_manual(guide = 'none', values=c("#000000", "darkred")) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.text = element_text(margin = margin(t=7,b=7, l=7, unit = "pt")))

p <- p1 / p2  + plot_layout(guides = "collect")
ggsave(p,
       file = "./figures/finemapping2.png",
       width = 17,
       dpi = 1500,
       height = 14)

tab <- out %>%
  dplyr::ungroup() %>%
  dplyr::mutate(BF = (mean(fPIi) / (1-mean(fPIi))) * ((1-flfdri) / flfdri)) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(PIP = BF / sum(BF)) %>%
  dplyr::arrange(-1*PIP) %>%
  dplyr::mutate(CD = cumsum(PIP)) %>%
  dplyr::left_join(overlap) %>%
  dplyr::group_by(BP38, CHR38, gene) %>%
  dplyr::mutate(shape = ifelse(is.na(shape), 16, shape)) %>%
  dplyr::ungroup() %>%
  group_by(gene) %>%
  dplyr::mutate(color = c(rep(TRUE, length(CD[CD < 0.95]) + 1), rep(FALSE, length(CD[CD > 0.95])-1))) %>%
  dplyr::group_by(CHR38, BP38) %>%
  dplyr::mutate(color= ifelse(gene == "TSLP", ifelse(PIP >0.9, TRUE, FALSE), color)) %>%
  dplyr::ungroup() %>% filter(color == TRUE) %>% dplyr::select(CHR38, BP38, ALT, SNPID, gene, P, ASTAO, ASTCO, EOSC, flfdri, fPIi, BF, PIP) %>%
  distinct()  %>% dplyr::rename(pi0 = fPIi,
                         lfdr = flfdri)
out %>%
    dplyr::ungroup() %>%
    dplyr::mutate(BF = (mean(fPIi) / (1-mean(fPIi))) * ((1-flfdri) / flfdri)) %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(PIP = BF / sum(BF)) %>%
    dplyr::arrange(-1*PIP) %>%
    dplyr::mutate(CD = cumsum(PIP)) %>%
    dplyr::left_join(overlap) %>%
    dplyr::group_by(BP38, CHR38, gene) %>%
    dplyr::mutate(shape = ifelse(is.na(shape), 16, shape)) %>%
    dplyr::ungroup() %>%
    group_by(gene) %>%
    dplyr::mutate(color = c(rep(TRUE, length(CD[CD < 0.95]) + 1), rep(FALSE, length(CD[CD > 0.95])-1))) %>%
    dplyr::group_by(CHR38, BP38) %>%
    dplyr::mutate(color= ifelse(gene == "TSLP", ifelse(PIP >0.9, TRUE, FALSE), color)) %>%
    dplyr::ungroup() %>% filter(color == TRUE) %>% group_by(gene) %>% dplyr::select(CHR38, BP38, ALT, SNPID, gene, P, ASTAO, ASTCO, EOSC, flfdri, fPIi, BF, PIP) %>% distinct() %>%  summarise(length(SNPID))
out %>%
  dplyr::ungroup() %>%
  dplyr::mutate(BF = (mean(fPIi) / (1-mean(fPIi))) * ((1-flfdri) / flfdri)) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(PIP = BF / sum(BF)) %>%
  dplyr::arrange(-1*PIP) %>%
  dplyr::mutate(CD = cumsum(PIP)) %>%
  dplyr::left_join(overlap) %>%
  dplyr::group_by(BP38, CHR38, gene) %>%
  dplyr::mutate(shape = ifelse(is.na(shape), 16, shape)) %>%
  dplyr::ungroup() %>%
  group_by(gene) %>%
  dplyr::mutate(color = c(rep(TRUE, length(CD[CD < 0.95]) + 1), rep(FALSE, length(CD[CD > 0.95])-1))) %>%
  dplyr::group_by(CHR38, BP38) %>%
  dplyr::mutate(color= ifelse(gene == "TSLP", ifelse(PIP >0.9, TRUE, FALSE), color)) %>%
  dplyr::ungroup() %>% filter(color == TRUE) %>% group_by(gene) %>% dplyr::select(CHR38, BP38, ALT, SNPID, gene, P, ASTAO, ASTCO, EOSC, flfdri, fPIi, BF, PIP) %>% distinct() %>%  summarise(length(SNPID))
saveRDS(tab, file = "EGPA_CI.rds")

p_visual <- out %>% filter(gene == "RUNX1") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(BF = (mean(fPIi) / (1-mean(fPIi))) * ((1-flfdri) / flfdri)) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(PIP = BF / sum(BF)) %>%
  dplyr::arrange(-1*PIP) %>%
  dplyr::mutate(CD = cumsum(PIP)) %>%
  dplyr::left_join(overlap) %>%
  dplyr::group_by(BP38, CHR38, gene) %>%
  dplyr::mutate(shape = ifelse(is.na(shape), 16, shape)) %>%
  dplyr::ungroup() %>%
  group_by(gene) %>%
  dplyr::mutate(color = c(rep(TRUE, length(CD[CD < 0.95]) + 1), rep(FALSE, length(CD[CD > 0.95])-1))) %>%
  dplyr::group_by(CHR38, BP38) %>%
  dplyr::ungroup() %>%
  ggplot(aes(x = BP38, y = -log10(fPi))) +
  geom_point(size =1.5, alpha = 0.7) +
  xlab("Position") +
  ylab("PP") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# CI overlap with summary statistics
dt <- out %>%
  dplyr::ungroup() %>%
  dplyr::mutate(BF = (mean(fPIi) / (1-mean(fPIi))) * ((1-flfdri) / flfdri)) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(PIP = BF / sum(BF)) %>%
  dplyr::arrange(-1*PIP) %>%
  dplyr::mutate(CD = cumsum(PIP)) %>%
  dplyr::left_join(overlap) %>%
  dplyr::group_by(BP38, CHR38, gene) %>%
  dplyr::mutate(shape = ifelse(is.na(shape), 16, shape)) %>%
  dplyr::ungroup() %>%
  group_by(gene) %>%
  dplyr::mutate(color = c(rep(TRUE, length(CD[CD < 0.95]) + 1), rep(FALSE, length(CD[CD > 0.95])-1))) %>%
  group_by(CHR38, BP38) %>%
  dplyr::mutate(color= ifelse(gene == "TSLP", ifelse(PIP >0.9, TRUE, FALSE), color)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(color == TRUE) %>%
  dplyr::ungroup() %>% dplyr::group_by(gene, shape) %>%
  dplyr::summarise(ll = sum(PIP[shape != 16])) %>%spread(shape,ll)
colnames(dt) <- c("Gene", "ASTAO", "NA", "ASTCO", "EOSC")
xtable(dt[,-3], digits = 3, type = "latex")

######## Statistics for main table in paper #####
df <- readRDS(paste0("./egpa_analysis/sffdr/EGPA"))
chr <- c(5, 2, 5, 10, 6, 21, 3, 17, 12, 11)
bp <- c(111066174, 111150990, 132454924, 9010073,
        90116186, 35365944,188730411, 49370984,
        56007301, 76590331)
gene <- c("TSLP", "BCL2L11", "IRF1", "GATA3", "BACH2",
          "RUNX1",  "TPRG1",  "ZNF652", "IKZF4", "LRRC32")
sdf <- data.frame(CHR38 = chr,
                  BP38 = bp, gene = gene)
tab <- sdf %>% left_join(df) %>% dplyr::select(CHR38, BP38, gene, SNPID, P, fQi)
files <- list.files("./egpa_analysis/primary/")[1]

# Merge with original effect sizes
col_select <- c("SNPID", "INFO", "CHR38", "BP38", "ALT", "REF", "ALTFREQ","P",  "OR",  "BETA", "SE")
tmp <- read_tsv(paste0("./egpa_analysis/primary/EGPA_Lyons_31719529_1-hg38.tsv.gz"),  col_select = col_select)
tmp$CHR38 <- as.numeric(tmp$CHR38)
tmp$BP38 <- as.numeric(tmp$BP38)
tab <- tab %>% left_join(tmp %>% dplyr::select(CHR38, BP38,ALT,REF, BETA,P, OR)%>% dplyr::rename( EGPA_REF = REF, EGPA_P = P, EGPA_OR =OR,EGPA_ALT = ALT, EGPA_BETA = BETA))

col_select <- c("SNPID", "INFO", "CHR38", "BP38", "ALT", "REF",  "OR",  "P",  "BETA", "SE")
tmp <- read_tsv(paste0("./egpa_analysis/conditional/ASTAO_Ferreira_30929738_1-hg38.tsv.gz"),  col_select = col_select)
tmp$CHR38 <- as.numeric(tmp$CHR38)
tmp$BP38 <- as.numeric(tmp$BP38)
tab <- tab %>% left_join(tmp %>% dplyr::select(CHR38, BP38,ALT,REF,P, OR, BETA)%>% dplyr::rename(ASTAO_REF = REF, ASTAO_P = P, ASTAO_OR = OR, ASTAO_ALT = ALT, ASTAO_BETA = BETA))

col_select <- c("SNPID", "INFO", "CHR38", "BP38", "ALT", "REF", "P", "OR",  "BETA", "SE")
tmp <- read_tsv(paste0("./egpa_analysis/conditional/ASTCO_Ferreira_30929738_1-hg38.tsv.gz"),  col_select = col_select)
tmp$CHR38 <- as.numeric(tmp$CHR38)
tmp$BP38 <- as.numeric(tmp$BP38)
tab <- tab %>% left_join(tmp %>% dplyr::select(CHR38, BP38,ALT,REF,P, OR, BETA)%>% dplyr::rename(ASTCO_REF = REF, ASTCO_ALT = ALT, ASTCO_P = P, ASTCO_OR = OR, ASTCO_BETA = BETA))

col_select <- c("SNPID",   "CHR38", "BP38", "ALT", "REF", "direction", "P",  "BETA", "SE")
tmp <- read_tsv(paste0("./egpa_analysis/conditional/EOSC_Astle_27863252_1-hg38.tsv.gz"),  col_select = col_select)
tmp$CHR38 <- as.numeric(tmp$CHR38)
tmp$BP38 <- as.numeric(tmp$BP38)
tab <- tab %>% left_join(tmp %>% dplyr::select(CHR38, BP38,ALT,REF, BETA, P, direction)%>% dplyr::rename(EOSC_REF = REF, EOSC_ALT = ALT,  EOSC_P = P, EOSC_BETA = BETA))

######## Compare EGPA p-values to published #####
col_select <- c("SNPID", "INFO", "CHR38", "BP38", "ALT", "REF", "ALTFREQ",   "BETA", "SE", "P")
tmp <- read_tsv(paste0("./egpa_analysis/primary/EGPA_Lyons_31719529_1-hg38.tsv.gz"), col_select = col_select)

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
  mutate(Pnew = 2*pnorm(abs(BETA / SE), lower.tail = F)) %>%
  select(-BETA, -SE)

p <- tmp %>% filter(!is.na(CHR38)) %>%
  ggplot(aes(x = -log10(Pnew), y = -log10(P))) +
  geom_point(size = 0.05, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + xlab("-log10(P) [Bolt LMM]") + ylab("-log10(P) [Published]")# +
ggsave(p,
       file = "./figures/PubPversusTransP.png",
       width = 4,
       dpi = 1000,
       height = 4)

