rm(list = ls())
library(VariantAnnotation)
library(updog)
library(tidyverse)
library(ldsep)
library(ggplot2)
library(ASRgenomics)
library(caret)
library(ggfortify)

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/1_updog/")
D1 <- read.table("DAl21_poly.txt", sep = "\t", row.names = 1, header = T)
D1[1:5,1:5] 
dim(D1)
# 1502 2477 DAl21_poly.txt
# 1479 4452 DAl21_poly_ori.txt
# 1502 2953 DAl21_poly1.txt


# class(D1)
# NZV <- nearZeroVar(D1)
# D1 <- D1[,-NZV]

D1 <- as.matrix(D1)

# ~~~~~~~~~~~
setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/1_updog/")
load("mout_2.Rdata")

summary(mout_2[["snpdf"]][["bias"]])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.01379  0.57766  0.91835  1.16047  1.14082 69.67090 

summary(mout_2[["snpdf"]][["od"]])
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0002674 0.0042218 0.0079175 0.0180830 0.0181405 0.5195386 

summary(mout_2[["snpdf"]][["prop_mis"]])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.01026 0.02480 0.05636 0.06948 0.71452 
hist(mout_2[["snpdf"]][["prop_mis"]])

# mout2 <- filter_snp(mout_2, bias < 1.14082 & od < 0.0181405 & pmax(Pr_0, Pr_1, Pr_2, Pr_3, Pr_4) < 0.95 & prop_mis < 0.1) # 1229 1502
# 
# mout2 <- filter_snp(mout_2, bias < 1.14082 & od < 0.0181405 & pmax(Pr_0, Pr_1, Pr_2, Pr_3, Pr_4) < 0.95) # 1387 1502

genomat1 <- format_multidog(mout_2, varname = "geno")
dim(genomat1) # 2953 1502
genomat1[1:5,1:7]

# snp.pruning -------------------------------------------------------------
setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/1_updog/")
c1 <- read.table("DAl21_poly.txt", sep = "\t", header = T, check.names = F, row.names = 1)
c1[1:5,1:5]
dim(c1) # 1502 2477

c1 <- as.matrix(c1)
geno.ps <- c1 / 4
geno.ps[1:5,1:5]
Mpr <- snp.pruning(M = geno.ps, pruning.thr = 0.9, window.n = 100, by.chrom = F, overlap.n = 10, seed = 1208, iterations = 15)

genomat5 <- Mpr$Mpruned * 4
genomat5[1:5,1:5]
dim(genomat5) # 1502 2455

NZV <- nearZeroVar(genomat5)
genomat5 <- genomat5[,-NZV]
dim(genomat5) # 1502 2434

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/1_updog/")
write.csv(genomat5, "DAl21_pruned.txt", row.names = T, quote = F)


# comboInfo <- findLinearCombos(genomat5)
# genomat7 <- genomat5[, -comboInfo$remove]


# df2 = cor(genomat5)
# hc = findCorrelation(df2, cutoff=0.9)
# hc = sort(hc)
# genomat7 = genomat5[,-c(hc)]


genomat5 <- as.data.frame(genomat5)
rownames(genomat5)


lev0 <- subset(rownames(genomat5), grepl("^UMN3097", rownames(genomat5)))
genomat6 <- genomat5 %>% dplyr::filter(!rownames(genomat5) %in% lev0) # no C0

# GWASpoly format ---------------------------------------------------------

genomat3 <- t(genomat5)
genomat3[1:5,1:7]
marks.4 <- as.data.frame(genomat3)
marks.4 <- marks.4 %>% tibble::rownames_to_column("snp_name") %>% separate(col = 1, into = c("Chrom", "Position"), remove = T, sep = "_") 
marks.4[1:5,1:7]

Marker <- seq(1:nrow(marks.4))
G6 <- cbind(Marker, marks.4)
G6$Chrom <- as.factor(G6$Chrom)

G6$Chrom <- recode_factor(G6$Chrom,
                          chr1.1 = "1", chr2.1 = "2", chr3.1 = "3",
                          chr4.1 = "4", chr5.1 = "5", chr6.1 = "6",
                          chr7.1 = "7", chr8.1 = "8")
summary(G6$Chrom)
length(Marker)
dim(G6) # 2554 1482
G6[1:5,1:7]

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/3_GWASpoly/")
write.csv(G6, "DAl21_GWAS.txt", row.names = F, quote = F)

# pca ---------------------------------------------------------------------

# genomat5[1:5,1:5]
# genomat6 <- as.data.frame(genomat5)
# geno.scale <- scale(genomat5, center = T, scale = T)
# svdgeno <- svd(geno.scale)
# PCA <- geno.scale %*% svdgeno$v
# PCA[1:5,1:5]
# 
# PCA4 <- PCA[,c(1:5)]
# head(PCA4)
# class(PCA4)
# colnames(PCA4) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
# PCA4 <- as.data.frame(PCA4)
# PCA4 <- PCA4 %>% rownames_to_column(var = "gen")
# 
# setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/2_pheno/2_BLUPs/")
# write.csv(PCA4, "pca_1.csv", quote = F, row.names = F)
# 
# geno.scale1 <- prcomp(genomat5, center = T, scale = T)
# attributes(geno.scale1)
# print(geno.scale1)
# 
# # Proportion of variance explained by PCA1, PCA2, PCA3
# PCA1 <- 100*round((svdgeno$d[1])^2 / sum((svdgeno$d)^2), d = 3); PCA1
# PCA2 <- 100*round((svdgeno$d[2])^2 / sum((svdgeno$d)^2), d = 3); PCA2
# PCA3 <- 100*round((svdgeno$d[3])^2 / sum((svdgeno$d)^2), d = 3); PCA3
# 
# # screenplot to visualize the proportion of variance explained by PCA
# plot(round((svdgeno$d)^2 / sum((svdgeno$d)^2), d = 7)[1:10], type = 'o', xlab = 'PCA', ylab = '% variance')

# PCA
library(ggfortify)
library(RColorBrewer) 

class(D1)
D1[1:5,1:5]
genomat5 <- as.data.frame(D1)
genomat5 <- genomat5 %>% rownames_to_column("Taxa")
genomat5[1:5,1:5]
genomat5$Taxa

# %>% dplyr::filter(Taxa %in% c("ZG21", "ZG25"))

# marks6 <- marks.5 %>% dplyr::filter(Taxa %in% c("ZG20","ZG21","ZG23","ZG25","ZG9"))

genomat5 <- separate(genomat5, 1, c("pop", "id"), sep = "_", remove = F, convert = FALSE, extra = "warn")
genomat5$pop <- as.factor(genomat5$pop)
summary(genomat5$pop)
genomat5[1:5,1:5]

genomat5$pop <- recode_factor(genomat5$pop,
                              UMN3097 = "Start_0", UMN3355 = "HighDig_1",
                              UMN3358 = "LowDig_1", UMN4016 = "HighDig_2",
                              UMN4351 = "LowDig_2")

genomat5 <- separate(genomat5, 2, c("Direction", "Cycle"), sep = "_", remove = T, convert = FALSE, extra = "warn")
genomat5$Cycle <- as.factor(genomat5$Cycle)
genomat5$Direction <- as.factor(genomat5$Direction)
summary(genomat5$Cycle)
summary(genomat5$Direction)

marks.6 <- genomat5[,-c(1:4)]
marks.6[1:5,1:5]
class(marks.6)
pca_res <- prcomp(marks.6, scale. = F, center = T)

autoplot(pca_res)

pc1 <- autoplot(pca_res, data = genomat5, colour = 'Direction', shape = "Cycle", alpha = 0.5) + theme_bw(base_family = "Arial", base_size = 12) + scale_colour_brewer(palette = "Dark2") + theme(panel.grid.minor = element_blank())


ggsave(filename = "~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/4_plots/PCA.pdf", plot = pc1, dpi = 300, width = 6, height = 4, device = cairo_pdf)

pc2 <- plot(pca_res, type="l", title(main = NULL))

ggsave(filename = "~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/4_plots/PCA.pdf", plot = pc1, dpi = 300, width = 6, height = 4, device = cairo_pdf)


# lev1 <- unique(subset(PCA$Taxa,  grepl("ZG", PCA$Taxa)))
PCA$Taxa <- gsub("ZG", "ZG_", PCA$Taxa)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~