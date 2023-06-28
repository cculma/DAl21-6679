rm(list = ls())
library(VariantAnnotation)
library(updog)
library(tidyverse)
library(ldsep)
library(ggplot2)
library(ASRgenomics)
library(caret)
library(ggfortify)
library(factoextra)
library(FactoMineR)

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/1_updog/")
D1 <- read.table("DAl21_poly.txt", sep = "\t", row.names = 1, header = T)
D1[1:5,1:5]
dim(D1)
D1 <- as.matrix(D1)

class(D1)
D1[1:5,1:5]
genomat5 <- as.data.frame(D1)
genomat5 <- genomat5 %>% rownames_to_column("gen")
genomat5[1:5,1:5]

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
genomat5[1:5,1:5]


marks.6 <- genomat5[,-c(1:4)]
marks.6[1:5,1:5]
dim(marks.6)
class(marks.6)
pca_res <- prcomp(marks.6, scale. = F, center = T)

pc1 <- autoplot(pca_res, data = genomat5, colour = 'Direction', shape = "Cycle", alpha = 0.5) + theme_bw(base_family = "Arial", base_size = 12) + scale_colour_brewer(palette = "Dark2") + theme(panel.grid.minor = element_blank()) + scale_y_reverse() 
pc1

ggsave(filename = "~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/4_plots/PCA_GSDIG.pdf", plot = pc1, dpi = 300, width = 6, height = 4, device = cairo_pdf)



pc2 <- plot(pca_res, type="l", title(main = NULL))
pc2

# pheno

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/3_GWASpoly/")
pheno <- read.csv("Stem_strength.csv")
pheno1 <- pheno[,1, drop = F]

genomat6 <- inner_join(pheno1, genomat5, by = "gen")
dim(genomat6)
genomat6[1:5,1:5]

marks.7 <- genomat6[,-c(1:4)]
rownames(marks.7) <- genomat6$gen
dim(marks.7)

pca_res1 <- prcomp(marks.7, scale. = F, center = T)

pc2 <- autoplot(pca_res1, data = genomat6, colour = 'Direction', shape = "Cycle", alpha = 0.5) + theme_bw(base_family = "Arial", base_size = 12) + scale_colour_brewer(palette = "Dark2") + theme(panel.grid.minor = element_blank()) + scale_x_reverse() 

pc2

ggsave(filename = "~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/4_plots/PCA_500.pdf", plot = pc2, dpi = 300, width = 6, height = 4, device = cairo_pdf)


# PCA6 <- as.data.frame(pca_res$var$coord)
PCA6 <- as.data.frame(pca_res$x[,1:5])
colnames(PCA6) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
PCA6 <- PCA6 %>% rownames_to_column(var = "gen")

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/3_GWASpoly/")
pheno <- read.csv("Stem_strength.csv")
pheno1 <- pheno[,1, drop = F]
pheno2 <- pheno[,c(1:5)]
pheno2 <- inner_join(pheno2, PCA6, by = "gen")

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/3_GWASpoly/")
write.csv(pheno2, "Stem_strength_1.csv", quote = F, row.names = F)


var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
var_explained[1:5]
pc2 <- plot(pca_res, type="l", title(main = NULL))
pc2
# pca ---------------------------------------------------------------------

eig.val <- get_eigenvalue(pca_res)
pca.res <- get_pca_ind(pca_res)

fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 50)) # ok
fviz_pca_var(pca_res, col.var = "black")
fviz_cos2(pca_res, choice = "var", axes = 1:2)


marks.6[1:5,1:5]
geno.scale <- scale(marks.6, center = T, scale = T)
svdgeno <- svd(geno.scale)
PCA <- geno.scale %*% svdgeno$v
PCA[1:5,1:5]

PCA4 <- PCA[,c(1:5)]
head(PCA4)
class(PCA4)
colnames(PCA4) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
PCA4 <- as.data.frame(PCA4)
PCA4 <- PCA4 %>% rownames_to_column(var = "gen")
head(PCA4)

PCA4 # CF ST
PCA7 <- PCA4 # CT SF
PCA5 <- PCA4 # CT ST
PCA6 <- PCA4 # CF SF

cor(PCA6$PC1, PCA4$PC1)
cor(PCA6$PC1, PCA5$PC1)

cor(PCA6$PC1, pheno$PC1)

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/2_pheno/2_BLUPs/")
write.csv(PCA4, "pca_1.csv", quote = F, row.names = F)
