rm(list = ls())
library(VariantAnnotation)
library(updog)
library(tidyverse)
library(ldsep)
library(ggplot2)
library(ASRgenomics)
library(caret)
library(ggfortify)
library(data.table)
library(ggpubr)

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/1_updog/")


c1 <- read.table("DAl21_poly_ori.txt", sep = "\t", header = T, check.names = F, row.names = 1)

c1 <- read.table("DAl21_poly.txt", sep = "\t", header = T, check.names = F, row.names = 1)
c1[1:5,1:5]
dim(c1) # 1502 2477

c1 <- as.matrix(c1)
geno.ps <- c1 / 4
geno.ps[1:5,1:5]

Mpr <- snp.pruning(M = geno.ps, pruning.thr = 0.90, window.n = 100, by.chrom = F, overlap.n = 10, seed = 1208, iterations = 20)

genomat7 <- Mpr$Mpruned * 4
genomat7[1:5,1:5]

NZV <- nearZeroVar(genomat7)
genomat7 <- genomat7[,-NZV]
dim(genomat7) # 1502 2455
setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/1_updog/")
write.csv(genomat7, "DAl21_pruned.txt", row.names = T, quote = F)

genomat5 <- Mpr$Mpruned
genomat5[1:5,1:5]

dim(genomat5) # 1502 2476
genomat6 <- genomat5[,c(1:5)]
dim(genomat6)
genomat6[1:5,1:5]
rowSums(genomat6)
colSums(genomat6)

MAF <- as.data.frame((colSums(genomat5)) / 1502) %>% rownames_to_column("marker")
colnames(MAF)[2] <- "MAF"


MAF1 <- MAF %>% dplyr::filter(MAF > 0.5)
MAF2 <- MAF %>% dplyr::filter(MAF <= 0.5)
MAF1$MAF <- 1 - MAF1$MAF

MAF3 <- rbind(MAF1, MAF2)
min(MAF3$MAF) # 0.05009987
max(MAF3$MAF) # 0.5

MAF_01 <- MAF3 %>% dplyr::filter(MAF >= 0.05 & MAF < 0.1)
MAF_02 <- MAF3 %>% dplyr::filter(MAF >= 0.1 & MAF < 0.2)
MAF_03 <- MAF3 %>% dplyr::filter(MAF >= 0.2 & MAF < 0.3)
MAF_04 <- MAF3 %>% dplyr::filter(MAF >= 0.3 & MAF < 0.4)
MAF_05 <- MAF3 %>% dplyr::filter(MAF >= 0.4)
nrow(MAF_01) + nrow(MAF_02) + nrow(MAF_03) + nrow(MAF_04) + nrow(MAF_05)
hist(MAF3$MAF)

c7 <- gghistogram(MAF3, "MAF", fill = "#00AFBB", rug = F, add_density = F, bins = 50) + theme_classic(base_family = "Arial", base_size = 12)

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/4_plots/")
ggsave(filename = "MAF1.pdf", plot = c7, dpi = 300, width = 4, height = 3, device = cairo_pdf)


ggplot(MAF3, aes(x=MAF)) + geom_histogram(color="black", fill="white")

ggplot(MAF3, aes(x=MAF)) + geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(alpha=.2, fill="#FF6666") 

lev4 <- rownames(genomat5)

lev5 <- subset(lev4, grepl("^UMN3097", lev4))
lev6 <- subset(lev4, grepl("^UMN3355", lev4))
lev7 <- subset(lev4, grepl("^UMN3358", lev4))
lev8 <- subset(lev4, grepl("^UMN4016", lev4))
lev9 <- subset(lev4, grepl("^UMN4351", lev4))

length(lev5) # 298
length(lev6) # 363
length(lev7) # 141
length(lev8) # 358
length(lev9) # 342
genomat6 <- genomat5
genomat6 <- genomat6 %>% dplyr::filter(rownames(genomat6) %in% lev3[2])

genomat6 <- genomat5[row.names(genomat5) %in% lev3[2], ]
genomat6 <- subset(genomat5, rownames(genomat5) %in% lev3[2])
genomat6 <- genomat5[c("UMN3097_248","UMN3097_257"), ]
genomat6 <- genomat5[lev3[[1]], ]
dim(genomat6)

lev3 <- list(lev4,lev5,lev6,lev7,lev8,lev9)
lev2 <- list()
for (i in 1:length(lev3)) {
  
  genomat6 <- genomat5[lev3[[i]], ]

  MAF <- as.data.frame((colSums(genomat6)) / length(lev3[[i]])) %>% rownames_to_column("marker")
  colnames(MAF)[2] <- "MAF"
  MAF1 <- MAF %>% dplyr::filter(MAF > 0.5)
  MAF2 <- MAF %>% dplyr::filter(MAF <= 0.5)
  MAF1$MAF <- 1 - MAF1$MAF
  MAF3 <- rbind(MAF1, MAF2)
  lev2[[length(lev2)+1]] <- MAF3
}
names(lev2) <- c("all","UMN3097","UMN3355","UMN3358","UMN4016","UMN4351")
Y2 <- rbindlist(lev2, use.names=TRUE, fill=TRUE, idcol="pop")

ggplot(Y2, aes(x=MAF, color=pop)) + geom_histogram(fill="white")


gghistogram(Y2, x = "MAF", rug = F, fill = "pop") + facet_grid(pop ~ .)


MAF_01 <- Y2 %>% dplyr::filter(MAF >= 0.05 & MAF < 0.1)
MAF_02 <- Y2 %>% dplyr::filter(MAF >= 0.1 & MAF < 0.2)
MAF_03 <- Y2 %>% dplyr::filter(MAF >= 0.2 & MAF < 0.3)
MAF_04 <- Y2 %>% dplyr::filter(MAF >= 0.3 & MAF < 0.4)
MAF_05 <- Y2 %>% dplyr::filter(MAF >= 0.4)

c1 <- count(MAF_01, pop)
c2 <- count(MAF_02, pop)
c3 <- count(MAF_03, pop)
c4 <- count(MAF_04, pop)
c5 <- count(MAF_05, pop)

c1$MAF <- "0.05-0.1"
c2$MAF <- "0.1-0.2"
c3$MAF <- "0.2-0.3"
c4$MAF <- "0.3-0.4"
c5$MAF <- "0.4-0.5"
c6 <- rbind(c1,c2,c3,c4,c5)


c7 <- ggbarplot(c6, "MAF", "n",
          fill = "pop", color = "pop", palette = "Paired",
          label = F,
          position = position_dodge(0.9)) + theme_classic(base_family = "Arial", base_size = 12)

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/4_plots/")
ggsave(filename = "MAF.pdf", plot = c7, dpi = 300, width = 5, height = 3, device = cairo_pdf)
