# install.packages(c("updog","ldsep"))
library(updog)
library(ldsep)

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/1_updog/")
load("mout_2.Rdata")

genomat1 <- format_multidog(mout_2, varname = "geno")
dim(genomat1) # 2953 1502
genomat1[1:5,1:7]
class(genomat1)
c1 <- t(genomat1)

geno.ps <- c1 / 4
geno.ps[1:5,1:5]
Mpr <- snp.pruning(M = geno.ps, pruning.thr = 0.9, window.n = 100, by.chrom = F, overlap.n = 10, seed = 1208, iterations = 15)

genomat5 <- Mpr$Mpruned * 4
genomat5[1:10,1:10]
dim(genomat5) # 1502 2946

NZV <- nearZeroVar(genomat5)
genomat5 <- genomat5[,-NZV]
dim(genomat5) # 1502 2425

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/1_updog/")
write.csv(genomat5, "DAl21_pruned.txt", row.names = T, quote = F)
A1 <- read.csv("DAl21_pruned.txt")
A1[1:5,1:5]

genomat3 <- t(genomat7)

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
dim(G6) # 2425 1505
G6[1:5,1:7]

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/3_GWASpoly/")
write.csv(G6, "DAl21_GWAS_2.txt", row.names = F, quote = F)

