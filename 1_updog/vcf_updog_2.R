library(VariantAnnotation)
library(updog)
library(vcfR)

## Read in data
setwd("/home/samac/medin297/msi/4_GSDIG")

## Read in data

vcf <- read.vcfR("DAl21-6679_Allele_match_counts_collapsed_plates1to16_sorted_geno.vcf")

ADmat <- extract.gt(vcf, "AD")
refmat <- strsplit(ADmat, ",")
refmat <- sapply(refmat, "[[", 1)
refmat <- matrix(refmat, nrow=nrow(ADmat))
refmat <- apply(refmat, 2, as.numeric)

altmat <- strsplit(ADmat, ",")
altmat <- sapply(altmat, "[[", 2)
altmat <- matrix(altmat, nrow=nrow(ADmat))
altmat <- apply(altmat, 2, as.numeric)

colnames(refmat) <- colnames(ADmat)
colnames(altmat) <- colnames(ADmat)
rownames(refmat) <- rownames(ADmat)
rownames(altmat) <- rownames(ADmat)
sizemat <- refmat + altmat

sizemat[1:5,1:5]
refmat[1:5,1:5]

mout_2 <- multidog(refmat = refmat,
                 sizemat = sizemat,
                 ploidy = 4,
                 model = "norm",
                 nc = 30)

save(mout_2, file = "mout_2.Rdata")
