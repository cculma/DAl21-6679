# Calculations for expected values of r2 under drift equilibrium [Hill and Weir (1988)
# implemented in Remington, et al. (2001)]

# Calculate and save LD statistics from TASSEL using sliding window of 50 SNPs and maf = 0.05
# read in table with NaN and distance calculated from TASSEL or wherever


###############################################################################
### 		For more detail, please see following reference
###############################################################################
## Remington, D. L., Thornsberry, J. M., Matsuoka, Y.,
## Wilson, L. M., Whitt, S. R., Doebley, J., ... & Buckler, E. S. (2001). Structure of linkage disequilibrium
## and phenotypic associations in the maize genome. 
## Proceedings of the national academy of sciences, 98(20), 11479-11484. 
## https://doi.org/10.1073/pnas.201394398
###############################################################################

library(ldsep)
library(updog)
library(tidyverse)

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/7_LD/")
load("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/1_updog/mout_2.Rdata")

nrow(mout_2$snpdf) # 2953
mout_2$snpdf[1:5,1:5]
class(mout_2) 
summary(mout_2[["snpdf"]][["od"]]) #0.0181405
summary(mout_2[["snpdf"]][["bias"]]) # 1.14082
msub <- filter_snp(mout_2, bias < 1.14082 & od < 0.0181405 & prop_mis < 0.5) 


msub$snpdf[1:5,1:5]
dim(msub$snpdf) # 1761   20

gp <- format_multidog(x = mout_2, varname = paste0("Pr_", 0:ploidy))
class(gp)
dim(gp)

ldout <- ldfast(gp = gp, type = "r2")

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/7_LD/")
# save(ldout, file = "ldout.RData")

# load("ldout.Rdata")

dim(ldout$ldmat)
ldout$ldmat[1:5,1:5]



lev1 <- mout_2$snpdf$snp
length(lev1)
lev1[1:5]

lev2 <- as.data.frame(lev1)
lev2 <- lev2 %>% separate(1, c("Chr", "position"), sep = "_", remove = F, convert = FALSE, extra = "merge")
lev2$Chr <- as.factor(lev2$Chr)
levels(lev2$Chr)


# lev2$marker <- gsub(":", "_", lev2$marker)
# lev2 <- lev2 %>% separate(2, c("Chr", "position"), sep = "_", remove = F, convert = FALSE, extra = "merge")

C1 <- lev2 %>% dplyr::filter(Chr == "chr1.1")
C2 <- lev2 %>% dplyr::filter(Chr == "chr2.1")
C3 <- lev2 %>% dplyr::filter(Chr == "chr3.1")
C4 <- lev2 %>% dplyr::filter(Chr == "chr4.1")
C5 <- lev2 %>% dplyr::filter(Chr == "chr5.1")
C6 <- lev2 %>% dplyr::filter(Chr == "chr6.1")
C7 <- lev2 %>% dplyr::filter(Chr == "chr7.1")
C8 <- lev2 %>% dplyr::filter(Chr == "chr8.1")


R2 <- ldout$ldmat
dim(R2)
length(lev2$lev1)
colnames(R2) <- rownames(R2) <- lev2$lev1

R2[1:5,1:5]
R2.1 <- R2[C1$lev1, C1$lev1]
R2.2 <- R2[C2$lev1, C2$lev1]
R2.3 <- R2[C3$lev1, C3$lev1]
R2.4 <- R2[C4$lev1, C4$lev1]
R2.5 <- R2[C5$lev1, C5$lev1]
R2.6 <- R2[C6$lev1, C6$lev1]
R2.7 <- R2[C7$lev1, C7$lev1]
R2.8 <- R2[C8$lev1, C8$lev1]

R2.1[upper.tri(R2.1)] <- NA
diag(R2.1) <- NA
R2.1 <- reshape2::melt(R2.1, na.rm = T)

R2.2[upper.tri(R2.2)] <- NA
diag(R2.2) <- NA
R2.2 <- reshape2::melt(R2.2, na.rm = T)

R2.3[upper.tri(R2.3)] <- NA
diag(R2.3) <- NA
R2.3 <- reshape2::melt(R2.3, na.rm = T)

R2.4[upper.tri(R2.4)] <- NA
diag(R2.4) <- NA
R2.4 <- reshape2::melt(R2.4, na.rm = T)

R2.5[upper.tri(R2.5)] <- NA
diag(R2.5) <- NA
R2.5 <- reshape2::melt(R2.5, na.rm = T)

R2.6[upper.tri(R2.6)] <- NA
diag(R2.6) <- NA
R2.6 <- reshape2::melt(R2.6, na.rm = T)

R2.7[upper.tri(R2.7)] <- NA
diag(R2.7) <- NA
R2.7 <- reshape2::melt(R2.7, na.rm = T)

R2.8[upper.tri(R2.8)] <- NA
diag(R2.8) <- NA
R2.8 <- reshape2::melt(R2.8, na.rm = T)

R3 <- rbind(R2.1,R2.2,R2.3,R2.4,R2.5,R2.6,R2.7,R2.8)
head(R3)

R3 <- R3 %>% separate(2, c("Locus2", "Position2"), sep = "_", remove = T, convert = FALSE, extra = "merge") %>% separate(1, c("Locus1", "Position1"), sep = "_", remove = T, convert = FALSE, extra = "merge")

R3$Position1 <- as.numeric(R3$Position1)
R3$Position2 <- as.numeric(R3$Position2)
R3$dist <- as.numeric(R3$Position1 - R3$Position2)
str(R3)
dim(R3)
colnames(R3)[5] <- "rsq"
R3 <- R3[R3$dist != "NaN",]


N = 1502 * 2
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), data=R3, start=Cstart, control=nls.control(maxiter=100))

rho <- summary(modelC)$parameters[1]

newrsq <- ( (10+rho*R3$dist) / ( (2+rho*R3$dist) * (11+rho*R3$dist) ) ) *
  ( 1 + ( (3+rho * R3$dist) * (12+12*rho*R3$dist + (rho*R3$dist)^2) ) / 
      (2*N*(2+rho*R3$dist) * (11+rho*R3$dist) ) )

newfile <- data.frame(R3$dist, newrsq)
head(newfile)

maxld <- max(newfile$newrsq,na.rm=TRUE) #using max LD value from adjusted data
halfdecay = maxld*0.5
halfdecaydist <- newfile$R3.dist[which.min(abs(newfile$newrsq-halfdecay))]
newfile <- newfile[order(newfile$R3.dist),]
newfile1 <- newfile


R3$dist
R3$rsq

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/7_LD/")

pdf("LD_DAI1.pdf", height=5, width = 5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(R3$dist, R3$rsq, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
lines(newfile$R3.dist, newfile$newrsq, col="red", lwd=2)
abline(h=0.1, col="blue") # if you need to add horizontal line
abline(v=halfdecaydist, col="green")
mtext(round(halfdecaydist,2), side=1, line=0.05, at=halfdecaydist, cex=0.75, col="green")
dev.off()

head(newfile)
newfile$Dist_kb <- newfile$R3.dist/1000
head(R3)

ggplot(newfile, aes(x = R3.dist, y = newrsq)) + geom_line() + xlim(0, 1000000) + theme_bw() + geom_hline(yintercept=0.1, linetype="dashed", color = "gray") + geom_vline(xintercept = halfdecaydist)

plot1 <- ggplot(newfile, aes(x = Dist_kb, y = newrsq)) + geom_line() + ylim(0, 0.2) + scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, by = 200)) + theme_classic(base_family = "Arial", base_size = 12) + geom_hline(yintercept=0.1, linetype="dashed", color = "gray") + geom_vline(xintercept = (halfdecaydist/1000), linetype="dashed", color = "blue") + labs(y = expression(LD ~ (r^2)), x = "Distance (kb)")

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/7_LD/")
ggsave(filename = "LD_DAI1.jpg", plot = plot1, width = 4, height = 4)

ggsave(filename = "LD_DAI1.pdf", plot = plot1, width = 4, height = 4, device = cairo_pdf)



lev3 <- mout_2$inddf$ind
