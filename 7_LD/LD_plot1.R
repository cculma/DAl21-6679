# LD by population

library(plyr)
library(ldsep)
library(updog)
library(tidyverse)

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/7_LD/")

load("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/1_updog/mout_2.Rdata")

gp <- format_multidog(x = mout_2, varname = paste0("Pr_", 0:ploidy))
class(gp)
dim(gp)


lev4 <- unique(mout_2$inddf$ind) 
lev5 <- subset(lev4, grepl("^UMN3097", lev4))
lev6 <- subset(lev4, grepl("^UMN3355", lev4))
lev7 <- subset(lev4, grepl("^UMN3358", lev4))
lev8 <- subset(lev4, grepl("^UMN4016", lev4))
lev9 <- subset(lev4, grepl("^UMN4351", lev4))

length(lev7)
lev12 <- list(lev5, lev6, lev7, lev8, lev9)
names(lev12) <- c("UMN3097","UMN3355","UMN3358","UMN4016","UMN4351")
as.data.frame(lev5)

for (i in 1:length(lev12)) {
  
}

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/10_subset_vcf/")
for (i in 1:length(lev12)) {
  data1 <- as.data.frame(lev12[[i]])
  write.table(data1, paste0(names(lev12[i]), ".txt"), col.names = F, row.names = F, quote = F, sep = "\t")
}



gp <- format_multidog(x = mout_2, varname = paste0("Pr_", 0:ploidy))
class(gp)
dim(gp)

gp
arr2 <- alply(arr1, 2, .dims = T)
arr3 <- array( unlist(arr2) , dim = c( 4 ,3 , 2 ))
class(arr3)
dim(arr3)

arr2 <- alply(gp, 3,  .dims = T)

# lev5
lev11 <- list()
for (i in 1:length(arr2)) {
  lev10 <- as.data.frame(arr2[[i]])
  lev10 <- lev10 %>% dplyr::select(all_of(lev5))
  lev11[[length(lev11)+1]] <- lev10
}
names(lev11) <- names(arr2)
arr3 <- array( unlist(lev11) , dim = c(2953 ,length(lev5) ,5))
ldout_UMN3097 <- ldfast(gp = arr3, type = "r2")

# lev6
lev11 <- list()
for (i in 1:length(arr2)) {
  lev10 <- as.data.frame(arr2[[i]])
  lev10 <- lev10 %>% dplyr::select(all_of(lev6))
  lev11[[length(lev11)+1]] <- lev10
}
names(lev11) <- names(arr2)
arr3 <- array( unlist(lev11) , dim = c(2953 ,length(lev6) ,5))
ldout_UMN3355 <- ldfast(gp = arr3, type = "r2")

# lev7
lev11 <- list()
for (i in 1:length(arr2)) {
  lev10 <- as.data.frame(arr2[[i]])
  lev10 <- lev10 %>% dplyr::select(all_of(lev7))
  lev10 <- na.omit(lev10)
  lev11[[length(lev11)+1]] <- lev10
}
names(lev11) <- names(arr2)
arr3 <- array( unlist(lev11) , dim = c(2953 ,length(lev7) ,5))
ldout_UMN3358 <- ldfast(gp = arr3, type = "r2", se = F)
dim(arr3)

# lev8
lev11 <- list()
for (i in 1:length(arr2)) {
  lev10 <- as.data.frame(arr2[[i]])
  lev10 <- lev10 %>% dplyr::select(all_of(lev8))
  lev11[[length(lev11)+1]] <- lev10
}
names(lev11) <- names(arr2)
arr3 <- array( unlist(lev11) , dim = c(2953 ,length(lev8) ,5))
ldout_UMN4016 <- ldfast(gp = arr3, type = "r2")

# lev9
lev11 <- list()
for (i in 1:length(arr2)) {
  lev10 <- as.data.frame(arr2[[i]])
  lev10 <- lev10 %>% dplyr::select(all_of(lev9))
  lev11[[length(lev11)+1]] <- lev10
}
names(lev11) <- names(arr2)
arr3 <- array( unlist(lev11) , dim = c(2953 ,length(lev9) ,5))
ldout_UMN4351 <- ldfast(gp = arr3, type = "r2")

####
# plot
arr4 <- list(ldout_UMN3097, ldout_UMN3355, ldout_UMN3358, ldout_UMN4016, ldout_UMN4351)

names(arr4) <- c("UMN3097","UMN3355","UMN3358","UMN4016","UMN4351")

lev1 <- mout_2$snpdf$snp
length(lev1)
lev1[1:5]

lev2 <- as.data.frame(lev1)
lev2 <- lev2 %>% separate(1, c("Chr", "position"), sep = "_", remove = F, convert = FALSE, extra = "merge")
lev2$Chr <- as.factor(lev2$Chr)
levels(lev2$Chr)

arr5 <- list()
for (i in 1:length(arr4)){
  
  R2 <- arr4[[i]]$ldmat
  dim(R2)
  colnames(R2) <- rownames(R2) <- lev2$lev1
  
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
  colnames(R3)[5] <- "rsq"
  R3 <- R3[R3$dist != "NaN",]
  arr5[[length(arr5)+1]] <- R3
}
names(arr5) <- names(arr4)

length(lev12[[1]])

arr6 <- list()
arr7 <- list() 
rm(N)
for (i in 1:length(arr5)) {
  R3 <- arr5[[3]]
  N = length(lev12[[3]]) * 4
  Cstart <- c(C=0.1)
  modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                  ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), data=R3, start=Cstart, control=nls.control(maxiter=100))
  rho <- summary(modelC)$parameters[1]
  newrsq <- ( (10+rho*R3$dist) / ( (2+rho*R3$dist) * (11+rho*R3$dist) ) ) *
    ( 1 + ( (3+rho * R3$dist) * (12+12*rho*R3$dist + (rho*R3$dist)^2) ) / 
        (2*N*(2+rho*R3$dist) * (11+rho*R3$dist) ) )
  
  newfile <- data.frame(R3$dist, newrsq)

  maxld <- max(newfile$newrsq,na.rm=TRUE) #using max LD value from adj data
  halfdecay = maxld*0.5
  halfdecaydist <- newfile$R3.dist[which.min(abs(newfile$newrsq-halfdecay))]
  newfile <- newfile[order(newfile$R3.dist),]
  arr6[[length(arr6)+1]] <- newfile
  arr7[[length(arr7)+1]] <- halfdecaydist

}

newfile1 <- newfile
