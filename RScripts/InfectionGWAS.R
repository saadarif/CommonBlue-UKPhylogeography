# marker assocation between ifnected and uninfected individuals in the north
# April 2020 SA

library("adegenet")
library(vcfR)
library(RColorBrewer)
library(poppr)
library(pegas)
library("factoextra")
library(reshape2)
library("ggplot2")
library(tidyverse)
library(pscl)

#perform the PCA first
#use the to 50% missing threshold to include more individuals
VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/populations.snps.filter2.0.5.recode.vcf")
z <- vcfR2genlight(VCF)
popnames <- z@ind.names
pops <- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[1])}))
pops<- str_replace_all(pops, "RNL", "CFW")
pop(z) = pops
ploidy(z) <- 2
sex<- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[2])}))
ind <- cbind( z@ind.names, pops, sex)
colnames(ind)[1] <- "ind"
z@other$ind.metrics <- ind

w <- tab(z, freq = TRUE, NA.method = "mean")

#perform PCA
w.pca <- dudi.pca(w, scannf = F, scale=F, nf=4)

#make a datafame
dudi.pca <- tibble(pc1=w.pca$li[,1], pc2=w.pca$li[,2], pc3=w.pca$li[,3], pc4=w.pca$li[,4],Sample=z@ind.names, sex=sex, pop=pops)

#read in infection data and joins with PCA and marker data
wolbachia <- read.csv("/media/data_disk/PROJECTS/Saad/CommonBlue/Wolbachia/Wolbachia_infection.csv", sep="\t", stringsAsFactors = F)
wol <- wolbachia %>% mutate(infected=ifelse(log(percentage_classified_wolcbachia)>1,1,0))
names(wol)[1] <- "Sample"

#join data sets
infectionP <- inner_join(wol, dudi.pca, by="Sample")

#keep only pops in the north
#infectionP <- infectionP %>% filter(POP %in% c("BER", "TUL", "DGC", "MLG", "OBN", "RVS", "RHD"))

#get snp MAT
mat <- as.matrix(z)
#subset the indivduals
mat <- mat[rownames(mat) %in% infectionP$Sample,]
#reorder snps to match phenotypes
mat <- mat[order(rownames(mat)),]

#logistic regression the long way with PCAs as base
marker_results=data.frame(locus=character(length=dim(mat)[2]), aic=double(length=dim(mat)[2]),
                           lrt_pval=double(length=dim(mat)[2]),
                           n_snps=integer(length=dim(mat)[2]), pseudoR2= double(length=dim(mat)[2]),stringsAsFactors = F)

for (i in 1:dim(mat)[2]){
  #remove missing vals
  miss = is.na(mat[,i])
  baseModel<- glm(infectionP$infected[!miss] ~ infectionP$pc1[!miss] + infectionP$pc2[!miss] , family=binomial(link=logit))
  logmod <- glm(infectionP$infected[!miss] ~ infectionP$pc1[!miss] + infectionP$pc2[!miss]+ mat[!miss,i], family=binomial(link=logit))
  marker_results$locus[i] <- colnames(mat)[i]
  marker_results$aic[i] <-  logmod$aic
  x <- anova(baseModel, logmod, test="Chisq")
  marker_results$lrt_pval[i] <- x$`Pr(>Chi)`[2]
  
  
  #number of snps at locus
  mark = str_split(colnames(mat)[i], "_")[[1]][1]
  marker_results$n_snps[i] <- length(which(VCF@fix[,1] == mark))
  
  #mcfadden PseudoR2
  marker_results$pseudoR2[i] <- pR2(logmod)[4]
}

#plot pvalues
plot(-1*log10(marker_results$lrt_pval))
marker_results[marker_results$lrt_pval<0.05/dim(mat[2]),]

xtabs(~infectionP$infected+mat[,"23250_55"])


#association between the two infection
infectionP$strain <- ifelse(infectionP$POP %in% c("BER", "TUL"), 1, 0)
#also
infectionP$strain[infectionP$Sample %in% c("MLG_f010", "MLG_m002", "OBN_m_110")] <- 1


marker_results2=data.frame(locus=character(length=dim(mat)[2]), aic=double(length=dim(mat)[2]),
                          lrt_pval=double(length=dim(mat)[2]),
                          n_snps=integer(length=dim(mat)[2]), pseudoR2= double(length=dim(mat)[2]),stringsAsFactors = F)

for (i in 1:dim(mat)[2]){
  #remove missing vals
  miss = is.na(mat[,i])
  baseModel<- glm(infectionP$strain[!miss] ~ infectionP$pc1[!miss] + infectionP$pc2[!miss] , family=binomial(link=logit))
  logmod <- glm(infectionP$strain[!miss] ~ infectionP$pc1[!miss] + infectionP$pc2[!miss]+ mat[!miss,i], family=binomial(link=logit))
  marker_results2$locus[i] <- colnames(mat)[i]
  marker_results2$aic[i] <-  logmod$aic
  x <- anova(baseModel, logmod, test="Chisq")
  marker_results2$lrt_pval[i] <- x$`Pr(>Chi)`[2]
  
  
  #number of snps at locus
  mark = str_split(colnames(mat)[i], "_")[[1]][1]
  marker_results2$n_snps[i] <- length(which(VCF@fix[,1] == mark))
  
  #mcfadden PseudoR2
  marker_results2$pseudoR2[i] <- pR2(logmod)[4]
}

#plot pvalues
plot(-1*log10(marker_results2$lrt_pval))
marker_results[marker_results2$lrt_pval<0.05/dim(mat[2]),]
xtabs(~infectionP$strain+infectionP$infected+mat[,"14392_78"])
