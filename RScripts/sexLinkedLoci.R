#detetcing sex linked loxi
#SA MARCH 2020


#Do the PCA first

library("adegenet")
library(vcfR)
library("ade4")
#library(ape)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library("factoextra")
library(pscl) #pseudoR2

setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/scripts/FigScripts")

#Long/Lats for populations
#Levels: MLG MDC BER TUL DGC RVS OBN
#Long/Late in same order to nearest site
MLG=c(-5.834874, 56.991962)# check field notes 16
MDC=c(-3.843251, 53.295675) #12
BER=c(-7.213534, 57.713854)#16
TUL=c(-7.025334, 58.185530)# check field notes 16
DGC=c(-4.016982, 57.878096)#14
RVS=c(-2.596053, 56.023971)#6
OBN=c(-5.482266, 56.433569)#14
BWD=c(-3.898247, 50.579755) #devon 7
FRN=c(1.982669, 44.157335) #Will's site in France 6
PCP=c(-4.308992, 51.674219) #Pembrey county park 13
ETB=c( 0.243747, 50.769187) #eastbourne 14
MMS=c(1.255096, 52.168021) #martin's meadows in suffolk 14
CFW=c(-0.281520, 53.254429) #should be Chamber's Farm wood CFW , lincolnshire 10
RHD=c(-1.477433, 54.713907) #Raisby hill durham 13
BMD=c(-1.125439, 51.782563) #Oxford 16


coords <- rbind(MLG, MDC, BER, TUL, DGC, RVS, OBN, BWD, FRN, PCP, ETB, MMS, CFW, RHD, BMD)
colnames(coords) <- c("lon", "lat")

#calculating distance from latlon
#Dgeo <- dist(as.matrix(cbind(coords[,1], coords[,2])))
VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/populations.snps.filter2.0.5.recode.vcf")
z <- vcfR2genlight(VCF)

#get pops and label 
popnames <- z@ind.names
pops <- sapply(strsplit(popnames, '_'), function(x){paste0(x[1])})
pops<- str_replace_all(pops, "RNL", "CFW")
#replace error sex in ETB_f_188
z@ind.names <- str_replace_all(z@ind.names, "ETB_f_188", "ETB_m_188")
z@ind.names <- str_replace_all(z@ind.names, "PCP_f_161", "PCP_m_161")
pop(z) = as.factor(pops) 
ploidy(z) <- 2

#add latlong information to genlight object
latlong = as.data.frame( z@pop)
latlong$lat <- NULL
latlong$lon <- NULL

for (i in 1:dim(latlong)[1]){
  latlong$lat[i] <- coords[rownames(coords)==latlong[i,1],1]
  latlong$lon[i] <- coords[rownames(coords)==latlong[i,1],2]
  
}

colnames(latlong) <- c("POP","lat", "lon")
#add lat long to gen light object
z@other$latlong <- latlong[,2:3]
#should now be able to use the dartR IBD function
#add sex information as well
popnames <- z@ind.names 
sex<- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[2])}))

ind <- cbind( z@ind.names, pops, sex, latlong[,2:3])
colnames(ind)[1] <- "ind"
z@other$ind.metrics <- ind
w <- tab(z, freq = TRUE, NA.method = "mean")

#perform PCA
w.pca <- dudi.pca(w, scannf = F, scale=F, nf=4)

dudi.pca <- tibble(pc1=w.pca$li[,1], pc2=w.pca$li[,2], pc3=w.pca$li[,3], pc4=w.pca$li[,4],lat=z@other$latlong[,2], pop=pops, ind=z@other$ind.metrics$ind)


#--------------------------------------------------------------------------------------
#find sex linked loci by association
mat <- as.matrix(z)

marker_results=data.frame(locus=character(length=dim(mat)[2]), aic=double(length=dim(mat)[2]), no_hets=integer(length=dim(mat)[2]), stringsAsFactors = F)

for (i in 1:dim(mat)[2]){
  cross<-xtabs(~mat[,i]+z@other$ind.metrics$sex)
  markers_d=as.factor(rownames(cross))
  logmod <- glm(cross ~ markers_d, family=binomial(link=logit))
  marker_results$locus[i] <- colnames(mat)[i]
  marker_results$aic[i] <-  logmod$aic
}

#for which loci are males always homozygotes
for (i in 1:dim(mat)[2]){
  cross<-xtabs(~mat[,i]+z@other$ind.metrics$sex)
  markers_d=as.factor(rownames(cross))
  no_hets = sum(cross[2:dim(cross)[1],2])
  marker_results$no_hets[i] <-  no_hets
}

hist(marker_results$aic)

marker_results[marker_results$aic<10,]



#see stats for marker
grab=which(VCF@fix[,1] == "9544")
#extract info
gt<- sapply(strsplit(VCF@gt[grab[1],2:length(VCF@gt[grab[1],])], ':'), function(x){paste0(x[1])})
gt[gt=="./."] <- NA
dp <- as.numeric(sapply(strsplit(VCF@gt[grab[1],2:length(VCF@gt[grab[1],])], ':'), function(x){paste0(x[2])}))
ad <- (sapply(strsplit(VCF@gt[grab[1],2:length(VCF@gt[grab[1],])], ':'), function(x){paste0(x[3])}))
ad[ad=="."] <- NA
#get sex and pop
pop <- (sapply(strsplit(names(gt), '_'), function(x){paste0(x[1])}))
#replace error sex in ETB_f_188
names(gt)<- str_replace_all(names(gt), "ETB_f_188", "ETB_m_188")
names(gt)<- str_replace_all(names(gt), "PCP_f_161", "PCP_m_161")
sex <- (sapply(strsplit(names(gt), '_'), function(x){paste0(x[2])}))
geno <- data.frame( pop=pop, sex=sex, gt=gt, readepth=dp, alleledp=ad)
#averageref <- mean(as.numeric(sapply(strsplit(as.character(geno$alleledp[geno$sex=="f"]), ','), function(x){paste0(x[1])})),  na.rm=T)
#averagesnp <- mean(as.numeric(sapply(strsplit(as.character(geno$alleledp[geno$sex=="f"]), ','), function(x){paste0(x[2])})),  na.rm=T)
xtabs(~geno$pop+geno$sex+geno$gt)

#results

#1591   9681_65 8.360417       0  single snp 2 FRNs female
#1609   9781_55 8.216726       0  single snp good
#1819  11011_27 8.140440       0  single snp good
#3727  22073_63 8.286326       0  1 FRN female 1 BMD female
#3728  22073_76 8.348842       0  1 FRN female 1 BMD female
#4168  24861_17 8.053044       0  good identical on all snps
#4169  24861_38 8.053044       0  good
#4170  24861_43 8.053044       0  good 
#4298  25851_91 8.355518       0  2 FRN females
#5548 396673_60 8.595931       0  1 BMD 1 DGC 2 FRN 3 MDC

#5177 222601_60 7.957453       0  multiple SNPs
#5477 356589_17 8.276631       0  multiple SNPS
#2682  16136_81 8.342999       0  multiple snps on locys
#3593  21117_46 9.029789      15  male hets multiple snps on locus
#958    6424_29 9.043133      11  male hets, multiple snps on locus
#1569   9544_25 7.847342       0  multiple snps on locud

#logistic regression the long way
marker_results2=data.frame(locus=character(length=dim(mat)[2]), aic=double(length=dim(mat)[2]), pseudoR2=double(length=dim(mat)[2]),
                           n_snps=integer(length=dim(mat)[2]),no_hets=integer(length=dim(mat)[2]), stringsAsFactors = F)

for (i in 1:dim(mat)[2]){
  
  logmod <- glm(z@other$ind.metrics$sex ~ mat[,i], family=binomial(link=logit))
  marker_results2$locus[i] <- colnames(mat)[i]
  #number of snps at locus
  mark = str_split(colnames(mat)[i], "_")[[1]][1]
  marker_results2$n_snps[i] <- length(which(VCF@fix[,1] == mark))
  #aic
  marker_results2$aic[i] <-  logmod$aic
  #pseudoR2
  marker_results2$pseudoR2[i] <-pR2(logmod)[4]
  
  #for which loci are males always homozygotes
  cross<-xtabs(~mat[,i]+z@other$ind.metrics$sex)
  markers_d=as.factor(rownames(cross))
  no_hets = sum(cross[2:dim(cross)[1],2])
  marker_results2$no_hets[i] <-  no_hets
}


hist(marker_results2$aic)
marker_results2[marker_results2$aic<80,]
marker_results[marker_results$aic<10,]
head(marker_results2[order(marker_results2$aic, decreasing = F),],10)

#logistic regression the long way with PCAs as base
marker_results3=data.frame(locus=character(length=dim(mat)[2]), aic=double(length=dim(mat)[2]), no_f_hets=integer(length=dim(mat)[2]),
                           lrt_pval=double(length=dim(mat)[2]), no_m_hets=integer(length=dim(mat)[2]), no_genotypes=integer(length=dim(mat)[2]),
                           n_snps=integer(length=dim(mat)[2]), pseudoR2= double(length=dim(mat)[2]),stringsAsFactors = F)

for (i in 1:dim(mat)[2]){
  #remove missing vals
  miss = is.na(mat[,i])
  baseModel<- glm(z@other$ind.metrics$sex[!miss] ~ dudi.pca$pc1[!miss] + dudi.pca$pc2[!miss] , family=binomial(link=logit))
  logmod <- glm(z@other$ind.metrics$sex[!miss] ~  dudi.pca$pc1[!miss] + dudi.pca$pc2[!miss]+ mat[!miss,i], family=binomial(link=logit))
  marker_results3$locus[i] <- colnames(mat)[i]
  marker_results3$aic[i] <-  logmod$aic
  x <- anova(baseModel, logmod, test="Chisq")
  marker_results3$lrt_pval[i] <- x$`Pr(>Chi)`[2]
  
  #for which loci are males always homozygotes
  cross<-xtabs(~mat[,i]+z@other$ind.metrics$sex)
  markers_d=as.factor(rownames(cross))
  no_hets = sum(cross[2:dim(cross)[1],2])
  marker_results3$no_m_hets[i] <-  no_hets
  marker_results3$no_f_hets[i] <- sum(cross[2:dim(cross)[1],1])
  #number of genotypes should be two
  marker_results3$no_genotypes[i] <-  dim(cross)[1]
  
  #number of snps at locus
  mark = str_split(colnames(mat)[i], "_")[[1]][1]
  marker_results3$n_snps[i] <- length(which(VCF@fix[,1] == mark))
  
  #mcfadden PseudoR2
  marker_results3$pseudoR2[i] <- pR2(logmod)[4]
}

#plot(-1*log10(marker_results3$lrt_pval))
#marker_results3[-1*log10(marker_results3$lrt_pval)>10,]

#filter

sexmarkers <- marker_results3 %>% filter(lrt_pval < 0.05/5592 & no_genotypes==2 & no_m_hets==0 ) #

#keep only makers with consitent snps

sexmarkers_s <- sexmarkers %>% filter(locus!="9544_25" & locus != "16136_81" & locus != "222601_60" & locus != "356589_17"  )

#generate tidy table of sex-specific marker genotypes

#get unique loci numbers
uloci=unique(sapply(sexmarkers_s$locus, function(x){str_split(x, "_")[[1]][1]}))

#Get infection information
wolbachia <- read.csv("../../Wolbachia/Wolbachia_infection.csv", sep="\t")
wol <- wolbachia %>% mutate(infected=ifelse(log(percentage_classified_wolcbachia)>1,"wIca2","Uninfected"))
wol$infected[wol$POP=="BER"  | wol$POP=="TUL"] <- "wIca1"
wol$infected[wol$X.sample=="MLG_f_010" | wol$X.sample=="MLG_m_002" | wol$X.sample=="OBN_m_110"] <- "wIca1"

#correct erroneous female assignments
wol$X.sample <- str_replace_all(wol$X.sample, "ETB_f_188", "ETB_m_188")
wol$X.sample<- str_replace_all(wol$X.sample, "PCP_f_161", "PCP_m_161")

#data frame to store data
genotypes=data.frame()

for (marker in uloci){
  #get the locus
  grab=which(VCF@fix[,1] == marker)
  for (i in 1:length(grab)){
    
    gt<- sapply(strsplit(VCF@gt[grab[i],2:length(VCF@gt[grab[i],])], ':'), function(x){paste0(x[1])})
    gt[gt=="./."] <- NA # missing genotypes
    #get sex and pop
    pop <- (sapply(strsplit(names(gt), '_'), function(x){paste0(x[1])}))
    #replace error sex in ETB_f_188 and PCP
    names(gt)<- str_replace_all(names(gt), "ETB_f_188", "ETB_m_188")
    names(gt)<- str_replace_all(names(gt), "PCP_f_161", "PCP_m_161")
    sex <- (sapply(strsplit(names(gt), '_'), function(x){paste0(x[2])}))
    geno <- data.frame(names=names(gt), pop=pop, sex=sex, gt=gt)
    #lexicogrpahical order data frame to it matches wiht the infection
    geno <- geno[order(geno$names),]
    geno$inf <- wol$infected[wol$X.sample %in% geno$names]
    genoT <- xtabs(~geno$inf+geno$sex+geno$gt)
    genotypes <- rbind(genotypes, genoT[,,1])
    genotypes <- rbind(genotypes, genoT[,,2])
    }
}

genotypes$geno <- rep(c(rep("Homozygote", 6), rep("Heterozygote", 6)), length(sexmarkers_s$locus))
genotypes$markers <- as.vector(sapply(sexmarkers_s$locus, function(x) {rep(x,12)}))

genoWide<- reshape(genotypes, idvar = c( "geno.sex",  "geno", "markers"), timevar="geno.inf", direction="wide")

#grouping by Outer Hebrides, Northern Mainland Scotland and the rest for plotting
#groupedGeno <- data.frame(Marker=genoWide$markers, M.Sex= genoWide$geno.sex, Geno= genoWide$geno, Outer_Hebrides=genoWide$Freq.BER+ genoWide$Freq.TUL, 
#                          NScot_mainland= genoWide$Freq.DGC+genoWide$Freq.MLG+genoWide$Freq.OBN+genoWide$Freq.RHD,
#                          Uninfected = genoWide$Freq.RVS+genoWide$Freq.BMD+genoWide$Freq.BWD+genoWide$Freq.ETB+genoWide$Freq.FRN+genoWide$Freq.MDC+genoWide$Freq.MMS+genoWide$Freq.PCP+genoWide$Freq.RNL)
#convert to long format for saving
#groupedGenolong <- gather(groupedGeno, Population, Counts, Outer_Hebrides:Uninfected , factor_key = TRUE)


write.csv(genoWide, "TableS7_sex_specific.csv", row.names = F, quote=F)
write.csv(genotypes, "sexmarkerdata.csv", row.names = F, quote=F)
#rename and reorder columns



#without accounting for pop stucture

marker_results4=data.frame(locus=character(length=dim(mat)[2]), aic=double(length=dim(mat)[2]), 
                           lrt_pval=double(length=dim(mat)[2]), no_m_hets=integer(length=dim(mat)[2]), no_f_hets=integer(length=dim(mat)[2]),
                           no_genotypes=integer(length=dim(mat)[2]),
                           n_snps=integer(length=dim(mat)[2]), pseudoR2= double(length=dim(mat)[2]),stringsAsFactors = F)

for (i in 1:dim(mat)[2]){
  #remove missing vals
  miss = is.na(mat[,i])
  baseModel<- glm(z@other$ind.metrics$sex[!miss] ~ 1 , family=binomial(link=logit))
  logmod <- glm(z@other$ind.metrics$sex[!miss] ~ mat[!miss,i], family=binomial(link=logit))
  marker_results4$locus[i] <- colnames(mat)[i]
  marker_results4$aic[i] <-  logmod$aic
  x <- anova(baseModel, logmod, test="Chisq")
  marker_results4$lrt_pval[i] <- x$`Pr(>Chi)`[2]
  
  #for which loci are males always homozygotes
  cross<-xtabs(~mat[,i]+z@other$ind.metrics$sex)
  #markers_d=as.factor(rownames(cross))
  no_hets = sum(cross[2:dim(cross)[1],2])
  marker_results4$no_m_hets[i] <-  no_hets
  
  #number of genotypes should be two
  marker_results4$no_genotypes[i] <-  dim(cross)[1]
  
  #number of snps at locus
  mark = str_split(colnames(mat)[i], "_")[[1]][1]
  marker_results4$n_snps[i] <- length(which(VCF@fix[,1] == mark))
  
  #mcfadden PseudoR2
  marker_results4$pseudoR2[i] <- pR2(logmod)[4]
}

plot(-1*log10(marker_results4$lrt_pval))
#marker_results3[-1*log10(marker_results3$lrt_pval)>10,]

#filter

sexmarkers4 <- marker_results4 %>% filter(lrt_pval < 0.05/5592 & no_genotypes==2 & no_hets==0 )

