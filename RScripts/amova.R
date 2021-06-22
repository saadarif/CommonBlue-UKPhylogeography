#AMOVA
#DEC 2020

library(adegenet)
library(poppr) ## Ignore 'OMP parallel support' info...
library(ape)
library(vcfR)
library(dartR)
library(tidyverse)

#read in the snp data 
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
coords <- as.data.frame(coords)
#convert rownames to column
coords <- tibble::rownames_to_column(coords, "POP")

VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/BayeScan/p15r50miss25/noFRNRVS/populations.snps.filter2.0.25.recode_outlier.recode.vcf")
z <- vcfR2genlight(VCF)

#get pops and label 
popnames <- z@ind.names
pops <- sapply(strsplit(popnames, '_'), function(x){paste0(x[1])})
pops<- str_replace_all(pops, "RNL", "CFW")
pop(z) = as.factor(pops) 
ploidy(z) <- 2

#add latlong information to genlight object
latlong = as.data.frame( z@pop)
latlong$lat <- NULL
latlong$lon <- NULL

for (i in 1:dim(latlong)[1]){
  latlong$lat[i] <- coords[coords$POP==latlong[i,1],2]
  latlong$lon[i] <- coords[coords$POP==latlong[i,1],3]
  
}

names(latlong) <- c("POP","lat", "lon")
#add lat long to gen light object
z@other$latlong <- latlong[,2:3]
#add sex information as well
popnames <- z@ind.names 
sex<- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[2])}))
ind <- cbind( z@ind.names, pops, sex, latlong[,2:3])
ind$`z@ind.names` <- as.character(ind$`z@ind.names` )
colnames(ind)[1] <- "ind"
z@other$ind.metrics <- ind

#read in the co1 data
mitotypes <- read.csv("/media/data_disk/PROJECTS/Saad/CommonBlue/scripts/FigScripts/Tables/TableS3_RADstats_short.csv", header=T, na.strings=c("", "NA"), row.names = NULL, stringsAsFactors = F)
#get rid of the NAs
mitotypes <- mitotypes[!is.na(mitotypes$CO1.Haplotype),]
names(mitotypes)[1] <- "ind"
mitotypes$CO1.hap2[mitotypes$CO1.Haplotype=="Alicante-Provence" | mitotypes$CO1.Haplotype=="Iberia-Italy" ] <- "southernRef"
mitotypes$CO1.hap2[mitotypes$CO1.Haplotype=="Palearctic" ] <- "northerbRef"
#add the haplogroup classification to individual attributes
z@other$ind.metrics<- left_join(z@other$ind.metrics, mitotypes)

#--------------------------
#basic summaries and ftats
library(hierfstat)
geni <- gl2gi(z)
pop(geni)=pops
summary(geni)
basic.stats(geni)

#by population
n.pop <- seppop(geni) 

mean.hobs <- do.call("c", lapply(n.pop, function(x) mean(summary(x)$Hobs))) 
mean.hobs[is.nan(mean.hobs)] <- NA 

mean.hs <- do.call("c", lapply(n.pop, function(x) mean(summary(x)$Hexp)))
mean.hs[is.nan(mean.hs)] <- NA
#---------------------------------------------
#AMOVA
#RETAIN INDVIDUALS WITH mtDNA data only
#get available data
zmt <- z[!is.na(z@other$ind.metrics$CO1.Haplotype), ]

#-------------------------------------------------------
#AMOVA
#remove the french population for this analysis
zmt1 <- zmt[zmt$pop!="FRN"]
zmt1$strata <- zmt1$other$ind.metrics[,c(2,10,15)]

table(strata(zmt1, ~CO1.Haplotype, combine = FALSE))

CBamova <- poppr.amova(zmt1,~CO1.Haplotype, within )=
CBamovatest <- randtest(CBamova, nrepet=9999)
Cplot(CBamovatest)
CBamovatest

#---------------------------------------------------------
#grouping ib-it and a-p totheret

table(strata(zmt1, ~CO1.hap2, combine = FALSE))
CBamova <- poppr.amova(zmt1,~CO1.hap2, within=F )

#------- 
#removign potentially discordant populations
zmt2 <- zmt1[zmt1$pop!="CFW"]
zmt2 <- zmt2[zmt2$pop !="MDC"]
zmt2 <- zmt2[zmt2$pop != "MMS"]
