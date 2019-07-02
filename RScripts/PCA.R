#PCA for genepop files from stacks
#SA 28th March 2018

library("adegenet")

#Descend into appropriate directory
setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks/stacks.denovoM3n3/populations.r70.p7_MAF05_MOH065_NOTSINGLESNP_noContam//")
#Make sure the .gen suffix is added to genepop formatted files from stacks
x =read.genepop('populations.snps.gen')

#Extract population labels
pop(x) =sapply(strsplit(indNames(x), '_'), function(x){x[1]})

#Fill in missing values for snp PCA
y <- tab(x, freq = TRUE, NA.method = "mean")

#perform PCA
y.pca <- dudi.pca(y, scannf = F, scale=F)

#scree plot
barplot(y.pca$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

#PCA plot
s.class(y.pca$li, pop(x), col=rainbow(nPop(x)))
#add scree to plot
add.scatter.eig(y.pca$eig[1:10], xax=1, yax=2)
toto <- summary(x)
barplot(toto$Hexp - toto$Hobs, main = "Heterozygosity: expected-observed",
         ylab = "Hexp - Hobs")
#Look at IBD and sPCA tutorials
#------------------------------------------------------------------------------------------------------------------
#Long/Lats for populations
pop(x)
#Levels: MLG MDC BER TUL DGC RVS OBN
#Long/Late in same order to nearest site
MLG=c(-5.834874, 56.991962)# check field notes 16
MDC=c(-3.843251, 53.295675) #12
BER=c(-7.213534, 57.713854)#16
TUL=c(-7.025334, 58.185530)# check field notes 16
DGC=c(-4.016982, 57.878096)#14
RVS=c(-2.596053, 56.023971)#6
OBN=c(-5.482266, 56.433569)#14
BMD=c(-1.125439, 51.782563) #Oxford 16
BWD=c(-3.898247, 50.579755) #devon 7
ETB=c( 0.243747, 50.769187) #eastbourne 14
FRN=c(1.982669, 44.157335) #Will's site in France 6
MMS=c(1.255096, 52.168021) #martin's meadows in suffolk 14
PCP=c(-4.308992, 51.674219) #Pembrey county park 13
RHD=c(-1.477433, 54.713907) #Raisby hill durham 13
RNL=c(-0.281520, 53.254429) #should be Chamber's Farm wood CFW , lincolnshire 10



coords <- rbind(MLG, MDC, BER, TUL,DGC,RVS,OBN, BMD, BWD, ETB, FRN, MMS, PCP, RHD, RNL)


longlat = data.frame(pop=character(length = length(pop(x))), long=numeric(length = length(pop(x))), 
                     lat=numeric(length=length(pop(x))), stringsAsFactors = F)

for (i in 1:length(indNames(x))){
  place=unlist(strsplit(indNames(x)[i], '_'))[1]
  index=which(rownames(coords)==place)
  longlat$pop[i]=indNames(x)[i]
  longlat$long[i]=coords[index,1]
  longlat$lat[i]=coords[index,2]
  
}

write.table(longlat, "/media/data_disk/PROJECTS/Saad/CommonBlue/eemsAnalysis/coords", quote=F, row.names = F, col.names = F)
