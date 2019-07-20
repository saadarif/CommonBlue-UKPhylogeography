#PCA for genepop files from stacks
#SA 28th March 2018

library("adegenet")

#Descend into appropriate directory
setwd('/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4/')

# Load the genotypes.
x = read.genepop('./populations.r50.p15_moh_0.65/populations.snps.gen')

# Load the popmap and update the genepop object.
popmap = read.table('../../INFO//popmap.tsv', col.names=c('sample', 'pop'), row.names='sample')
pop(x) = popmap[indNames(x), 'pop']

#Fill in missing values for snp PCA
y <- tab(x, freq = TRUE, NA.method = "mean")

#perform PCA
y.pca <- dudi.pca(y, scannf = F, scale=F, nf=3)

# Do the PCA alternate
#x.scaled = scaleGen(x, NA.method='mean')
#x.pca = dudi.pca(x.scaled, scannf=F, nf=10)


#scree plot
barplot(y.pca$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

#PCA plot
col <- funky(15)
s.class(y.pca$li, pop(x),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=2, clabel=.75, grid=FALSE)
#s.class(y.pca$li, pop(x), xax=1,yax=2, col=rainbow(nPop(x)),csub=2)
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
BWD=c(-3.898247, 50.579755) #devon 7
FRN=c(1.982669, 44.157335) #Will's site in France 6
PCP=c(-4.308992, 51.674219) #Pembrey county park 13
ETB=c( 0.243747, 50.769187) #eastbourne 14
MMS=c(1.255096, 52.168021) #martin's meadows in suffolk 14
RNL=c(-0.281520, 53.254429) #should be Chamber's Farm wood CFW , lincolnshire 10
RHD=c(-1.477433, 54.713907) #Raisby hill durham 13
BMD=c(-1.125439, 51.782563) #Oxford 16


coords <- rbind(MLG, MDC, BER, TUL, DGC, RVS, OBN, BWD, FRN, PCP, ETB, MMS, RNL, RHD, BMD)
colnames(coords) <- c("lat", "lon")

#calculating distance from latlon
Dgeo <- dist(as.matrix(cbind(coords[,1], coords[,2])))

#######################################################
#fwith vcf
VCF <- read.vcfR("./populations.r50.p15_moh_0.65/populations.snps.filter2.0.5.recode.vcf")
z <- vcfR2genind(VCF)

#get pops
popnames <- rownames(z@tab) 
pops <- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[1])}))
pop(z) = pops



w <- tab(z, freq = TRUE, NA.method = "mean")

#perform PCA
w.pca <- dudi.pca(w, scannf = F, scale=F, nf=3)


#PCA plot
col <- funky(15)
s.class(w.pca$li, pop(z) ,xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=2, clabel=.75, grid=FALSE)

toto <- summary(z)
########################################################
#writing out lat long data for use with eems


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
