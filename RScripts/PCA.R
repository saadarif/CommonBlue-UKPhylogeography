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
z <- vcfR2genlight(VCF)

#get pops
popnames <- z@ind.names 
pops <- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[1])}))
pop(z) = pops

w <- tab(z, freq = TRUE, NA.method = "mean")

#perform PCA
w.pca <- dudi.pca(w, scannf = F, scale=F, nf=3)

#PCA plot
col <- funky(15)
s.class(w.pca$li, pop(z) ,xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=2, clabel=.75, grid=FALSE)
#add scree to plot
add.scatter.eig(w.pca$eig[1:10], xax=1, yax=2)

toto <- summary(z)
########################################################
#add latlong information to genlight object

latlong = as.data.frame( z@pop)
latlong$lat <- NULL
latlong$lon <- NULL

for (i in 1:dim(latlong)[1]){
  latlong$lat[i] <- coords[rownames(coords)==latlong[i,1],1]
  latlong$lon[i] <- coords[rownames(coords)==latlong[i,1],2]
  
}

colnames(latlong) <- c("pop","lat", "lon")
#add lat long to gen light object
z@other$latlong <- latlong[,2:3]
#should now be able to use the dartR IBD function
#add sex information as well
popnames <- z@ind.names 
sex<- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[2])}))
ind <- cbind( z@ind.names, pops, sex, latlong[,2:3])
colnames(ind)[1] <- "ind"
z@other$ind.metrics <- ind
#########################################################
#ibd with dartR
library(dartR)
ibdx <- gl.ibd(z, permutations=9999)

#reshape the dist matrices 
library(reshape2)
df <- melt(as.matrix(ibdx$Dgen), varnames = c("pop1", "pop2"))
colnames(df)[3] <- "dgen"
df2 <- melt(as.matrix(ibdx$Dgeo), varnames = c("pop1", "pop2"))
colnames(df2)[3] <- "dgeo"
df$dgeo <- df2$dgeo
df[which(df$dgen>0.1 & df$dgeo>13.0),]
#remove any 0 values
ibddf <- df[-(df$row==df$col),]
#plot
plot(1, type="n", xlab="Log Dist", ylab=expression("F"[st]*"/"*"(1-F"[st]*")"), xlim=c(10.5, 15), ylim=c(0, 0.16))
#plot bernerary and tul in same colours
points(ibddf$dgeo[ibddf$col=="BER" ],ibddf$dgen[ibddf$col=="BER" ], pch=19,col="red" )
points(ibddf$dgeo[ibddf$col=="TUL" ],ibddf$dgen[ ibddf$col=="TUL" ], pch=19,col="blue" )
points(ibddf$dgeo[ibddf$col=="MLG" ],ibddf$dgen[ ibddf$col=="MLG" ], pch=19,col="Khaki" )
points(ibddf$dgeo[ibddf$col=="OBN" ],ibddf$dgen[ ibddf$col=="OBN" ], pch=19,col="salmon" )
points(ibddf$dgeo[ibddf$col=="DGC" ],ibddf$dgen[ ibddf$col=="DGC" ], pch=19,col="black" )
#everything else
points(ibddf$dgeo[ibddf$row!="TUL" & ibddf$col!="TUL" & ibddf$row!="BER" & ibddf$col!="BER"],ibddf$dgen[ibddf$row!="TUL" & ibddf$col!="TUL"& ibddf$row!="BER" & ibddf$col!="BER"], pch=19,col="grey80" )
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
