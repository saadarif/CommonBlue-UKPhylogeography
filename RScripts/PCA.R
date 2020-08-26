#PCA for genepop files from stacks
#SA 28th March 2018

library("adegenet")
library(vcfR)
library(dartR)
library("ade4")
library(ape)
library(RColorBrewer)
library(poppr)
library(pegas)
library("factoextra")
library(reshape2)
library("ggplot2")
library(tidyverse)

#Descend into appropriate directory
setwd('/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m5_M4_n4_new//populations.r50.p15_moh_0.65/')

# Load the genotypes.
#x = read.genepop('./populations.r50.p15_moh_0.65/populations.snps.gen')
# Load the popmap and update the genepop object.
#popmap = read.table('../../INFO//popmap.tsv', col.names=c('sample', 'pop'), row.names='sample')
#pop(x) = popmap[indNames(x), 'pop']
#Fill in missing values for snp PCA
#y <- tab(x, freq = TRUE, NA.method = "mean")
#perform PCA
#y.pca <- dudi.pca(y, scannf = F, scale=F, nf=3)
# Do the PCA alternate
#x.scaled = scaleGen(x, NA.method='mean')
#x.pca = dudi.pca(x.scaled, scannf=F, nf=10)
#scree plot
#barplot(y.pca$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))
#PCA plot
#col <- funky(15)
#s.class(y.pca$li, pop(x),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
#        cstar=0, cpoint=2, clabel=.75, grid=FALSE)
#s.class(y.pca$li, pop(x), xax=1,yax=2, col=rainbow(nPop(x)),csub=2)
#add scree to plot
#add.scatter.eig(y.pca$eig[1:10], xax=1, yax=2)

#toto <- summary(x)
#barplot(toto$Hexp - toto$Hobs, main = "Heterozygosity: expected-observed",
#         ylab = "Hexp - Hobs")
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
CFW=c(-0.281520, 53.254429) #should be Chamber's Farm wood CFW , lincolnshire 10
RHD=c(-1.477433, 54.713907) #Raisby hill durham 13
BMD=c(-1.125439, 51.782563) #Oxford 16


coords <- rbind(MLG, MDC, BER, TUL, DGC, RVS, OBN, BWD, FRN, PCP, ETB, MMS, CFW, RHD, BMD)
colnames(coords) <- c("lon", "lat")

#calculating distance from latlon
Dgeo <- dist(as.matrix(cbind(coords[,1], coords[,2])))

#######################################################
#fwith vcf
VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m5_M4_n4_new/populations.r50.p15_moh_0.65/populations.snps.filter2.0.25.recode.vcf")
z <- vcfR2genlight(VCF)

#get pops
popnames <- z@ind.names
#popnames <- rownames(z@tab)
pops <- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[1])}))
pops<- str_replace_all(pops, "RNL", "CFW")
pop(z) = pops
ploidy(z) <- 2
#add latlong information to genlight object
latlong = data.frame( pop=z@pop, lat=numeric(length=length(z@pop)), lon=numeric(length=length(z@pop)))

for (i in 1:dim(latlong)[1]){
  print(i)
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

w <- tab(z, freq = TRUE, NA.method = "mean")

#perform PCA
w.pca <- dudi.pca(w, scannf = F, scale=F, nf=4)

#PCA plot
col <- funky(15)
s.class(w.pca$li, pop(z) ,xax=1,yax=2, col=transp(col,.7), axesell=F,
        cstar=0, cpoint=1.5, clabel=.75, grid=FALSE, addaxes = T) 
title(xlab= "PC1 (7%)", ylab="PC2 (2%)")

#add scree to plot
add.scatter.eig(w.pca$eig[1:10], xax=1, yax=2)
fviz_pca_var(w.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE ,    # Avoid text overlapping
             select.var = list(name = NULL, cos2 = 20, contrib = NULL),
             col.circle = "grey70"
)

#high loadding markers on pc1
quantile(w.pca$c1[,1], probs=c(0.01, 0.99))
markersPCA1 <-w.pca$c1[which(abs(w.pca$c1[,1]) >0.05),]
markersPCA1[order(markersPCA1$CS1),]
which(VCF@fix[,1]=="21281")
VCF@fix[VCF@fix[,1]=="25996",]
# VCF@gt[2066,]
#second axis
quantile(w.pca$c1[,2], probs=c(0.01, 0.99))
markersPCA2 <-w.pca$c1[which(abs(w.pca$c1[,2]) >0.05),]
markersPCA2[order(markersPCA2$CS2),]
which(VCF@fix[,1]=="23511")
VCF@fix[VCF@fix[,1]=="23511",]

#get genotypes in usable formate quick
mat <- as.matrix(z)
mat[,c("25996_6")]
xtabs(~mat[,c("21281_19")]+pop(z))

#w.pca$c1[which(abs(w.pca$c1[,1]) > 0.035
toto <- summary(z)
########################################################
#neighbour joinging tree
gl.tree.nj(z, type="phylogram", outgroup = "FRN")

#UPGMA tree with bootstrap support
tree <- aboot(z, tree = "upgma", distance=bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
cols <- brewer.pal(n = nPop(z), name = "rainbow")
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(z)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
#legend('topleft', legend = c("CA","OR","WA"), fill = cols, border = FALSE, bty = "n", cex = 2)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")
##########################################################
#fixed differences between pops
library(dartR)
#drop pops wiht < 4 inds
glnew3 <- gl.drop.pop(z, pop.list=c("RVS"))
fd <- gl.fixed.diff(glnew3, tloc=0.05, test=TRUE, delta=0.02, reps=100, v=1 )
#########################################################
#ibd with dartR
library(dartR)
ibdx <- gl.ibd(z, permutations=99)

#reshape the dist matrices 
library(reshape2)
df <- melt(as.matrix(ibdx$Dgen), varnames = c("pop1", "pop2"))
colnames(df)[3] <- "dgen"
df2 <- melt(as.matrix(ibdx$Dgeo), varnames = c("pop1", "pop2"))
colnames(df2)[3] <- "dgeo"
df$dgeo <- df2$dgeo
df[which(df$dgen>0.1 & df$dgeo>13.0),]
#remove redudnacies
ibddf <- df[as.numeric(df$pop2)>as.numeric(df$pop1),]
#delete duplicate comparisons


ibddf$comp <- paste0(ibddf$pop1, "-", ibddf$pop2)
#plot
plot(1, type="n", xlab="Log Dist", ylab=expression("F"[st]*"/"*"(1-F"[st]*")"), xlim=c(10.5, 15), ylim=c(0, 0.16))
#plot bernerary and tul in same colours
#plot all points first
points(ibddf$dgeo, ibddf$dgen)
abline(lm(ibddf$dgen~ibddf$dgeo))

#hilight north south compairsons
pattern=c("TUL", "MLG", "DGC", "OBN", "BER")
pattern1=c("CFW", "MDC", "MMS", "PCP", "BMD", "ETB", "BWD", "FRN")

points(ibddf$dgeo[grepl(paste(pattern,collapse="|"), ibddf$comp) &  grepl(paste(pattern1,collapse="|"), ibddf$comp)],ibddf$dgen[grepl(paste(pattern,collapse="|"), ibddf$comp) &  grepl(paste(pattern1,collapse="|"), ibddf$comp)], pch=19,col="black" )


#highlight OH vs mainland comparisons

pattern=c("MLG", "DGC", "OBN")
pattern1= c("BER", "TUL")
points(ibddf$dgeo[grepl(paste(pattern,collapse="|"), ibddf$comp) &  grepl(paste(pattern1,collapse="|"), ibddf$comp)],ibddf$dgen[grepl(paste(pattern,collapse="|"), ibddf$comp) &  grepl(paste(pattern1,collapse="|"), ibddf$comp)], pch=19,col="grey80" )

points(ibddf$dgeo[ grepl("TUL", ibddf$comp) & grepl("BER", ibddf$comp)],ibddf$dgen[grepl("TUL", ibddf$comp) & grepl("BER", ibddf$comp)], pch=19,col="red" )



#what happens to IBD without North-South Comparisons
pattern=c("TUL", "MLG", "DGC", "OBN", "BER")
pattern1=c("CFW", "MDC", "MMS", "PCP", "BMD", "ETB", "BWD", "FRN")
ibddf$NS <- ifelse(grepl(paste(pattern,collapse="|"), ibddf$comp) &  grepl(paste(pattern1,collapse="|"), ibddf$comp), 1, 0)
ibd_sub <- ibddf[ibddf$NS==0,]


plot(1, type="n", xlab="Log Dist", ylab=expression("F"[st]*"/"*"(1-F"[st]*")"), xlim=c(10.5, 15), ylim=c(0, 0.16))
#plot bernerary and tul in same colours
#plot all points first
points(ibd_sub$dgeo, ibd_sub$dgen)
abline(lm(ibd_sub$dgen~ibd_sub$dgeo))

#---------------------------------------------------
#hierfstats
library(hierfstat)
geni <- gl2gi(z)
pop(geni) <-  pops

#Testing for sex bias
#covnert genind to hierfstat data frame
hf<-genind2hierfstat(geni,pop=pop(geni))
sextestVaic<-sexbias.test(hf, z@other$ind.metrics$sex, test="vAIc", nperm=9999)
sextestMaic<-sexbias.test(hf, z@other$ind.metrics$sex, test="mAIc", nperm=9999)
sextestFst<-sexbias.test(hf, z@other$ind.metrics$sex, test="FST", nperm=9999)

#change pops beteeen north and south
#remove frn next time
hf$pop <- ifelse(hf$pop %in% c("BER", "TUL", "DGC", "MLG", "OBN", "RVS"), "North", "South")
sextestVaic<-sexbias.test(hf, z@other$ind.metrics$sex, test="vAIc", nperm=9999)
sextestMaic<-sexbias.test(hf, z@other$ind.metrics$sex, test="mAIc", nperm=9999)
sextestFst<-sexbias.test(hf, z@other$ind.metrics$sex, test="FST", nperm=999)

sexaic<-AIc(hf)
boxplot(sexaic ~ z@other$ind.metrics$sex)
boxplot(aic[which((hf$pop != "TUL") & (hf$pop != "BER"))] ~  z@other$ind.metrics$sex[which((hf$pop != "TUL") & (hf$pop != "BER"))])
var.test(aic ~ z@other$ind.metrics$sex, 
         alternative = "two.sided")

#manual permutation test

test_stat = numeric(length=1000)
for (i in 1:1000 )
{
  perm = sample(z@other$ind.metrics$sex, replace=FALSE)
  test_stat[i] <-sexbias.test(hf, perm, test="vAIc")$statistic
}

#violinplot
aic_p <- data.frame(aic=aic, sex=z@other$ind.metrics$sex, pop=z@other$ind.metrics$pops)
library(ggplot2)
p <- ggplot(aic_p, aes(x=sex, y=aic)) + 
  geom_violin()
p + geom_jitter(shape=19, position=position_jitter(0.2))
#remove hebrides
hfsub <- subset(hf, (hf$pop != "TUL") & (hf$pop != "BER"))
aics<-AIc(hf[which((hf$pop != "TUL") & (hf$pop != "BER")),])
sexbias.test(hf[which((hf$pop != "TUL") & (hf$pop != "BER"))], z@other$ind.metrics$sex[which((hf$pop != "TUL") & (hf$pop != "BER"))])
boxplot(aics ~  z@other$ind.metrics$sex[which((hf$pop != "TUL") & (hf$pop != "BER"))])

#WC fst 
fst_out <- boot.ppfst(dat=hf,diploid=TRUE, nboot=1000)
#using dartR
res<-gl.fst.pop(z, nboots = 1000, nclusters = 2)


#plptting the Fst's
fstMat <- as.matrix(as.dist(res$Fsts)) 
# Use fst between variables as distance to reorder
dd <- as.dist(fstMat)
hc <- hclust(dd)
fstMat <-fstMat[hc$order, hc$order]
#reorder the p-value matrix in the same way
pfstMat <- as.matrix(as.dist(res$Pvalues)) 
pfstMat <-pfstMat[(rownames(fstMat)), (colnames(fstMat))]
#reshape the fst matrix for use with ggplot2
fstMat[lower.tri(fstMat)] <- NA
diag(fstMat) <- NA
pfstMat[lower.tri(pfstMat)] <- NA
diag(pfstMat) <- NA

fstmat <- melt(fstMat, na.rm = TRUE)
names(fstmat) <- c("Pop1", "Pop2", "Fst")
pfstmat <- melt(pfstMat, na.rm=T)
fstmat$pval <- pfstmat$value

# Create a ggheatmap
ggheatmap <- ggplot(fstmat, aes(Pop2, Pop1, fill = Fst))+
  geom_tile(color = "white")+
  scale_fill_distiller(palette = "YlOrRd", direction = 1, name=expression("F"[st])) +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)


ggheatmap + 
  geom_text(aes(Pop2, Pop1, label = ifelse(pval < 0.05/105, round(Fst, 3), "") ), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 13, barheight = 1,
                               title.position = "top", title.hjust = 0.5))




#locus specific fst
locusfst<-Fst(as.loci(geni), pop=pop(geni))



#---------------------------------------------------------------
#identifying sex linked markers
#get sex
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

#are there any markers where female are always homozygotes?

for (i in 1:dim(mat)[2]){
  cross<-xtabs(~mat[,i]+z@other$ind.metrics$sex)
  
  if (dim(cross)[1] == 2){
    if (cross[1,2] == 0)
    {print(cross)}
  }
 
}
#females are never homozygoyous with putative sex linked markers
#look at coverage at some of these markers
#excample
marker_results$locus[marker_results$aic<11]
#get marker row

#based on 25% missing data filter
"""locus      aic no_hets
185    2388_17 8.054117       0 # 3 snps with same strong signal, ETB 1 MLG 1
257    3377_91 8.358014       0 # 1 snp with strong signal, ETB 1 MLG 1, FRN 2
447    5968_27 8.140440       0 # 1 snp with strong signal ETB 1 MLG 1
791   10262_29 9.038965      12 # has 2 snps with inconsistent signal
1109  13518_65 8.361586       0 # 1 snp with strong signal ETB 1 MLG 1, FRN 2
1120  13618_55 8.216726       0 # 1 snp with strong signal ETB 1 MLG 1
2072  26044_63 8.286326       0 # 2 snps with strong consitent signal ETB 1, FRN1, BMD 1, 1 MLG
2497 222956_60 7.958471       0 # 2 snps with inconsistent signal
2671 396798_60 8.595931       0 # 1 snp with some signal ETB 1 BMD1, DGC 1, FRN2, MDC3, MLG1
    13381                       #2 snps with inconsistent signal
    20107                       # 3 snps with inconsistent signal
    25088"""                    # 4 snps with incosnsitent signal

grab=which(VCF@fix[,1] == "13518")
#extract info
gt<- sapply(strsplit(VCF@gt[grab[1],2:length(VCF@gt[grab[1],])], ':'), function(x){paste0(x[1])})
gt[gt=="./."] <- NA
dp <- as.numeric(sapply(strsplit(VCF@gt[grab[2],2:length(VCF@gt[grab[2],])], ':'), function(x){paste0(x[2])}))
ad <- (sapply(strsplit(VCF@gt[grab[1],2:length(VCF@gt[grab[1],])], ':'), function(x){paste0(x[3])}))
ad[ad=="."] <- NA
#get sex and pop
pop <- (sapply(strsplit(names(gt), '_'), function(x){paste0(x[1])}))
sex <- (sapply(strsplit(names(gt), '_'), function(x){paste0(x[2])}))
geno <- data.frame( pop=pop, sex=sex, gt=gt, readepth=dp, alleledp=ad)
averageref <- mean(as.numeric(sapply(strsplit(as.character(geno$alleledp[geno$sex=="f"]), ','), function(x){paste0(x[1])})),  na.rm=T)
averagesnp <- mean(as.numeric(sapply(strsplit(as.character(geno$alleledp[geno$sex=="f"]), ','), function(x){paste0(x[2])})),  na.rm=T)
xtabs(~geno$pop+geno$sex+geno$gt)

markers= c("2388_17", "2388_38" , "2388_43", "3377_91", "5968_27", "13518_65", "13618_55", "26044_63", "26044_76", "396798_60")
sexgts <- mat[,markers]

#repeat analysis with 50% missing data

#----------------------------------------------------------------
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

#----------------------------------------------
#pcadapt
install.packages("pcadapt")
library(pcadapt)

filename <- read.pcadapt("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/populations.snps.filter2.0.25.recode.vcf", type="vcf")
x <- pcadapt(filename, K=10)
plot(x, option = "screeplot")
x <- pcadapt(filename, K=3)
summary(x)
plot(x, option = "scores", i=1,j=3, pop = pops)
#outlier detection using q-values
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha=0.01
outliersPCadapt <- which(qval < alpha)

#--------------------------------------
#pca for males only

#subset genlight object for males only
males <- z[z@other$ind.metrics$sex=="m"]

w <- tab(males, freq = TRUE, NA.method = "mean")

#perform PCA
w.pca <- dudi.pca(w, scannf = F, scale=F, nf=3)

#PCA plot
col <- funky(15)
s.class(w.pca$li, pop(males) ,xax=1,yax=2, col=transp(col,.7), axesell=F,
        cstar=0, cpoint=1.5, clabel=.75, grid=FALSE, addaxes = T) 
title(xlab= "PC1 (7%)", ylab="PC2 (2%)")
