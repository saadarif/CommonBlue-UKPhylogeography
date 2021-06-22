#Bayescan, PCADAPT, OutFlank outliers
#April 2020

library("adegenet")
library(vcfR)
library(RColorBrewer)
library(poppr)
library(pegas)
library("factoextra")
library(reshape2)
library("ggplot2")
library(coda)
library(tidyverse)
library(ggrepel)
library(dartR)
library(grid)
library(OutFLANK)
#Descend into appropriate directory

#Outgroup FRN and RVS (4 inds) were removed
setwd('/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new//populations.r50.p15_moh_0.65/BayeScan/p15r50miss25/noFRNRVS/')

#checking for convergence
chain1 <- read.table("bayescanLongRunOD100.sel", header = T)
chain1 <- chain1[-c(1)]

#create mcmc objetc
chain1 <- mcmc(chain1,thin=10)

#checcking convergence
plot(chain1)
summary(chain1)
autocorr.diag(chain1)
effectiveSize(chain1)

#tests of chain convergence
geweke.diag(chain1, frac1=0.1, frac2 = 0.5)
heidel.diag(chain1, eps=.1, pvalue = 0.05)

#multiple chains
chain2 <- read.table("bayescanLongRunOD100_chain2.sel", header = T)
chain2 <- chain2[-c(1)]
chain2 <- mcmc(chain2,thin=10)
combined = mcmc.list(chain1,chain2)
plot(combined)
gelman.diag(combined)
gelman.plot(combined,ask)

#info from http://evomics.org/wp-content/uploads/2016/01/BayeScan_BayeScEnv_exercises.pdf


#read in Fsts and check outliers
bayesFst <-  read.table("bayescanLongRunOD100_fst.txt", header = T)
source("/opt/BayeScan2.1/R functions/plot_R.r")
plot_bayescan(bayesFst)
#filter 1% fdr
outliers <- bayesFst[bayesFst$qval<0.01,]
#check if outliers are same as PCA outliers



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



VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/populations.snps.filter2.0.25.recode.vcf")
z <- vcfR2genlight(VCF)
popnames <- z@ind.names
pops <- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[1])}))
pops<- str_replace_all(pops, "RNL", "CFW")
pop(z) = pops
ploidy(z) <- 2
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

dudi.pca <- tibble(pc1=w.pca$li[,1], pc2=w.pca$li[,2], pc3=w.pca$li[,3], pc4=w.pca$li[,4],lat=z@other$latlong[,2], pop=pops)
#PCA plot
pc1per <- round(w.pca$eig[1]/sum(w.pca$eig) *100, 1)
pc2per <- round(w.pca$eig[2]/sum(w.pca$eig) *100, 1)
pc3per <- round(w.pca$eig[3]/sum(w.pca$eig) *100, 1)
pc4per <- round(w.pca$eig[4]/sum(w.pca$eig) *100, 1)
col <- funky(15)

labels <- dudi.pca %>% distinct(pop, .keep_all=T)

#current plot with pop colours at random
p <- ggplot(dudi.pca,aes(x=pc1, y=pc2, colour=pop, group=pop)) + geom_point() +stat_ellipse(level=0.9) + scale_color_manual(name="Pop",values=col) + coord_fixed(xlim=c(-7, 6), ylim=c(-6.5,6)) +
  xlab(paste0("PC1 ", pc1per, "%")) + ylab(paste0("PC2 ", pc2per, "%")) + theme_bw() + theme(legend.position = "none",panel.grid = element_blank(), axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"), 
                                                                                             plot.margin = unit(c(1,1,1,1), "lines"))
pall <- p + geom_label_repel(data=labels, aes(x=pc1, y=pc2, label=pop), size=3)

library(maps)

UK <- map_data("world") %>% filter(region=="UK")
coords <- as.data.frame(coords)
p1 <- ggplot() +
  geom_polygon(data = UK, aes(x=long, y = lat, group=group), colour="grey",alpha=0.1) +  
  geom_label(data=coords, aes(x=lon, y=lat, label=rownames(coords)), size=3, hjust = 0,label.size =NA, fill=NA) + #scale_color_gradient(low = "red", high = "blue", name="Latitude", guide="none")+
  theme_void() + ylim(50,59) + xlim(NA,3) + coord_map() 

map_grob <- ggplotGrob(p1)
pall <- pall + annotation_custom(grob=map_grob, xmin=-6, xmax=2, ymin=-1, ymax=6.5)
#save the plots
vp <- viewport(width = .4, height = .4,  x = 0.49, y = 0.72)
print(pall)
print(p1, vp = vp)
#snps correlated with pcs
#fviz_pca_var(w.pca,
#             col.var = "contrib", # Color by contributions to the PC
#             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#             repel = TRUE ,    # Avoid text overlapping
#             select.var = list(name = NULL, cos2 = 20, contrib = NULL),
#             col.circle = "grey70"
#

#)




#get genotypes in usable formate quick
mat <- as.matrix(z)
dim(mat)#use row ids of bayes outliers to select the columnss
outliermat <- mat[,as.numeric(rownames(outliers))]
#append marker information to bayescan outlier for sorting
outliers$marker <- colnames(outliermat)
outliers[order(outliers$fst),]
#
xtabs(~mat[,c("21281_19")]+pop(z))

#outlier detetction using PCAdapt
#library(pcadapt)
#filename <- read.pcadapt("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/populations.snps.filter2.0.25.recode.vcf", type="vcf")
#x <- pcadapt(filename, K=10)
#plot(x, option = "screeplot")
#x <- pcadapt(filename, K=2)
#summary(x)
#plot(x, option = "scores", i=1,j=3, pop = pops)
##outlier detection using q-values
#library(qvalue)
#qval <- qvalue(x$pvalues)$qvalues
#alpha=0.05
#outliersPCadapt <- which(qval < alpha)

#outlierPCad <- colnames(mat)[outliersPCadapt]


#--------------------------
#USING OUTFLANK

#remove FRN and RVS
z1 <- z[pop(z)!="FRN" & pop(z)!="RVS",]

#outliers using outFLANK
#make a snpwise fst matirx
#code missing data as 9

mat <- as.matrix(z1)
mat[is.na(mat)] <- 9
my_fst <- MakeDiploidFSTMat(mat, locusNames = colnames(mat), popNames = pop(z1))

#run outflank
OF_res<- OutFLANK(my_fst, NumberOfSamples = 13)
#plot
OutFLANKResultsPlotter(OF_res,withOutliers =F, Zoom=T)
length(OF_res$results$LocusName[OF_res$results$OutlierFlag==T])
OF_out <- as.character(OF_res$results$LocusName[OF_res$results$OutlierFlag==T])
#right right trim
#OF_res2<- OutFLANK(my_fst, RightTrimFraction = 0.1, NumberOfSamples = 13)
#OutFLANKResultsPlotter(OF_res2,withOutliers = F, Zoom=T)
#length(OF_res2$results$LocusName[OF_res2$results$OutlierFlag==T])
#OF2_out <- as.character(OF_res2$results$LocusName[OF_res2$results$OutlierFlag==T])
#check right tail pvals, should be flat but for a bump on selected sites
#hist(OF_res2$results$pvaluesRightTail)

#intersection between and bayescan, outlfank
allU <- intersect(outliers$marker,as.character(OF_res$results$LocusName[OF_res$results$OutlierFlag==T]))
#all unique markers in both
all <- unique(c(outliers$marker, as.character(OF_res$results$LocusName[OF_res$results$OutlierFlag==T])))

#write out the outliers locui to remove from file

out_write <- str_split(all, "_", simplify = T)[,1]
#retain only unique loci
out_write <- unique(out_write)
write.table(out_write, "outlier_locus", quote=F, row.names=F, col.names = F)

#write out the outliers locus and postion to file
out_pos <- data.frame(chr=str_split(all, "_", simplify = T)[,1], pos=str_split(all, "_", simplify = T)[,2])
write.table(out_pos, "outlier_position", quote=F, row.names=F, col.names = F)

#get rid of Outflank NA coding
mat <- as.matrix(z)
#drawing allele frequencies across lats, pops
snps <- data.frame(POP=pop(z))
snps <- cbind(snps, mat[,all])

#reorder pops in lat order
snps$POP <- factor(snps$POP, c("FRN", "ETB", "BWD", "BMD", "PCP", "MMS", "CFW", "MDC", "RHD", "RVS", "OBN", "MLG", "DGC", "BER", "TUL"))

allelefreqs <- data.frame(Pops=levels((snps$POP)))
allelefreqs$Pops <- factor(allelefreqs$Pops, c("FRN", "ETB", "BWD", "BMD", "PCP", "MMS", "CFW", "MDC", "RHD", "RVS", "OBN", "MLG",  "BER",  "DGC","TUL"))
for (i in 2:dim(snps)[2]){
  #get table of allele frequencies

  x<-xtabs(~snps[,i]+snps$POP)
  #pick the most common allele
  #if there are only two classeses
  if (dim(x)[1]==3) { 
    #most common allele
    if (sum(x[1,])>sum(x[3,])){
      mca_ind = 1
      other_ind =3
    }
    else {mca_ind=3; other_ind=1 }
    }
  
  
  if (dim(x)[1]==2){
    print( colnames(snps)[i])
    freq <- (2*x[1,]+x[2,])/(2*(x[1,]+x[2,]))
    allelefreqs$marker <- freq
    names(allelefreqs)[i] <- names(snps)[i]
  }
  else {
  #print( colnames(outliermat)[i])
  freq <- (2*x[mca_ind,]+x[2,])/(2*(x[mca_ind,]+x[2,]+x[other_ind,]))
  allelefreqs$marker <- freq
  names(allelefreqs)[i] <- names(snps)[i]}
}


#plot all the outliers
library(reshape2)
d <- melt(allelefreqs, id.vars = "Pops")
#filter France and other small pops
#d <- d %>% filter(Pops!="FRN" & Pops!="RVS" & Pops!="BWD")
#clines
ggplot(d , aes(x=Pops, y=value, group=variable)) + geom_line(alpha=0.2)
#merge pops
d$region <- "South"
d$region[d$Pops=="BER" | d$Pops=="TUL" | d$Pops=="DGC" | d$Pops=="MLG" | d$Pops=="OBN" |  d$Pops=="RHD"] <- "North"
d %>% group_by(region, variable) %>% summarise(mean=mean(value)) %>% 
ggplot(aes(x=mean, fill=region)) + geom_density(alpha=0.2) #+ geom_histogram()
plot(d$value[d$region=="North"], d$value[d$region=="South"])

#---------------------------------

#hierarchical cluster of outlier alleles across pops
#get martix
hclust_matrix <- allelefreqs %>% select(-Pops) %>% as.matrix() 
#assignrownames
rownames(hclust_matrix) <- allelefreqs$Pops 
#no need to scale calculate dist
locus_dist <- dist(t(hclust_matrix))
#cluster
locus_hclust <- hclust(locus_dist, method="complete")

plot(locus_hclust, labels=F)

#cut the groups
locus_cluster <- cutree(locus_hclust, k=3) %>% enframe() %>% rename(SNP= name, cluster=value)#

#join cluser membership
d <- d %>% rename(SNP=variable)
dclust <- d %>% inner_join(locus_cluster, by = "SNP")

dclust %>% ggplot(aes(Pops, value)) + geom_line(aes(group=SNP), alpha=0.5) + 
  geom_line(stat="summary", fun.y="mean", colour="Red", size=2, aes(group=1)) +
  facet_grid(cols=vars(cluster))
#-------------------------------------------------------------------------------------------------------------------
#PCA of neutral vs outlier

#Use the saved VCF file for this as all neutral rad loci have been remoced from those files
VCFn <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/populations.snps.filter2.0.25.recode.neutral.vcf")
zneut <- vcfR2genlight(VCFn)
popnames <- zneut@ind.names
pops <- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[1])}))
pops<- str_replace_all(pops, "RNL", "CFW")
pop(zneut) = pops
ploidy(zneut) <- 2
#add lat long to gen light object
zneut@other$latlong <- latlong[,2:3]
#should now be able to use the dartR IBD function
#add sex information as well
popnames <- zneut@ind.names 
sex<- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[2])}))
ind <- cbind( zneut@ind.names, pops, sex, latlong[,2:3])
colnames(ind)[1] <- "ind"
zneut@other$ind.metrics <- ind

#Do the pca
w1 <- tab(zneut, freq = TRUE, NA.method = "mean")

#perform PCA
w1.pca <- dudi.pca(w1, scannf = F, scale=F, nf=4)

dudi.pca1 <- tibble(pc1=w1.pca$li[,1], pc2=w1.pca$li[,2], pc3=w1.pca$li[,3], pc4=w1.pca$li[,4],lat=z@other$latlong[,2], pop=pops)
#PCA plot
pc1per <- round(w1.pca$eig[1]/sum(w1.pca$eig) *100, 1)
pc2per <- round(w1.pca$eig[2]/sum(w1.pca$eig) *100, 1)
pc3per <- round(w1.pca$eig[3]/sum(w1.pca$eig) *100, 1)
pc4per <- round(w1.pca$eig[4]/sum(w1.pca$eig) *100, 1)
col <- funky(15)

labels <- dudi.pca1 %>% distinct(pop, .keep_all=T)

#current plot with pop colours at random
p <- ggplot(dudi.pca1,aes(x=pc1, y=pc2, colour=pop, group=pop)) + geom_point() +stat_ellipse(level=0.9) + scale_color_manual(name="Pop",values=col) + coord_fixed(xlim=c(-7, 6), ylim=c(-6.5,6)) +
  xlab(paste0("PC1 ", pc1per, "%")) + ylab(paste0("PC2 ", pc2per, "%")) + theme_bw() + theme(legend.position = "none",panel.grid = element_blank(), axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"), 
                                                                                             plot.margin = unit(c(1,1,1,1), "lines"))
pneut <- p + geom_label_repel(data=labels, aes(x=pc1, y=pc2, label=pop), size=3)



fviz_pca_var(w1.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE ,    # Avoid text overlapping
             select.var = list(name = NULL, cos2 = 20, contrib = NULL),
             col.circle = "grey70"
)


#FST plot for neutral loci
#pca with outlier loci
VCFo <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/populations.snps.filter2.0.25.recode.outlier.vcf")
zout <- vcfR2genlight(VCFo)
popnames <- zout@ind.names
pops <- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[1])}))
pops<- str_replace_all(pops, "RNL", "CFW")
pop(zout) = pops
ploidy(zout) <- 2
#add lat long to gen light object
zout@other$latlong <- latlong[,2:3]
#should now be able to use the dartR IBD function
#add sex information as well
popnames <- zout@ind.names 
sex<- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[2])}))
ind <- cbind( zout@ind.names, pops, sex, latlong[,2:3])
colnames(ind)[1] <- "ind"
zout@other$ind.metrics <- ind
w2 <- tab(zout, freq = TRUE, NA.method = "mean")

#perform PCA
w2.pca <- dudi.pca(w2, scannf = F, scale=F, nf=4)

dudi.pca2 <- tibble(pc1=w2.pca$li[,1], pc2=w2.pca$li[,2], pc3=w2.pca$li[,3], pc4=w2.pca$li[,4],lat=z@other$latlong[,2], pop=pops)
#PCA plot
pc1per <- round(w2.pca$eig[1]/sum(w2.pca$eig) *100, 1)
pc2per <- round(w2.pca$eig[2]/sum(w2.pca$eig) *100, 1)
pc3per <- round(w2.pca$eig[3]/sum(w2.pca$eig) *100, 1)
pc4per <- round(w2.pca$eig[4]/sum(w2.pca$eig) *100, 1)
col <- funky(15)

labels <- dudi.pca2 %>% distinct(pop, .keep_all=T)

#current plot with pop colours at random
p <- ggplot(dudi.pca2,aes(x=pc1, y=pc2, colour=pop, group=pop)) + geom_point() +stat_ellipse(level=0.9) + scale_color_manual(name="Pop",values=col) + coord_fixed(xlim=c(-4, 4), ylim=c(-4,4)) +
  xlab(paste0("PC1 ", pc1per, "%")) + ylab(paste0("PC2 ", pc2per, "%")) + theme_bw() + theme(legend.position = "none",panel.grid = element_blank(), axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"), 
                                                                                             plot.margin = unit(c(1,1,1,1), "lines"))
pout <- p + geom_label_repel(data=labels, aes(x=pc1, y=pc2, label=pop), size=3)

fviz_pca_var(w2.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE ,    # Avoid text overlapping
             select.var = list(name = NULL, cos2 = 20, contrib = NULL),
             col.circle = "grey70"
)


#--------------------------------------------------------------------------------------------------------------------
#FST matrix
#FST plot for neutral loci

#using dartR
resn<-gl.fst.pop(zneut, nboots = 10, nclusters = 6)


#plptting the Fst's
fstMat <- as.matrix(as.dist(resn$Fsts)) 
# Use fst between variables as distance to reorder
dd <- as.dist(fstMat)
hc <- hclust(dd)
fstMat <-fstMat[hc$order, hc$order]
#reorder the p-value matrix in the same way
pfstMat <- as.matrix(as.dist(resn$Pvalues)) 
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
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size=12))+
  coord_fixed()


neutFST <- ggheatmap + 
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
    legend.direction = "horizontal",
    legend.title = element_text( size = 12),
    legend.text = element_text( size = 12))+
  guides(fill = guide_colorbar(barwidth = 15, barheight = 2,
                               title.position = "top", title.hjust = 0.5))



#-----
#FST for outlier loci
reso<-gl.fst.pop(zout, nboots = 10, nclusters = 6)


#plptting the Fst's
fstMat <- as.matrix(as.dist(reso$Fsts)) 
# Use fst between variables as distance to reorder
dd <- as.dist(fstMat)
hc <- hclust(dd)
fstMat <-fstMat[hc$order, hc$order]
#reorder the p-value matrix in the same way
pfstMat <- as.matrix(as.dist(reso$Pvalues)) 
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
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size=12))+
  coord_fixed()


outlierFST <- ggheatmap + 
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
    legend.direction = "horizontal",
    legend.title = element_text( size = 12),
    legend.text = element_text( size = 12))+
  guides(fill = guide_colorbar(barwidth = 15, barheight = 2,
                               title.position = "top", title.hjust = 0.5))




#------------------------------------------------------------------------------------------------------------
#IBD for neutral and outlier loci
ibdneut <- gl.ibd(zneut, permutations=999, plot = T)

#reshape the dist matrices 
library(reshape2)
df <- melt(as.matrix(ibdneut$Dgen), varnames = c("pop1", "pop2"))
colnames(df)[3] <- "dgen"
df2 <- melt(as.matrix(ibdneut$Dgeo), varnames = c("pop1", "pop2"))
colnames(df2)[3] <- "dgeo"
df$dgeo <- df2$dgeo
#remove redudnacies
ibneutdf <- df[as.numeric(df$pop2)>as.numeric(df$pop1),]
#delete duplicate comparisons
ibneutdf$comp <- paste0(ibneutdf$pop1, "-", ibneutdf$pop2)
#add columms for comparisons
#north vs south
pattern=c("TUL", "MLG", "DGC", "OBN", "BER")
pattern1=c("CFW", "MDC", "MMS", "PCP", "BMD", "ETB", "BWD")
ibneutdf$region <- ifelse(grepl(paste(pattern,collapse="|"), ibneutdf$comp) &  grepl(paste(pattern1,collapse="|"), ibneutdf$comp),"North vs. South", "other")
#highlight OH vs mainland comparisons

#pattern=c("MLG", "DGC", "OBN")
#pattern1= c("BER", "TUL")
#ibneutdf$region[grepl(paste(pattern,collapse="|"), ibneutdf$comp) &  grepl(paste(pattern1,collapse="|"), ibneutdf$comp)] <-  "OH-Main"

#French comparisons 
ibneutdf$region[grepl("FRN", ibneutdf$comp)] <-  "France vs. the Rest"

#plot
#plot(1, type="n", xlab="Log Dist", ylab=expression("F"[st]*"/"*"(1-F"[st]*")"), xlim=c(10.5, 15), ylim=c(0, 0.16))
#plot bernerary and tul in same colours
#plot all points first
#points(ibneutdf$dgeo, ibneutdf$dgen)
#abline(lm(ibneutdf$dgen~ibneutdf$dgeo))

#hilight north south compairsons
pattern=c("TUL", "MLG", "DGC", "OBN", "BER")
pattern1=c("CFW", "MDC", "MMS", "PCP", "BMD", "ETB", "BWD", "FRN")

points(ibneutdf$dgeo[grepl(paste(pattern,collapse="|"), ibneutdf$comp) &  grepl(paste(pattern1,collapse="|"), ibneutdf$comp)],ibneutdf$dgen[grepl(paste(pattern,collapse="|"), ibneutdf$comp) &  grepl(paste(pattern1,collapse="|"), ibneutdf$comp)], pch=19,col="black" )


#highlight OH vs mainland comparisons

pattern=c("MLG", "DGC", "OBN")
pattern1= c("BER", "TUL")
points(ibneutdf$dgeo[grepl(paste(pattern,collapse="|"), ibneutdf$comp) &  grepl(paste(pattern1,collapse="|"), ibneutdf$comp)],ibneutdf$dgen[grepl(paste(pattern,collapse="|"), ibneutdf$comp) &  grepl(paste(pattern1,collapse="|"), ibneutdf$comp)], pch=19,col="grey80" )

points(ibneutdf$dgeo[ grepl("TUL", ibneutdf$comp) & grepl("BER", ibneutdf$comp)],ibneutdf$dgen[grepl("TUL", ibneutdf$comp) & grepl("BER", ibneutdf$comp)], pch=19,col="red" )



#----
ibdout <- gl.ibd(zout, permutations=2, plot=T)

#reshape the dist matrices 
df <- melt(as.matrix(ibdout$Dgen), varnames = c("pop1", "pop2"))
colnames(df)[3] <- "dgen"
df2 <- melt(as.matrix(ibdout$Dgeo), varnames = c("pop1", "pop2"))
colnames(df2)[3] <- "dgeo"
df$dgeo <- df2$dgeo
df[which(df$dgen>0.1 & df$dgeo>13.0),]
#remove redudnacies
iboutdf <- df[as.numeric(df$pop2)>as.numeric(df$pop1),]
#delete duplicate comparisons
iboutdf$comp <- paste0(iboutdf$pop1, "-", iboutdf$pop2)
iboutdf$region <- ibneutdf$region


#join the two dfs 
ibneutdf$data <- "neutral"
iboutdf$data <- "outlier"
ibd <- rbind(ibneutdf, iboutdf)


ibdplots <- ggplot(ibd, aes(x=dgeo, y=dgen, group=data, shape=region, colour=data)) + geom_point(size=3, alpha=0.6) + theme_bw() + geom_smooth(method="lm", se=F, weight=0.5)+
    xlab("Log Distance")+ ylab(expression("F"[st]*"/"*"(1-F"[st]*")"))+ ylim(c(-.1,1.3)) + scale_colour_manual(values=c("grey", "black")) +
    labs(colour="SNP Type", shape="Regional Comparison")+
  theme(legend.position=c(0.2,0.75),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(1,1,1,1), "lines"),
        legend.title = element_text( size = 10),
        #legend.text = element_text( size = 9),
        #legend.spacing.y = unit(-0.5, "mm"),
        #legend.key.size = unit(2,"line"),
        #legend.text = element_text(margin = margin(r = 8, unit = "pt")),
        #legend.box.background = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



ggarrange(pall,pneut, pout, ibdplots, labels = c("A", "B", "C", "D"), ncol=2, nrow=2,font.label = list(size = 20, color = "black", face = "bold", family = NULL)) %>% 
  ggexport(filename="test.pdf", height =12, width=16)

outlierfst <- ibd%>% filter(data=="neutral")


ggplot(outlierfst, aes(x=dgeo, y=dgen,  shape=region)) + theme_bw() +
  xlab("Log Distance")+ ylab(expression("F"[st]*"/"*"(1-F"[st]*")"))+ ylim(c(-.01,.15)) +
  labs( shape="Regional Comparison")+ scale_shape_manual(values=c(15,16,17))+
  theme(legend.position=c(0.2,0.75),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(1,1,1,1), "lines"),
        legend.title = element_text( size = 10),
        #legend.text = element_text( size = 9),
        #legend.spacing.y = unit(-0.5, "mm"),
        #legend.key.size = unit(2,"line"),
        #legend.text = element_text(margin = margin(r = 8, unit = "pt")),
        #legend.box.background = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + geom_density_2d(na.rm=T, size=0.25) + #stat_density_2d(geom="polygon", aes(fill= ..level..)) +
        geom_point(size=3, alpha=0.6, colour="black")  +
        geom_smooth(method="lm", aes(group=1), se=F, colour="black", na.rm=T, size=1)


#"4155_29"   "4155_42"   "19539_66"  "20988_40"  "21281_19"  "34704_16"  "363255_87"
