#Fig3 PCA for populations of DDRadseq data
#Feb 2020

library("adegenet")
library(vcfR)
library("ade4")
#library(ape)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
require(maps)
library(grid)
library(ggrepel)
library("factoextra")

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
VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/populations.snps.filter2.0.25.recode.singlesnp.vcf")
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

#perform pca with pcaadapt
#library(pcadapt)

#filename <- read.pcadapt("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/populations.snps.filter2.0.25.recode.vcf", type="vcf")
#x <- pcadapt(filename, K=10)
#plot(x, option = "screeplot")
#3 princpal components capture should be enough
#x <- pcadapt(filename, K=3)
#summary(x)
#the third principal component differentiates france from the GB populations
#plot(x, option = "scores", i=1,j=2, pop = pops)

#get data for plot
#pc_data <- tibble(pc1=x$scores[,1], pc2=x$scores[,2], pc3=x$scores[,3],pop=pops)

#ggplot(pc_data,aes(x=pc1, y=pc2, colour=pop)) + geom_point() +stat_ellipse(level=0.9)


## pca using adegent
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
p <- ggplot(dudi.pca,aes(x=pc1, y=pc2, colour=pop, group=pop)) + geom_point() +stat_ellipse(level=0.9) + scale_color_manual(name="Pop",values=col) + coord_fixed(xlim=c(-6, 5), ylim=c(-5,6)) +
  xlab(paste0("PC1 ", pc1per, "%")) + ylab(paste0("PC2 ", pc2per, "%")) + theme_bw() + theme(legend.position = "none",panel.grid = element_blank(), axis.text=element_text(size=14, face="bold"),axis.title=element_text(size=14,face="bold"), 
                                                                                             plot.margin = unit(c(1,1,1,1), "lines"))
p <- p + geom_label_repel(data=labels, aes(x=pc1, y=pc2, label=pop), size=3)

#pop colour by latitude
p8 <- ggplot(dudi.pca,aes(x=pc1, y=pc2, colour=lat, group=pop)) + geom_point() +stat_ellipse(level=0.9) + scale_color_gradient(low = "red", high = "blue", name="Latitude") + coord_fixed(xlim=c(-7, 6), ylim=c(-6,7)) +
  xlab(paste0("PC1 ", pc1per, "%")) + ylab(paste0("PC2 ", pc2per, "%")) + theme_bw() +theme(legend.position = c(0.88,.18),panel.grid = element_blank(), legend.title = element_text(size = 10)) 
p8 <- p8 + geom_label_repel(data=labels, aes(x=pc1, y=pc2, label=pop), size=3) 
p8
#add map

#Get the region
UK <- map_data("world") %>% filter(region=="UK")
#cover coords to data frame
c1 <- as.data.frame(coords)
p1 <- ggplot() +
  geom_polygon(data = UK, aes(x=long, y = lat, group=group), fill="grey",alpha=0.3) +  
  geom_label(data=c1, aes(x=lon, y=lat, label=rownames(c1)), size=2.5, hjust = 0,label.size = 0.15 ) + #scale_color_gradient(low = "red", high = "blue", name="Latitude", guide="none")+
  theme_void() + ylim(50,59) + xlim(NA,2.25) + coord_map() 

#save the plots
vp <- viewport(width = .34, height = .34,  x = 0.5, y = 0.69)
pdf(file = "Fig3.pdf", width = 8, height = 11)
print(p)
print(p1, vp = vp)
dev.off()


#suppl figures
#skree plot
skree <- data.frame(PC=1:10, eigen=w.pca$eig[1:10])

sk <- ggplot(skree, aes(x=PC, y=eigen)) + geom_line() + geom_point()+theme_bw() + scale_x_continuous(breaks=1:10) +
  xlab("Number of PCS") + ylab("Eigenvalues") +
  theme(axis.text=element_text(size=14, face="bold"),axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(1,1,1,1), "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

fviz_eig(w.pca)
#plot pca 3 vs 4

pc34 <- dudi.pca %>% select(c(pc3:pop))

labels34 <- pc34 %>% distinct(pop, .keep_all=T)
p2 <- ggplot(pc34 ,aes(x=pc3, y=pc4, colour=pop, group=pop)) + geom_point() +stat_ellipse(level=0.9) + scale_color_manual(name="Pop",values=col) + coord_fixed() +
  xlab(paste0("PC3 ", pc3per, "%")) + ylab(paste0("PC4 ", pc4per, "%")) + theme_bw() + theme(legend.position = "none",panel.grid = element_blank(), axis.text=element_text(size=14, face="bold"),axis.title=element_text(size=14,face="bold"), 
                                                                                             plot.margin = unit(c(1,1,1,1), "lines"))
p2 <- p2 + geom_label_repel(data=labels34, aes(x=pc3, y=pc4, label=pop), size=3)
p2

pdf(file = "FigS4.pdf", width = 14, height = 8)
ggarrange(sk,p2, labels = c("A", "B"), font.label = list(size = 20, color = "black", face = "bold", family = NULL))
dev.off()



p12 <- ggplot(pc34,aes(x=pc3, y=pc4, colour=lat, group=pop)) + geom_point() +stat_ellipse(level=0.9) + scale_color_gradient(low = "red", high = "blue", name="Latitude") + coord_fixed() +
  xlab(paste0("PC3 ", pc3per, "%")) + ylab(paste0("PC4 ", pc4per, "%")) + theme_bw() +theme(legend.position = c(0.88,.22),panel.grid = element_blank(), legend.title = element_text(size = 10)) 
p12 <- p12 + geom_label_repel(data=labels34, aes(x=pc3, y=pc4, label=pop), size=3) 
p12


#supplemental figure PCA fo PC3 and 4
#PCAs of PC1 AND 2 at different missigness values

#-----------------------------------------------------------------------------------------------------
#supplemental figure S3
VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/populations.snps.filter2.0.5.recode.vcf")
z <- vcfR2genlight(VCF)

#get pops and label 
popnames <- z@ind.names
pops <- sapply(strsplit(popnames, '_'), function(x){paste0(x[1])})
pops<- str_replace_all(pops, "RNL", "CFW")
pop(z) = as.factor(pops) 
ploidy(z) <- 2
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
p50 <- ggplot(dudi.pca,aes(x=pc1, y=pc2, colour=pop, group=pop)) + geom_point() +stat_ellipse(level=0.9) + scale_color_manual(name="Pop",values=col) + coord_fixed(xlim=c(-6, 5), ylim=c(-5,5)) +
  xlab(paste0("PC1 ", pc1per, "%")) + ylab(paste0("PC2 ", pc2per, "%")) + theme_bw() + ggtitle("m:4-M:4-n:4, r-50%, 176 individuals, 5592 SNPs")+
  theme(legend.position = "none",panel.grid = element_blank(), axis.text=element_text(size=14, face="bold"),axis.title=element_text(size=14,face="bold"), 
                                                                                             plot.margin = unit(c(1,1,1,1), "lines"))
p50 <- p50 + geom_label_repel(data=labels, aes(x=pc1, y=pc2, label=pop), size=3)

VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/populations.snps.filter2.0.25.recode.vcf")
z <- vcfR2genlight(VCF)

#get pops and label 
popnames <- z@ind.names
pops <- sapply(strsplit(popnames, '_'), function(x){paste0(x[1])})
pops<- str_replace_all(pops, "RNL", "CFW")
pop(z) = as.factor(pops) 
ploidy(z) <- 2
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
p75 <- ggplot(dudi.pca,aes(x=pc1, y=pc2, colour=pop, group=pop)) + geom_point() +stat_ellipse(level=0.9) + scale_color_manual(name="Pop",values=col) + coord_fixed(xlim=c(-6, 5), ylim=c(-5,5)) +
  xlab(paste0("PC1 ", pc1per, "%")) + ylab(paste0("PC2 ", pc2per, "%")) + theme_bw() + ggtitle("m:4-M:4-n:4, r-50%, 148 individuals, 5592 SNPs")+
  theme(legend.position = "none",panel.grid = element_blank(), axis.text=element_text(size=14, face="bold"),axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(1,1,1,1), "lines"))
p75 <- p75 + geom_label_repel(data=labels, aes(x=pc1, y=pc2, label=pop), size=3)

pdf(file = "FigS3.pdf", width = 14, height = 8)
ggarrange(p50,p75, labels = c("A", "B"), font.label = list(size = 20, color = "black", face = "bold", family = NULL))
dev.off()

#-------------------------------------------
#PCAs testing effects of missing loci

VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r70.p15_moh_0.65/populations.snps.filter2.0.25.recode.vcf")
z <- vcfR2genlight(VCF)

#get pops and label 
popnames <- z@ind.names
pops <- sapply(strsplit(popnames, '_'), function(x){paste0(x[1])})
pops<- str_replace_all(pops, "RNL", "CFW")
pop(z) = as.factor(pops) 
ploidy(z) <- 2
w <- tab(z, freq = TRUE, NA.method = "mean")

#perform PCA
w.pca <- dudi.pca(w, scannf = F, scale=F, nf=4)

dudi.pca <- tibble(pc1=w.pca$li[,1], pc2=w.pca$li[,2], pc3=w.pca$li[,3], pc4=w.pca$li[,4], pop=pops)
#PCA plot
pc1per <- round(w.pca$eig[1]/sum(w.pca$eig) *100, 1)
pc2per <- round(w.pca$eig[2]/sum(w.pca$eig) *100, 1)
pc3per <- round(w.pca$eig[3]/sum(w.pca$eig) *100, 1)
pc4per <- round(w.pca$eig[4]/sum(w.pca$eig) *100, 1)
col <- funky(15)

labels <- dudi.pca %>% distinct(pop, .keep_all=T)

#current plot with pop colours at random
pr70 <- ggplot(dudi.pca,aes(x=pc1, y=pc2, colour=pop, group=pop)) + geom_point() +stat_ellipse(level=0.9) + scale_color_manual(name="Pop",values=col) + coord_fixed(xlim=c(-4, 4), ylim=c(-4,4)) +
  xlab(paste0("PC1 ", pc1per, "%")) + ylab(paste0("PC2 ", pc2per, "%")) + theme_bw() + ggtitle("m:4-M:4-n:4, r-70%, 177 individuals, 742 SNPs")+
  theme(legend.position = "none",panel.grid = element_blank(), axis.text=element_text(size=14, face="bold"),axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(1,1,1,1), "lines"))
pr70 <- pr70 + geom_label_repel(data=labels, aes(x=pc1, y=pc2, label=pop), size=3)

VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r60.p15_moh_0.65/populations.snps.filter2.0.25.recode.vcf")
z <- vcfR2genlight(VCF)

#get pops and label 
popnames <- z@ind.names
pops <- sapply(strsplit(popnames, '_'), function(x){paste0(x[1])})
pops<- str_replace_all(pops, "RNL", "CFW")
pop(z) = as.factor(pops) 
ploidy(z) <- 2
w <- tab(z, freq = TRUE, NA.method = "mean")

#perform PCA
w.pca <- dudi.pca(w, scannf = F, scale=F, nf=4)

dudi.pca <- tibble(pc1=w.pca$li[,1], pc2=w.pca$li[,2], pc3=w.pca$li[,3], pc4=w.pca$li[,4], pop=pops)
#PCA plot
pc1per <- round(w.pca$eig[1]/sum(w.pca$eig) *100, 1)
pc2per <- round(w.pca$eig[2]/sum(w.pca$eig) *100, 1)
pc3per <- round(w.pca$eig[3]/sum(w.pca$eig) *100, 1)
pc4per <- round(w.pca$eig[4]/sum(w.pca$eig) *100, 1)
col <- funky(15)

labels <- dudi.pca %>% distinct(pop, .keep_all=T)

#current plot with pop colours at random
pr60 <- ggplot(dudi.pca,aes(x=pc1, y=pc2, colour=pop, group=pop)) + geom_point() +stat_ellipse(level=0.9) + scale_color_manual(name="Pop",values=col) + coord_fixed(xlim=c(-6, 5), ylim=c(-5,5)) +
  xlab(paste0("PC1 ", pc1per, "%")) + ylab(paste0("PC2 ", pc2per, "%")) + theme_bw() + ggtitle("m:4-M:4-n:4, r-60%, 163 individuals, 2203 SNPs")+
  theme(legend.position = "none",panel.grid = element_blank(), axis.text=element_text(size=14, face="bold"),axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(1,1,1,1), "lines"))
pr60 <- pr60 + geom_label_repel(data=labels, aes(x=pc1, y=pc2, label=pop), size=3)

pdf(file = "FigS3.pdf", width = 14, height = 10)
ggarrange(p50,p75, pr60, pr70, labels = c("A", "B", "C", "D"), ncol=2, nrow=2,font.label = list(size = 20, color = "black", face = "bold", family = NULL))
dev.off()

