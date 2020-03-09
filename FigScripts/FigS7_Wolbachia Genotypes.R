#PCA for Wolbachia genotypes

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

setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/scripts/FigScripts")

VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/Wolbachia//stacks_m4_M4_n4/populations.r50.p3/populations.snps.filter2.0.5.recode.vcf")
z <- vcfR2genlight(VCF)

#get pops and label 
popnames <- z@ind.names
pops <- sapply(strsplit(popnames, '_'), function(x){paste0(x[1])})
pops<- str_replace_all(pops, "RNL", "CFW")
pop(z) = as.factor(pops) 
ploidy(z) <- 1
w <- tab(z, freq = TRUE, NA.method = "mean")

#perform PCA
w.pca <- dudi.pca(w, scannf = F, scale=F, nf=4)

dudi.pca <- tibble(pc1=w.pca$li[,1], pc2=w.pca$li[,2], pc3=w.pca$li[,3], pc4=w.pca$li[,4], pop=pops)
#PCA plot
pc1per <- round(w.pca$eig[1]/sum(w.pca$eig) *100, 1)
pc2per <- round(w.pca$eig[2]/sum(w.pca$eig) *100, 1)
pc3per <- round(w.pca$eig[3]/sum(w.pca$eig) *100, 1)
pc4per <- round(w.pca$eig[4]/sum(w.pca$eig) *100, 1)
col <- c("#A6CEE3", "#4F9F3B","#FDAC4F", "#D1AAB7", "#A99099","#B15928"  )

labels <- dudi.pca %>% distinct(pop, .keep_all=T)

#current plot with pop colours at random
p <- ggplot(dudi.pca,aes(x=pc1, y=pc2, colour=pop, group=pop)) + geom_point(size=2)  + scale_color_manual(name="Pop",values=col) + coord_fixed() +
  xlab(paste0("PC1 ", pc1per, "%")) + ylab(paste0("PC2 ", pc2per, "%")) + theme_bw() + theme(legend.position = "none",panel.grid = element_blank(), axis.text=element_text(size=14, face="bold"),axis.title=element_text(size=14,face="bold"), 
                                                                                             plot.margin = unit(c(1,1,1,1), "lines"))
p <- p + geom_label_repel(data=labels, aes(x=pc1, y=pc2, label=pop), size=3)
p <- p + annotate("text", -4.2, 2, label="MLGm002", size=3)
p <- p + annotate("text", -5.1, 0.9, label="OBNm110", size=3)

#skree plot
skree <- data.frame(PC=1:10, eigen=w.pca$eig[1:10])

sk <- ggplot(skree, aes(x=PC, y=eigen)) + geom_line() + geom_point()+theme_bw() + scale_x_continuous(breaks=1:10) +
  xlab("Number of PCS") + ylab("Eigenvalues") +
  theme(axis.text=element_text(size=14, face="bold"),axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(1,1,1,1), "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


pdf(file = "FigS7.pdf", width = 12, height = 6)
p
dev.off()

