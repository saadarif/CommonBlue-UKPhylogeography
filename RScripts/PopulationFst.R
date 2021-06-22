#Population Fst
#April 2020

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


VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/populations.snps.filter2.0.25.recode.vcf")
z <- vcfR2genlight(VCF)
popnames <- z@ind.names
pops <- sapply(strsplit(popnames, '_'), function(x){paste0(x[1])})
pops<- str_replace_all(pops, "RNL", "CFW")
pop(z) = as.factor(pops) 
ploidy(z) <- 2
popnames <- z@ind.names 
sex<- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[2])}))
ind <- cbind( z@ind.names, pops, sex)
colnames(ind)[1] <- "ind"
z@other$ind.metrics <- ind


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


