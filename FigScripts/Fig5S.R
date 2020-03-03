library(vcfR)
library(poppr)
library(ape)
library(adegenet)
library(RColorBrewer)
library(treeio)
library(ggtree)
library(tidyverse)
library(ggplot2)
library(dartR)

#Read in the wolbachia data
wol1.VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/Wolbachia/stacks_m4_M4_n4/populations.r50.p3/out.recode.vcf")
#convert to genlight
gl.wol <- vcfR2genlight(wol1.VCF)
#add pop labels
popnames <- gl.wol@ind.names
pops <- sapply(strsplit(popnames, '_'), function(x){paste0(x[1])})
pop(gl.wol) = as.factor(pops)
ploidy(gl.wol) <- 2
#make the distance tree
tree <- aboot(gl.wol, tree = "upgma", distance =bitwise.dist , sample = 100, showtree = F, cutoff = 50, quiet = T)
tree#plot the tree
cols <- brewer.pal(n = nPop(gl.wol), name = "Dark2")
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(gl.wol)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
