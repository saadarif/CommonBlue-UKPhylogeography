#Fig 5. phyologeny of wolbachia snps against CO1 sequences of given individuals
#Feb 2020 SA
setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/scripts/FigScripts")

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
wol.VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/Wolbachia/stacks_m4_M4_n4/populations.r50.p3/populations.snps.filter2.0.5.recode.vcf")

#convert to genlight
gl.wol <- vcfR2genlight(wol.VCF)



#add pop labels
popnames <- gl.wol@ind.names
pops <- sapply(strsplit(popnames, '_'), function(x){paste0(x[1])})
pop(gl.wol) = as.factor(pops)
ploidy(gl.wol) <- 2

#fixed diffs
fd <- gl.fixed.diff(gl.wol)
#make the distance tree
tree <- aboot(gl.wol, tree = "upgma", distance =bitwise.dist , sample = 1000, showtree = F, cutoff = 50, quiet = T)
tree#plot the tree
cols <- brewer.pal(n = nPop(gl.wol), name = "Dark2")
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(gl.wol)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
#legend('topleft', legend = c("CA","OR","WA"), fill = cols, border = FALSE, bty = "n", cex = 2)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")
#save tree in newick format
write.tree(tree, "woltree.nwk")
#tree for CO1 data
co1 <- read.dna("/media/data_disk/PROJECTS/Saad/CommonBlue/NEWCO1/NewCO1_infected.fasta", format = 'fasta')
#convert to genind
gi.co1 <- DNAbin2genind(co1)
#reset names and pops
indnames<- indNames(gi.co1)
shortname <- sapply(strsplit(indnames, ' '), function(x){paste0(x[1])})
shortname <- sapply(strsplit(shortname, '_'), function(x){paste0(x[1],"_", x[2], "_", x[3])})
#replacenames
indNames(gi.co1) <- shortname
#add pops
pop(gi.co1) <-  sapply(strsplit(shortname, '_'), function(x){paste0(x[1])})
#copnvert to gen light for bitwise dist
gl.co1<- gi2gl(gi.co1)
ploidy(gl.co1) <-2
co1tree <- aboot(gl.co1, tree = "nj", dist=bitwise.dist,sample = 1000, showtree = F, cutoff = 50, quiet = T)
plot.phylo(co1tree, cex = 0.8, font = 2, adj = 0,tip.color =  cols[pop(gi.co1)] )
nodelabels(co1tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
#legend('topleft', legend = c("CA","OR","WA"), fill = cols, border = FALSE, bty = "n", cex = 2)
axis(side = 1)
#title(xlab = "Genetic distance (proportion of loci that are different)")
#save tree in newick format
write.tree(co1tree, "co1tree.nwk")

#plot trees next to each other using ggtree
#convert phylo object to a tibble for use with ggtree
#tbl_woltree <- as_tibble(tree) %>% as.treedata()
#tbl_co1tree <- as_tibble(co1tree) %>% as.treedata()

#read in the saved trees in newick format
wol_tree <- read.newick("woltree.nwk", node.label = "support")
co1_tree <- read.newick("co1tree.nwk", node.label = "support")
p1 <- ggtree(wol_tree)
p2 <- ggtree(co1_tree)

#add popinformation
#get pop names
#get pop names
pops <- sapply(strsplit(wol_tree@phylo$tip.label, '_'), function(x){paste0(x[1])})
wol <- as_tibble(wol_tree)
dat <- tibble(label=as.character(wol_tree@phylo$tip.label), pop=pops) %>% mutate(geog=ifelse(pop=="BER" | pop =="TUL", "OH", "Main"))
y <- full_join(wol, dat, by = 'label')
wol_tree2 <- as.treedata(y)

#see node labels

p1  + geom_label2(aes(subset=!isTip, label=node), size=3, color="darkred", alpha=0.5) 

p1 <- ggtree(wol_tree2) + geom_text2(aes(subset =  (support > 70), label=support), hjust=-0.1)
p1 <- p1 + geom_tippoint(aes(colour=pop, group=pop,size=geog),fill="gray80", shape=19, show.legend =TRUE ) + scale_size_manual(values=c(3.5, 1.5)) +
      scale_colour_manual(values=c("#A6CEE3", "#4F9F3B","#FDAC4F", "#D1AAB7", "#A99099","#B15928"  ))
p1 + geom_treescale()

#construct the CO1 tree
pops <- sapply(strsplit(co1_tree@phylo$tip.label, '_'), function(x){paste0(x[1])})
co1<- as_tibble(co1_tree)
dat2 <- tibble(label=as.character(co1_tree@phylo$tip.label), pop=pops) %>% mutate(geog=ifelse(pop=="BER" | pop =="TUL", "OH", "Main"))
y <- full_join(co1, dat2, by = 'label')
co1_tree2 <- as.treedata(y)

#see node labels

p2  + geom_label2(aes(subset=!isTip, label=node), size=3, color="darkred", alpha=0.5) 
#make the tree
p2 <- ggtree(co1_tree2)  + geom_text2(aes(subset =  (support > 97) & !(isTip), label=support), hjust=1.3,  vjust=-.5)
p2 <- p2 + geom_tippoint(aes(colour=pop, group=pop,size=geog),fill="gray80", shape=19, show.legend =TRUE ) + scale_size_manual(values=c(3.5, 1.5)) +
  scale_colour_manual(values=c("#A6CEE3", "#4F9F3B","#FDAC4F", "#D1AAB7", "#A99099","#B15928"  )) 

  #geom_hilight_encircle(node=33, fill='yellow',alpha=0.6) +
  #geom_hilight_encircle(node=46, fill='purple', alpha=0.6) 
#p2 =revts(p2)



d1 <- p2$data
d2 <- p1$data

## reverse x-axis and 
## set offset to make the tree in the right hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1
#manually adjust tree positions
d2$x <-d2$x - 0.9


pp <- p2  + geom_tree(data=d2)   + geom_text2(aes(subset =  (support > 97) & !(isTip), label=support), data=d2, hjust=-.3,  vjust=-.5)+
  geom_tippoint(aes(size=geog, colour=pop, group=pop),data=d2,fill="gray80", shape=19, show.legend =TRUE )  #+
  #geom_hilight_encircle(node=33, data=d2,fill='yellow',alpha=0.6) +
  #geom_hilight_encircle(node=46, data=d2,fill='purple', alpha=0.6) 
  
dd <- bind_rows(d1, d2) %>% 
  filter(!is.na(label))

pp <- pp + geom_line(aes(x, y, group=label), data=dd, color='grey', alpha=0.5)

#save tree image
pdf(file = "Fig5.pdf", width = 9, height = 11)
#annoatate tree
pp + annotate("text", 0.38,28, hjust=0, size=6, label=expression(atop(italic("Wolbachia ") * "ddRAD Loci", "(37 loci, 68 SNPs)"))) +
  annotate("text", 0.0,19, hjust=0, size=6, label=expression(atop("Mitochondrial CO1", "(16 SNPs)"))) +
  geom_cladelabel(node=45, label="",offset=0.01,align=TRUE, colour="Purple" ,hjust=-0.1, barsize = 2, alpha=0.5)+ 
  geom_cladelabel(node=32,label="",offset=0.01,align=TRUE, colour="Yellow", angle=90 ,hjust=-0.1, barsize = 2, alpha=0.5)+ 
  #geom_cladelabel(node=46, label="Alicante-Sierra Nevada", offset=0.02,align=TRUE, color='black', angle=270, offset.text = 0.01)+ 
  #theme(legend.position = c(0.1,0.85))
  annotate("text", x= -0.01 ,y=42, hjust=0, size=6, label="Population", Parse=TRUE) +
  annotate("point", x= 0.0 ,y=39, hjust=0, size=2, shape=19, colour="#A6CEE3", Parse=TRUE) +
  annotate("text", x= 0.02 ,y=39, hjust=0, size=4, label="BER (Outer Hebrides)", Parse=TRUE) +
  annotate("point", x= 0.0 ,y=36, hjust=0, size=2, shape=19, colour="#B15928", Parse=TRUE) +
  annotate("text", x= 0.02 ,y=36, hjust=0, size=4, label="TUL (Outer Hebrides)", Parse=TRUE) +
  annotate("point", x= 0.0 ,y=33, hjust=0, size=4, shape=19, colour="#4F9F3B", Parse=TRUE) +
  annotate("text", x= 0.02 ,y=33, hjust=0, size=4, label="DGC (Scotland North East)", Parse=TRUE) +
  annotate("point", x= 0.0 ,y=30, hjust=0, size=4, shape=19, colour="#FDAC4F", Parse=TRUE) +
  annotate("text", x= 0.02 ,y=30, hjust=0, size=4, label="MLG (Scotland West Coast)", Parse=TRUE) +
  annotate("point", x= 0.0 ,y=27, hjust=0, size=4, shape=19, colour="#D1AAB7", Parse=TRUE) +
  annotate("text", x= 0.02 ,y=27, hjust=0, size=4, label="OBN (Scotland West Coast)", Parse=TRUE) +
  annotate("point", x= 0.0 ,y=24, hjust=0, size=4, shape=19, colour="#A99099", Parse=TRUE) +
  annotate("text", x= 0.02 ,y=24, hjust=0, size=4, label="RHD (Northern England)", Parse=TRUE) + 
  annotate("text", x=0.23, y=8, angle=270, size =4, label="Iberia-Italy")+
  annotate("text", x=0.23, y=23, angle=270, size =4, label="Alicante")+
  annotate("segment", x = 0.29, xend = 0.29, y = 21, yend = 42, colour="black", size=2, alpha=0.5) +
  annotate("segment", x = 0.29, xend = 0.29, y = 1, yend = 20, colour="black", size=2, alpha=0.5) +
  annotate("text", x=0.28, y=30, angle=90, size =4, label=expression(italic("w")*"Ica2"))+
  annotate("text", x=0.28, y=10, angle=90, size =4, label=expression(italic("w")*"Ica1"))+ 
  geom_treescale(x=0.45, y=36, fontsize=5) 
dev.off()


#----------------------------------------------------------------------------------------------
#supplemental figure
#same as figure 5 but with more samples and more missing data

wol.VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/Wolbachia/stacks_m4_M4_n4/populations.r50.p3/out.recode.vcf")

#convert to genlight
gl.wol <- vcfR2genlight(wol.VCF)
#add pop labels
popnames <- gl.wol@ind.names
pops <- sapply(strsplit(popnames, '_'), function(x){paste0(x[1])})
pop(gl.wol) = as.factor(pops)
ploidy(gl.wol) <- 2
#make the distance tree
tree <- aboot(gl.wol, tree = "upgma", distance =bitwise.dist , sample = 1000, showtree = F, cutoff = 50, quiet = T)
tree#plot the tree
cols <- brewer.pal(n = nPop(gl.wol), name = "Dark2")
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(gl.wol)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
#legend('topleft', legend = c("CA","OR","WA"), fill = cols, border = FALSE, bty = "n", cex = 2)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")
#save tree in newick format
write.tree(tree, "woltree_more.nwk")
#tree for CO1 data
co1 <- read.dna("/media/data_disk/PROJECTS/Saad/CommonBlue/NEWCO1/NewCO1_infected_all.fasta", format = 'fasta')
#convert to genind
gi.co1 <- DNAbin2genind(co1)
#reset names and pops
indnames<- indNames(gi.co1)
shortname <- sapply(strsplit(indnames, ' '), function(x){paste0(x[1])})
shortname <- sapply(strsplit(shortname, '_'), function(x){paste0(x[1],"_", x[2], "_", x[3])})
#replacenames
indNames(gi.co1) <- shortname
#add pops
pop(gi.co1) <-  sapply(strsplit(shortname, '_'), function(x){paste0(x[1])})
#copnvert to gen light for bitwise dist

gl.co1<- gi2gl(gi.co1, parallel = F)

co1tree <- aboot(gl.co1, tree = "upgma", dist=bitwise.dist,sample = 1000, showtree = F, cutoff = 50, quiet = T)
plot.phylo(co1tree, cex = 0.8, font = 2, adj = 0,tip.color =  cols[pop(gi.co1)] )
nodelabels(co1tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
#legend('topleft', legend = c("CA","OR","WA"), fill = cols, border = FALSE, bty = "n", cex = 2)
axis(side = 1)
#title(xlab = "Genetic distance (proportion of loci that are different)")
#save tree in newick format
write.tree(co1tree, "co1tree_more.nwk")

#plot trees next to each other using ggtree
#convert phylo object to a tibble for use with ggtree
#tbl_woltree <- as_tibble(tree) %>% as.treedata()
#tbl_co1tree <- as_tibble(co1tree) %>% as.treedata()

#read in the saved trees in newick format
wol_tree <- read.newick("woltree_more.nwk", node.label = "support")
co1_tree <- read.newick("co1tree_more.nwk", node.label = "support")
p1 <- ggtree(wol_tree)
p2 <- ggtree(co1_tree)

#add popinformation
#get pop names
#get pop names
pops <- sapply(strsplit(wol_tree@phylo$tip.label, '_'), function(x){paste0(x[1])})
wol <- as_tibble(wol_tree)
dat <- tibble(label=as.character(wol_tree@phylo$tip.label), pop=pops) %>% mutate(geog=ifelse(pop=="BER" | pop =="TUL", "OH", "Main"))
y <- full_join(wol, dat, by = 'label')
wol_tree2 <- as.treedata(y)

#see node labels

p1  + geom_label2(aes(subset=!isTip, label=node), size=3, color="darkred", alpha=0.5) 

p1 <- ggtree(wol_tree2) + geom_text2(aes(subset =  (support > 70), label=support), hjust=-0.1)
p1 <- p1 + geom_tippoint(aes(colour=pop, group=pop,size=geog),fill="gray80", shape=19, show.legend =TRUE ) + scale_size_manual(values=c(3.5, 1.5)) +
  scale_colour_manual(values=c("#A6CEE3", "#4F9F3B","#FDAC4F", "#D1AAB7", "#A99099","#B15928"  ))
p1 + geom_treescale()

#construct the CO1 tree
pops <- sapply(strsplit(co1_tree@phylo$tip.label, '_'), function(x){paste0(x[1])})
co1<- as_tibble(co1_tree)
dat2 <- tibble(label=as.character(co1_tree@phylo$tip.label), pop=pops) %>% mutate(geog=ifelse(pop=="BER" | pop =="TUL", "OH", "Main"))
y <- full_join(co1, dat2, by = 'label')
co1_tree2 <- as.treedata(y)

#see node labels

p2  + geom_label2(aes(subset=!isTip, label=node), size=3, color="darkred", alpha=0.5) 
#make the tree
p2 <- ggtree(co1_tree2)  + geom_text2(aes(subset =  (support > 97) & !(isTip), label=support), hjust=1.3,  vjust=-.5)
p2 <- p2 + geom_tippoint(aes(colour=pop, group=pop,size=geog),fill="gray80", shape=19, show.legend =TRUE ) + scale_size_manual(values=c(3.5, 1.5)) +
  scale_colour_manual(values=c("#A6CEE3", "#4F9F3B","#FDAC4F", "#D1AAB7", "#A99099","#B15928"  )) 

#geom_hilight_encircle(node=33, fill='yellow',alpha=0.6) +
#geom_hilight_encircle(node=46, fill='purple', alpha=0.6) 
#p2 =revts(p2)



d1 <- p2$data
d2 <- p1$data

## reverse x-axis and 
## set offset to make the tree in the right hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1
#manually adjust tree positions
d2$x <-d2$x - 0.9


pp <- p2  + geom_tree(data=d2)   + geom_text2(aes(subset =  (support > 70) & !(isTip), label=support), data=d2, hjust=-.3,  vjust=-.5)+
  geom_tippoint(aes(size=geog, colour=pop, group=pop),data=d2,fill="gray80", shape=19, show.legend =TRUE )  #+
#geom_hilight_encircle(node=33, data=d2,fill='yellow',alpha=0.6) +
#geom_hilight_encircle(node=46, data=d2,fill='purple', alpha=0.6) 

dd <- bind_rows(d1, d2) %>% 
  filter(!is.na(label))

pp <- pp + geom_line(aes(x, y, group=label), data=dd, color='grey', alpha=0.5)

#save tree image
pdf(file = "FigS8.pdf", width = 9, height = 11)
#annoatate tree
pp + annotate("text", 0.34,36, hjust=0, size=6, label=expression(atop(italic("Wolbachia ") * "ddRAD Loci", "(37 loci, 74 SNPs)"))) +
  annotate("text", 0.0,19, hjust=0, size=6, label=expression(atop("Mitochondrial CO1", "(16 SNPs)"))) +
  geom_cladelabel(node=46, label="",offset=0.01,align=TRUE, colour="Purple" ,hjust=-0.1, barsize = 2, alpha=0.5)+ 
  geom_cladelabel(node=33,label="",offset=0.01,align=TRUE, colour="Yellow", angle=90 ,hjust=-0.1, barsize = 2, alpha=0.5)+ 
  #geom_cladelabel(node=46, label="Alicante-Sierra Nevada", offset=0.02,align=TRUE, color='black', angle=270, offset.text = 0.01)+ 
  #theme(legend.position = c(0.1,0.85))
  annotate("text", x= -0.01 ,y=55, hjust=0, size=6, label="Population", Parse=TRUE) +
  annotate("point", x= 0.0 ,y=52, hjust=0, size=2, shape=19, colour="#A6CEE3", Parse=TRUE) +
  annotate("text", x= 0.02 ,y=52, hjust=0, size=4, label="BER (Outer Hebrides)", Parse=TRUE) +
  annotate("point", x= 0.0 ,y=49, hjust=0, size=2, shape=19, colour="#B15928", Parse=TRUE) +
  annotate("text", x= 0.02 ,y=49, hjust=0, size=4, label="TUL (Outer Hebrides)", Parse=TRUE) +
  annotate("point", x= 0.0 ,y=46, hjust=0, size=4, shape=19, colour="#4F9F3B", Parse=TRUE) +
  annotate("text", x= 0.02 ,y=46, hjust=0, size=4, label="DGC (Scotland North East)", Parse=TRUE) +
  annotate("point", x= 0.0 ,y=43, hjust=0, size=4, shape=19, colour="#FDAC4F", Parse=TRUE) +
  annotate("text", x= 0.02 ,y=43, hjust=0, size=4, label="MLG (Scotland West Coast)", Parse=TRUE) +
  annotate("point", x= 0.0 ,y=40, hjust=0, size=4, shape=19, colour="#D1AAB7", Parse=TRUE) +
  annotate("text", x= 0.02 ,y=40, hjust=0, size=4, label="OBN (Scotland West Coast)", Parse=TRUE) +
  annotate("point", x= 0.0 ,y=37, hjust=0, size=4, shape=19, colour="#A99099", Parse=TRUE) +
  annotate("text", x= 0.02 ,y=37, hjust=0, size=4, label="RHD (Northern England)", Parse=TRUE) + 
  annotate("text", x=0.23, y=8, angle=270, size =4, label="Iberia-Italy")+
  annotate("text", x=0.23, y=23, angle=270, size =4, label="Alicante")+
  annotate("segment", x = 0.29, xend = 0.29, y = 29, yend = 55, colour="black", size=2, alpha=0.5) +
  annotate("segment", x = 0.29, xend = 0.29, y = 1, yend = 27, colour="black", size=2, alpha=0.5) +
  annotate("text", x=0.28, y=41, angle=90, size =4, label=expression(italic("w")*"Ica2"))+
  annotate("text", x=0.28, y=14, angle=90, size =4, label=expression(italic("w")*"Ica1"))+ 
  geom_treescale(x=0.4, y=42, fontsize=5) 
dev.off()
