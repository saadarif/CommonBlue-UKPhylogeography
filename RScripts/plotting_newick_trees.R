library("treeio")
library("ggtree")
library("adegenet")
library("ape")

setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4/populations.r50.p15_moh_0.65/TreeBuilding/raxml-ng_out_allsnps_FRNoutg///")

tbeTree <-read.newick("T14_tbeSup.raxml.support", node.label='support')


#get pop names
pops <- sapply(strsplit(tbeTree@phylo$tip.label, '_'), function(x){paste0(x[1])})



dat <- data.frame(ind=as.character(tbeTree@phylo$tip.label), pop=pops, stringsAsFactors = F)

dat$colour <- "England&Wales"
dat$colour[dat$pop=="BER" | dat$pop=="TUL"] <- "Hebrides"
dat$colour[dat$pop=="FRN"] <- "France"
dat$colour[dat$pop=="MLG" | dat$pop=="DGC" | dat$pop=="OBN"]  <- "N.Scotland"
dat$colour[dat$pop=="RVS" | dat$pop=="RHD"]  <- "N.England"

populations <- split(dat$ind, dat$colour)

#add group info to tree

p <- ggtree(tbeTree) + theme_tree2() + geom_text2(aes(subset =  (support > .8), label="*"), size=7, col="black", nudge_x = -0.0001, nudge_y = -1, hjust = 1) 
p <- groupOTU(p, populations, 'Location') + aes(color=Location) +  scale_fill_discrete(name = "Location") + scale_color_manual(values = c("grey","black", "red", "khaki", "purple" ), na.value = "black", name = "Location") +theme(legend.position = c(0.1, .8),  legend.box.margin = margin(0, 0, 0, 0, "cm"), legend.key.size = unit(1.2, "cm") )  +
  guides(color = guide_legend(override.aes = list(size = 5, shape = 15)))

#flip branches
flip(p, 219, 153) %>% flip(192, 164)
#collapse a clade 
p <- scaleClade(p, 216, .2) %>% collapse(216, 'min', fill="purple") 
p <- collapse(p, node=214) 
#circuler tree
circ <- ggtree(tbeTree, layout="fan", open.angle=120)+ geom_text2(aes(subset =  (support > .60), label="*"), size=8, col="black", hjust = 1)
groupOTU(circ, populations, 'Pops') + aes(color=Pops) +
  theme(legend.position="right",  legend.key.size = unit(1., "lines") )
