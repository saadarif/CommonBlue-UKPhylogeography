library(pophelper)


setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4/populations.r50.p15_moh_0.65/Structure/StructureResults_K1_8/plots/")
sfiles <- list.files(path="/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4/populations.r50.p15_moh_0.65/Structure/StructureResults_K1_8/", pattern="*_f", full.names=T)

slist<-readQ(files=sfiles,filetype="structure")
out<-tabulateQ(slist, writetable=F)
out_s<-summariseQ(out, writetable=F)
evannoMethodStructure(data=out_s,exportplot=T) 

#plot k =2

pfiles <- "/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4/populations.r50.p15_moh_0.65/Structure/StructureResults_K1_8/str_K2_rep1_f"
pops <- read.delim("../../indfile", header=F, stringsAsFactors = F)
names(pops) <- c("ind", "pop", "order")

twolabset1 <- pops[,2:3,drop=FALSE]
twolabset1 <- twolabset1[order(twolabset1$order),]
onelabset1 <- twolabset1[,1,drop=FALSE]
plist <- readQ(files=pfiles,indlabfromfile=T)
#redorder the indiviuals accordign to desired plottign order
plist[[1]]<- plist[[1]][order(pops$order),]
#redorder the plot labeles

plotQ(plist[1],  clustercol=c("skyblue","darkblue"), splab = " K = 2",showlegend=T, ,returnplot=T,exportplot=T,quiet=T,basesize=5,
      grplab=onelabset1, grplabsize=0.9,linesize=0.8,pointsize=1, divsize=0.1)  #legendlab=