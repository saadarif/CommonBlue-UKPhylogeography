#Fig X Beast Time tree
library(tidyverse)
library(treeio)
library(ggplot2)
library(ggtree)
library(grid)
library(scatterpie)
x <- read.beast("/media/data_disk/PROJECTS/Saad/CommonBlue/CO1_BOLD/BEAST/modTest_CCS_sc_split_subst/combined.mcc.tree")
setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/scripts/FigScripts")

#Draw tree for editing for final figure

#rename taxa only to have smaller names
tipnames <- x@phylo$tip.label
boldids <- sapply(strsplit(tipnames, "|", fixed=TRUE),  "[", 1)
country <- sapply(strsplit(tipnames, "|", fixed=TRUE),  "[", 2)
#replace errror and correct formatting
country<- str_replace_all(country, "Polyommatus_icarus", "Italy")
country<- str_replace_all(country, "United_Kingdom", "United Kingdom")
d <- data.frame(label = tipnames,
                boldid = boldids, 
                country= country)
x2 <- full_join(x, d, by = "label")
#editing tree data
#remove hpd ranges for tips
#x@data[isTip(x, x@data$node),]$height_0.95_HPD <- NA
#tree plus time scale + units + space for tip labels

p= ggtree(x2) +geom_treescale(x=-4, y=120, fontsize=3)  + scale_x_continuous(name="MYA", labels = abs) + coord_cartesian(clip = 'off') + theme_tree2(plot.margin=margin(6, 120, 6, 6)) 
p1=revts(p)
#add uncertainty of time
p2= p1+geom_range('height_0.95_HPD', color='red', size=2, alpha=.3)
#add posterior probaility as number +geom_nodelab(aes(x=branch, label=round(posterior*100, 0)), vjust=-.5, size=3, nudge_x = 0.1)
p3= p2  + geom_tiplab(aes(label=paste(boldids,country)),size=2)
#add node labels for further editing
p3 + geom_label2(aes(subset=!isTip, label=node), size=3, color="darkred", alpha=0.5) 

#for removing bars for with support less than 0.7
x2@data$height_0.95_HPD[x2@data$posterior < 0.7] <- NA
#x2@data$height_0.95_HPD[x2@data$node==281] <- NA

#p= ggtree(x2) +geom_treescale(x=-4, y=100, fontsize=3)  + scale_x_continuous(name="MYA", labels = abs) + coord_cartesian(clip = 'off') + theme_tree2(plot.margin=margin(6, 6, 6, 6) ) 
p1 = ggtree(x2) + geom_treescale(x=-1.5, y=100, fontsize=4)  + scale_x_continuous(name="MYA", breaks=c(0,-0.5,-1,-1.5,-2,-2.5,-3),labels = abs, limits=c(NA,0.5)) + theme_tree2(plot.margin=margin(2, 2, 2, 2) )  
p1=revts(p1)
#add uncertainty of time
p2= p1+geom_range('height_0.95_HPD', color='red', size=2, alpha=.3)
#add posterior probaility as number
#p3= p2 +geom_nodelab(aes(x=branch, subset=!isTip & posterior > 0.7, label=round(posterior*100, 0)),hjust=1.5 size=2, nudge_x = 0.08) #+ geom_tiplab(aes(label=country),size=2)
p3 = p2 + geom_text2(aes(label=round(posterior*100, 0), subset=posterior>0.7), size=2.5, nudge_x = -0.07, vjust=-.3)

#flip nodes for arrangment
p3<- ggtree::flip(p3,189,216)  
p3<- ggtree::flip(p3,298,217)   
#p3<- ggtree::rotate(p3, 189)
#p3<- ggtree::rotate(p3, 298)
#p3<- ggtree::rotate(p3, 217)
#add more palatable tip labels
p3 + geom_tiplab(aes(label=paste(x2@extraInfo$boldid, x2@extraInfo$country)), size=1.5)


#p3 <- collapse(p3, 218, 'max', fill="yellow", col="black", alpha=0.25) 
#p3 <- collapse(p3, 281, 'max', fill="yellow", col="black", alpha=0.25) 
#p3 <- collapse(p3, 298, 'max', fill="green", col="black", alpha=0.25) 
#p3 <- collapse(p3, 195, 'max', fill="purple", col="black", alpha=0.25)


#add annoated to clade colours

p6 <- p3 + annotate("text", -3.6, c(180,176), hjust=0, fontface=2,size=4, label=c(expression(bold("CO1 Lineage based on")), expression(bold("Dinca " * italic('et al.') * " (2011)"))))+
  annotate("point",-3.5, 170, size=4, shape=21,fill="yellow",color="black") +
  annotate("point",-3.5, 160, size=4, shape=21,fill="green", color="black") +
  annotate("point",-3.5, 150, size=4, shape=21,fill="purple",color="black") +
  annotate("point",-3.5, 140, size=4, shape=21,fill="red",color="black") +
  annotate("text", -3.3,170, hjust=0, size=3.5, label="Iberia-Italy") + 
  annotate("text", -3.3,160, hjust=0, size=3.5,  label= "Palaeartic" )+
  annotate("text", -3.3,150, hjust=0, size=3.5, label="Alicante-Provence")+
  annotate("text", -3.3,140, hjust=0, size=3.5, label="Sierra Nevada")

#add clade labels with country locations
ibit= c(paste0("Italy"), paste0("Spain"),  paste0("Portugal"), paste0("S. France"))
ibitUK= c(paste0("Great Britain:"), paste0("West Scotland"),  paste0("Northern England"), paste0("South Central England"),paste0("South West England"), paste0("Finland"))

palaeartic=c(paste0("France, Germany"), paste0("Romania, Finland"), paste0("Norway, Netherlands"), paste0("Iran, Turkey"), paste0("Austria, Canada"),
             paste0("Spain, Italy"), paste0("Greece"))
palaearticUK= c(paste0("Great Britain:"), paste0("Wales"), paste0("South East England"), paste0("South West England"), paste0("South Central England") )

apsnc= c(paste0("Great Britain:"), paste0("Outer Hebrides"), paste0("Inner Hebrides"), paste0("West Scotland"), paste0("Alicante, Germany"), paste0("Norway, Kazakhastan"))

p7 <- p6  +  annotate("text", x= c(.05,.05,.05,.05) ,y=c(50,47,44,41), hjust=0, size=2.5, label=ibit, Parse=TRUE) +
      annotate("text", x= c(.05,0.09,0.09,0.09,.09,.05) ,y=c(81,78,75,72, 69, 66), hjust=0, size=2.5, label=ibitUK, Parse=TRUE) +
      annotate("text", x= 0.05 ,y=c(150,147,144,141,138,135,132), hjust=0, size=2.5, label=palaeartic, Parse=TRUE) +
      annotate("text", x= c(.05,0.09,.09,0.09,.09) ,y=c(120,117,114,111,108), hjust=0, size=2.5, label=palaearticUK, Parse=TRUE) +
      annotate("text", x= 0.05 ,y=187, hjust=0, size=2.5, label="Crete", Parse=TRUE) +
      annotate("text", x= 0.05 ,y=184, hjust=0, size=2.5, label="Sierra Nevada", Parse=TRUE) +
      annotate("text", x= c(.05,.09,.09,.09,.05, 0.05) ,y=c(178,175,172,169,166,163), hjust=0, size=2.5, label=apsnc, Parse=TRUE)
p7
p7 <- p7 + ggtitle("A") + theme(plot.title = element_text(size= 20, hjust=0.05, vjust=-6, margin = margin(b = 0, t = -1, l = 2, unit = "cm")))

P8 <- p7  + geom_point2(aes(subset=node==217), fill='yellow', colour="black", pch=21, size=3) +  geom_point2(aes(subset=node==298), fill='green', colour="black", pch=21, size=3)+
    geom_point2(aes(subset=node==195), fill='purple', colour="black", pch=21, size=3) +  geom_point2(aes(subset=node==191), fill='red', colour="black", pch=21, size=3)


#plotting BOLD sample lineages (based on Dinca et al on map)
bold <- read.csv("../../CO1_BOLD/trimmed_alignments/specimen_info.txt", header=T, na.strings=c("", "NA"), row.names = NULL, sep="\t")
#Keep only UK samples for plotting
UKdata1 <- bold %>% filter(country=="United Kingdom")
#group by lineage and then by lat long
UKdata <- UKdata1 %>% group_by(Lineage.Dinca.et.al..2011., lat,lon) %>%mutate(n= n()) %>% distinct( Lineage.Dinca.et.al..2011.,.keep_all=TRUE)


# Get the world polygon and extract UK
library(maps)
UK <- map_data("world") %>% filter(region=="UK")


p10 <- ggplot() +
  geom_polygon(data = UK, aes(x=long, y = lat, group = group), fill="grey", alpha=0.25) +theme_void() + ylim(50,59)   +
  geom_jitter( data=UKdata, aes(x=lon, y=lat, fill=Lineage.Dinca.et.al..2011., size=n) ,position = position_jitter(width = 0.2, height = 0.2), alpha=0.5, colour="black", pch=21) +
 theme_void() + ylim(50,59) + coord_map() + scale_fill_manual(values=c("purple", "yellow", "green")) + ggtitle("B") + scale_size_continuous(name="No. Specimens", range=c(3,7), breaks=c(1,3,5))+ 
 guides(fill=FALSE) + theme( legend.position = c(0.84, 0.8),legend.title=element_text(size=8),legend.text = element_text(size=7), plot.background = element_rect( colour = "black", size=1),
 plot.title = element_text(size= 20, hjust=0.05, vjust=-7, margin = margin(b = -0.1, t = -1, l = 2, unit = "cm")) ) 

    
#Plotting map as inset to tree
vp <- viewport(width = .5, height = .5,  x = 0.27, y = 0.3)
pdf(file = "Fig2AB.pdf", width = 8, height = 11)
print(P8)
print(p10, vp = vp)
dev.off()

#grouping by lat lon to count lineaage for piecharts
ukdataPie<-UKdata1 %>% group_by(Lineage.Dinca.et.al..2011.,lat,lon) %>% arrange(lat)  %>%  summarise(n=n())%>% spread(Lineage.Dinca.et.al..2011.,n)
names(ukdataPie) <- c("lat", "lon", "A", "IB", "P")
ukdataPie[is.na(ukdataPie)] <- 0
ukdataPie$radius <- rowSums(ukdataPie[,c(3:5)])


#plotting pie charts instead of bubbles
#not coord equal or coord fixed perserves the aspect ratio for the pie but alters the map

uk <- map_data("world") %>% filter(region=="UK")
ggplot(uk, aes(long, lat)) +
  geom_map(map=uk, aes(map_id=region), fill="grey97", color="grey") + ylim(50,59)   +
  geom_scatterpie(aes(x=lon, y=lat, r=radius/10), data=ukdataPie, cols=c("A", "IB" , "P"), color="black", alpha=.4 ) + coord_equal() +
  scale_fill_manual(values = c("A" = "Purple", "IB" = "yellow","P" = "green")) + guides(fill=F) + 
  theme_void() +
  geom_scatterpie_legend(ukdataPie$radius/10, 0,58, n=5 , labeller = function(x) x=sort(unique(ukdataPie$radius)))

df <- tibble(
  g = c(1, 1, 2, 2),
  x = c(1, 1, 2, 1)
) %>% group_by(g)
df %>% distinct()
df %>% distinct(x)


#----------------------------------------------------------------------------------------------------------
#Alternate layout for Fig2


p6 <- p3 + annotate("text", -3, c(90), hjust=0, fontface=2,size=4, label=c(expression(bold("CO1 Lineage based on ")* bold("Dinca " * italic('et al.') * " (2011)"))))+
  annotate("point",-2.9, 85, size=4, shape=21,fill="yellow",color="black") +
  annotate("point",-2.9, 75, size=4, shape=21,fill="green", color="black") +
  annotate("point",-2.9, 65, size=4, shape=21,fill="purple",color="black") +
  annotate("point",-2.9, 55, size=4, shape=21,fill="red",color="black") +
  annotate("text", -2.7,85, hjust=0, size=3.5, label="Iberia-Italy") + 
  annotate("text", -2.7,75, hjust=0, size=3.5,  label= "Palaeartic" )+
  annotate("text", -2.7,65, hjust=0, size=3.5, label="Alicante-Provence")+
  annotate("text", -2.7,55, hjust=0, size=3.5, label="Sierra Nevada")

#add clade labels with country locations
ibit= c(paste0("Italy"), paste0("Spain"),  paste0("Portugal"), paste0("S. France"))
ibitUK= c(paste0("Great Britain:"), paste0("West Scotland"),  paste0("Northern England"), paste0("South Central England"),paste0("South West England"), paste0("Finland"))

palaeartic=c(paste0("France, Germany"), paste0("Romania, Finland"), paste0("Norway, Netherlands"), paste0("Iran, Turkey"), paste0("Austria, Canada"),
             paste0("Spain, Italy"), paste0("Greece"))
palaearticUK= c(paste0("Great Britain:"), paste0("Wales"), paste0("South East England"), paste0("South West England"), paste0("South Central England") )

apsnc= c(paste0("Great Britain:"), paste0("Outer Hebrides"), paste0("Inner Hebrides"), paste0("West Scotland"), paste0("Alicante, Germany"), paste0("Norway, Kazakhastan"))

p7 <- p6  +  annotate("text", x= c(.05,.05,.05,.05) ,y=c(50,47,44,41), hjust=0, size=2.5, label=ibit, Parse=TRUE) +
  annotate("text", x= c(.05,0.09,0.09,0.09,.09,.05) ,y=c(81,78,75,72, 69, 66), hjust=0, size=2.5, label=ibitUK, Parse=TRUE) +
  annotate("text", x= 0.05 ,y=c(150,147,144,141,138,135,132), hjust=0, size=2.5, label=palaeartic, Parse=TRUE) +
  annotate("text", x= c(.05,0.09,.09,0.09,.09) ,y=c(120,117,114,111,108), hjust=0, size=2.5, label=palaearticUK, Parse=TRUE) +
  annotate("text", x= 0.05 ,y=187, hjust=0, size=2.5, label="Crete", Parse=TRUE) +
  annotate("text", x= 0.05 ,y=184, hjust=0, size=2.5, label="Sierra Nevada", Parse=TRUE) +
  annotate("text", x= c(.05,.09,.09,.09,.05, 0.05) ,y=c(178,175,172,169,166,163), hjust=0, size=2.5, label=apsnc, Parse=TRUE)
#p7
#p7 <- p7 + ggtitle("A") + theme(plot.title = element_text(size= 20, hjust=0.05, vjust=-6, margin = margin(b = 0, t = -1, l = 2, unit = "cm")))

P8 <- p7  + geom_point2(aes(subset=node==217), fill='yellow', colour="black", pch=21, size=3) +  geom_point2(aes(subset=node==298), fill='green', colour="black", pch=21, size=3)+
  geom_point2(aes(subset=node==195), fill='purple', colour="black", pch=21, size=3) +  geom_point2(aes(subset=node==191), fill='red', colour="black", pch=21, size=3)

p10 <- ggplot() +
  geom_polygon(data = UK, aes(x=long, y = lat, group = group), fill="grey", alpha=0.25) +theme_void() + ylim(50,59)   +
  geom_jitter( data=UKdata, aes(x=lon, y=lat, fill=Lineage.Dinca.et.al..2011., size=n) ,position = position_jitter(width = 0.2, height = 0.2), alpha=0.5, colour="black", pch=21) +
  theme_void() + ylim(50,59) + coord_map() + scale_fill_manual(values=c("purple", "yellow", "green"))  + scale_size_continuous(name="No. Specimens", range=c(3,7), breaks=c(1,3,5))+ 
  guides(fill=FALSE) + theme( legend.position = c(0.84, 0.8),legend.title=element_text(size=10),legend.text = element_text(size=8), plot.background = element_rect( colour = "black", size=1),
                              plot.title = element_text(size= 20, hjust=0.05, vjust=-7, margin = margin(b = -0.1, t = -1, l = 2, unit = "cm")) ) 

pdf(file = "Fig2Alt.pdf", width = 14, height = 8)
ggarrange(P8,p10, labels = c("A", "B"), font.label = list(size = 20, color = "black", face = "bold", family = NULL))
dev.off()

##---------------------------------------------------------------------------------------------
#layout with tip labels
library(tidyverse)
library(treeio)
library(ggplot2)
library(ggtree)
library(grid)
library(scatterpie)
x <- read.beast("/media/data_disk/PROJECTS/Saad/CommonBlue/CO1_BOLD/BEAST/modTest_CCS_sc_split_subst/combined.mcc.tree")
setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/scripts/FigScripts")
#rename taxa only to have smaller names
tipnames <- x@phylo$tip.label
boldids <- sapply(strsplit(tipnames, "|", fixed=TRUE),  "[", 1)
country <- sapply(strsplit(tipnames, "|", fixed=TRUE),  "[", 2)
#replace errror and correct formatting
country<- str_replace_all(country, "Polyommatus_icarus", "Italy")
country<- str_replace_all(country, "United_Kingdom", "United Kingdom")
d <- data.frame(label = tipnames,
                boldid = boldids, 
                country= country)
x2 <- full_join(x, d, by = "label")

#for removing bars for with support less than 0.8
x2@data$height_0.95_HPD[x2@data$posterior < 0.8] <- NA

p1 = ggtree(x2) + geom_treescale(x=-1.5, y=100, fontsize=4)  + scale_x_continuous(name="MYA", breaks=c(0,-0.5,-1,-1.5,-2,-2.5,-3),labels = abs, limits=c(NA,0.5)) + theme_tree2(plot.margin=margin(2, 2, 2, 2) )  
p1=revts(p1)
#add uncertainty of time
p2= p1+geom_range('height_0.95_HPD', color='black', size=2, alpha=.25)
#add posterior probaility as number
#p3= p2 +geom_nodelab(aes(x=branch, subset=!isTip & posterior > 0.7, label=round(posterior*100, 0)),hjust=1.5 size=2, nudge_x = 0.08) #+ geom_tiplab(aes(label=country),size=2)
p3 = p2 + geom_text2(aes(label=round(posterior, 1), subset=posterior>0.8), size=3, nudge_x = -0.05, vjust=-.3)

#flip nodes for arrangment
p3<- ggtree::flip(p3,189,216)  
p3<- ggtree::flip(p3,217,317)
p3<- ggtree::rotate(p3, 189) 

p3 <- p3 + geom_tiplab(aes(label=paste(x2@extraInfo$boldid, x2@extraInfo$country)), size=1.75)


#add annoated to clade colours

p6 <- p3 + annotate("text", -2.6, c(180), hjust=0, fontface=2,size=4, label=c(expression(bold("CO1 Lineage"))))+
  annotate("point",-2.5, 170, size=4, shape=21,fill="yellow",color="black") +
  annotate("point",-2.5, 160, size=4, shape=21,fill="green", color="black") +
  annotate("point",-2.5, 150, size=4, shape=21,fill="purple",color="black") +
  annotate("point",-2.5, 140, size=4, shape=21,fill="red",color="black") +
  annotate("text", -2.3,170, hjust=0, size=3.5, label="Iberia-Italy") + 
  annotate("text", -2.3,160, hjust=0, size=3.5,  label= "Palaeartic" )+
  annotate("text", -2.3,150, hjust=0, size=3.5, label="Alicante")+
  annotate("text", -2.3,140, hjust=0, size=3.5, label="Sierra Nevada")

P8 <- p6 + geom_hilight(218,fill = "yellow", alpha = 0.25) + geom_hilight(317,fill = "green", alpha = 0.25) + 
        geom_hilight(195,fill = "purple", alpha = 0.25) + geom_hilight(191,fill = "red", alpha = 0.25) 

P8 <- P8 + ggtitle("A") + theme(plot.title = element_text(size= 20, hjust=0.05, vjust=-6, margin = margin(b = 0, t = -1, l = 2, unit = "cm")))

P8

#plotting BOLD sample lineages (based on Dinca et al on map)
bold <- read.csv("../../CO1_BOLD/trimmed_alignments/specimen_info.txt", header=T, na.strings=c("", "NA"), row.names = NULL, sep="\t")
#Keep only UK samples for plotting
UKdata1 <- bold %>% filter(country=="United Kingdom")
#group by lineage and then by lat long
UKdata <- UKdata1 %>% group_by(Lineage.Dinca.et.al..2011., lat,lon) %>%mutate(n= n()) %>% distinct( Lineage.Dinca.et.al..2011.,.keep_all=TRUE)


# Get the world polygon and extract UK
library(maps)
UK <- map_data("world") %>% filter(region=="UK")


p10 <- ggplot() +
  geom_polygon(data = UK, aes(x=long, y = lat, group = group), fill="grey", alpha=0.3) +theme_void() + ylim(50,59)   +
  geom_jitter( data=UKdata, aes(x=lon, y=lat, fill=Lineage.Dinca.et.al..2011., size=n) ,position = position_jitter(width = 0.2, height = 0.2), alpha=0.5, colour="black", pch=21) +
  theme_void() + ylim(50,59) + coord_map() + scale_fill_manual(values=c("purple", "yellow", "green")) + ggtitle("B") + scale_size_continuous(name="No. Specimens", range=c(3,7), breaks=c(1,3,5))+ 
  guides(fill=FALSE) + theme( legend.position = c(0.84, 0.8),legend.title=element_text(size=8),legend.text = element_text(size=7), plot.background = element_rect( colour = "black", size=1),
                              plot.title = element_text(size= 20, hjust=0.05, vjust=-7, margin = margin(b = -0.1, t = -1, l = 2, unit = "cm")) ) 


#Plotting map as inset to tree
vp <- viewport(width = .35, height = .35,  x = 0.24, y = 0.22)
pdf(file = "Fig2AB.pdf", width = 8, height = 12)
print(P8)
print(p10, vp = vp)
dev.off()


#---------------------------------------------------------------------------------------------
#final alternate


library(tidyverse)
library(treeio)
library(ggplot2)
library(ggtree)
library(grid)
library(scatterpie)
x <- read.beast("/media/data_disk/PROJECTS/Saad/CommonBlue/CO1_BOLD/BEAST/modTest_CCS_sc_split_subst/combined.mcc.tree")
setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/scripts/FigScripts")
#rename taxa only to have smaller names
tipnames <- x@phylo$tip.label
boldids <- sapply(strsplit(tipnames, "|", fixed=TRUE),  "[", 1)
country <- sapply(strsplit(tipnames, "|", fixed=TRUE),  "[", 2)
#replace errror and correct formatting
country<- str_replace_all(country, "Polyommatus_icarus", "Italy")
country<- str_replace_all(country, "United_Kingdom", "United Kingdom")
d <- data.frame(label = tipnames,
                boldid = boldids, 
                country= country, stringsAsFactors = F)
d$country[d$boldid=="EULEP004-14"] <- "Crete" 

x2 <- full_join(x, d, by = "label")

#for removing bars for with support less than 0.8
x2@data$height_0.95_HPD[x2@data$posterior < 0.75] <- NA

p1 = ggtree(x2) + geom_treescale(x=-1.5, y=100, fontsize=4)  + scale_x_continuous(name="MYA", breaks=c(0,-0.5,-1,-1.5,-2,-2.5,-3),labels = abs, limits=c(NA,0.5)) + theme_tree2(plot.margin=margin(2, 2, 2, 2) )  
p1=revts(p1)
#add uncertainty of time
p2= p1+geom_range('height_0.95_HPD', color='black', size=2, alpha=.25)
#add posterior probaility as number
#p3= p2 +geom_nodelab(aes(x=branch, subset=!isTip & posterior > 0.7, label=round(posterior*100, 0)),hjust=1.5 size=2, nudge_x = 0.08) #+ geom_tiplab(aes(label=country),size=2)
p3 = p2 + geom_text2(aes(label=round(posterior, 1), subset=posterior>0.75), size=3, nudge_x = -0.05, vjust=-.3)

#flip nodes for arrangment
p3<- ggtree::flip(p3,189,216)  
p3<- ggtree::flip(p3,217,317)
p3<- ggtree::rotate(p3, 189) 

p3 <- p3 + geom_tiplab(aes(label=paste(x2@extraInfo$boldid, x2@extraInfo$country)), size=1.75)

#add annoated to clade colours

p6 <- p3 + annotate("text", -2.6, c(180), hjust=0, fontface=2,size=4, label=c(expression(bold("CO1 Lineage"))))+
  annotate("point",-2.5, 170, size=4, shape=21,fill="yellow",color="black") +
  annotate("point",-2.5, 160, size=4, shape=21,fill="green", color="black") +
  annotate("point",-2.5, 150, size=4, shape=21,fill="purple",color="black") +
  annotate("point",-2.5, 140, size=4, shape=21,fill="red",color="black") +
  annotate("text", -2.3,170, hjust=0, size=3.5, label="Iberia-Italy") + 
  annotate("text", -2.3,160, hjust=0, size=3.5,  label= "Palaeartic" )+
  annotate("text", -2.3,150, hjust=0, size=3.5, label="Alicante")+
  annotate("text", -2.3,140, hjust=0, size=3.5, label="Sierra Nevada")

#collapse clades
p6 <- collapse(p6, 195, 'max', fill="Purple", col="black", alpha=0.3) 
p6 <- collapse(p6, 219, 'max', fill="yellow", col="black", alpha=0.3) 
p6 <- collapse(p6, 282, 'max', fill="yellow", col="black", alpha=0.3)
p6 <- collapse(p6, 317, 'max', fill="green", col="black", alpha=0.3) 
p6 <- collapse(p6, 191, 'max', fill="Red", col="black", alpha=0.3)


#add clade labels with country locations
ibit= c(paste0("Italy"), paste0("Spain"),  paste0("Portugal"), paste0("S. France"))
ibitUK= c(paste0("Great Britain:"), paste0("Scotland (mainland)"),  paste0("Northern England"), paste0("Southern England"), paste0("Finland"))

palaeartic=c(paste0("France, Germany"), paste0("Romania, Finland"), paste0("Norway, Netherlands"), paste0("Iran, Turkey"), paste0("Austria, Canada"),
             paste0("Spain, Italy"), paste0("Greece"))
palaearticUK= c(paste0("Great Britain:"), paste0("Wales"), paste0("Southern England") )

apsnc= c(paste0("Great Britain:"), paste0("Outer Hebrides"), paste0("Scotland (mainland)"), paste0("Alicante, Germany"), paste0("Norway, Kazakhastan"))

p7 <- p6  +  annotate("text", x= c(.03,.03,.03,.03) ,y=c(78,75,72,69), hjust=0, size=2.5, label=ibit, Parse=TRUE) +
  annotate("text", x= c(.03,0.07,0.07,.07,.03) ,y=c(35,32,29,26,23), hjust=0, size=2.5, label=ibitUK, Parse=TRUE) +
  annotate("text", x= 0.03 ,y=c(150,147,144,141,138,135,132), hjust=0, size=2.5, label=palaeartic, Parse=TRUE) +
  annotate("text", x= c(.03,0.07,.07) ,y=c(120,117,114), hjust=0, size=2.5, label=palaearticUK, Parse=TRUE) +
  #annotate("text", x= 0.05 ,y=187, hjust=0, size=2.5, label="Crete", Parse=TRUE) +
  annotate("text", x= 0.03 ,y=184, hjust=0, size=2.5, label="Sierra Nevada", Parse=TRUE) +
  annotate("text", x= c(.03,.07,.07,.03, 0.03) ,y=c(178,175,172,169,166), hjust=0, size=2.5, label=apsnc, Parse=TRUE)



P8 <- p7 + ggtitle("A") + theme(plot.title = element_text(size= 20, hjust=0.05, vjust=-6, margin = margin(b = 0, t = -1, l = 2, unit = "cm")))

##-----------------------------------------------------------------
#supp figure with full tree

library(tidyverse)
library(treeio)
library(ggplot2)
library(ggtree)
library(grid)
library(scatterpie)
x <- read.beast("/media/data_disk/PROJECTS/Saad/CommonBlue/CO1_BOLD/BEAST/modTest_CCS_sc_split_subst/combined.mcc.tree")
setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/scripts/FigScripts")
#rename taxa only to have smaller names
tipnames <- x@phylo$tip.label
boldids <- sapply(strsplit(tipnames, "|", fixed=TRUE),  "[", 1)
country <- sapply(strsplit(tipnames, "|", fixed=TRUE),  "[", 2)
#replace errror and correct formatting
country<- str_replace_all(country, "Polyommatus_icarus", "Italy")
country<- str_replace_all(country, "United_Kingdom", "United Kingdom")
d <- data.frame(label = tipnames,
                boldid = boldids, 
                country= country, stringsAsFactors = F)
d$country[d$boldid=="EULEP004-14"] <- "Crete"
x2 <- full_join(x, d, by = "label")

#for removing bars for with support less than 0.8
x2@data$height_0.95_HPD[x2@data$posterior < 0.8] <- NA

p1 = ggtree(x2) + geom_treescale(x=-1.5, y=100, fontsize=4)  + scale_x_continuous(name="MYA", breaks=c(0,-0.5,-1,-1.5,-2,-2.5,-3),labels = abs, limits=c(NA,0.5)) + theme_tree2(plot.margin=margin(2, 2, 2, 2) )  
p1=revts(p1)
#add uncertainty of time
p2= p1+geom_range('height_0.95_HPD', color='black', size=2, alpha=.25)
#add posterior probaility as number
#p3= p2 +geom_nodelab(aes(x=branch, subset=!isTip & posterior > 0.7, label=round(posterior*100, 0)),hjust=1.5 size=2, nudge_x = 0.08) #+ geom_tiplab(aes(label=country),size=2)
p3 = p2 + geom_text2(aes(label=round(posterior, 1), subset=posterior>0.8), size=3, nudge_x = -0.05, vjust=-.3)

#flip nodes for arrangment
p3<- ggtree::flip(p3,189,216)  
p3<- ggtree::flip(p3,217,317)
p3<- ggtree::rotate(p3, 189) 

p3 <- p3 + geom_tiplab(aes(label=paste(x2@extraInfo$boldid, x2@extraInfo$country)), size=1.75)


#add annoated to clade colours

p6 <- p3 + annotate("text", -2.6, c(180), hjust=0, fontface=2,size=4, label=c(expression(bold("CO1 Lineage"))))+
  annotate("point",-2.5, 170, size=4, shape=21,fill="yellow",color="black") +
  annotate("point",-2.5, 160, size=4, shape=21,fill="green", color="black") +
  annotate("point",-2.5, 150, size=4, shape=21,fill="purple",color="black") +
  annotate("point",-2.5, 140, size=4, shape=21,fill="red",color="black") +
  annotate("text", -2.3,170, hjust=0, size=3.5, label="Iberia-Italy") + 
  annotate("text", -2.3,160, hjust=0, size=3.5,  label= "Palaeartic" )+
  annotate("text", -2.3,150, hjust=0, size=3.5, label="Alicante")+
  annotate("text", -2.3,140, hjust=0, size=3.5, label="Sierra Nevada")

P8 <- p6 + geom_hilight(218,fill = "yellow", alpha = 0.25) + geom_hilight(317,fill = "green", alpha = 0.25) + 
  geom_hilight(195,fill = "purple", alpha = 0.25) + geom_hilight(191,fill = "red", alpha = 0.25) 


pdf(file = "FigS2FullTree.pdf", width = 8, height = 12)
P8
dev.off()
