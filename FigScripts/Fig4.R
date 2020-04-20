#fig 4 ect. figures to wolbachia infection and distribution
#Febd 2020 SA

library(ggplot2)
library(tidyverse)
require(maps)
library(grid)
library(gridExtra)
library(ggpubr)

setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/scripts/FigScripts")

wolbachia <- read.csv("../../Wolbachia/Wolbachia_infection.csv", sep="\t")
ggplot(wolbachia, aes(x=log(percentage_classified_wolcbachia))) + geom_histogram()
wolbachia[log(wolbachia$percentage_classified_wolcbachia)>1,]

#create a binary variable for whether individual was infected or not based on threshold of log(percent)>1

wol <- wolbachia %>% mutate(infected=ifelse(log(percentage_classified_wolcbachia)>1,1,0))

#create column for male/females
wol <- wol %>% mutate(sex=str_split(X.sample,"_", simplify = T)[,2])

#replace RNL with CFW
wol$POP<- str_replace_all(wol$POP, "RNL", "CFW")


#make two plots 
#A: map with infection levels
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
coords <- as.data.frame(coords)
coords$POP <- rownames(coords)

#pairwise fishers tests for proportions
wol1<- subset(wol,  POP=="BER" | POP=="TUL" |POP=="DGC" |POP=="MLG" |POP=="OBN" |POP=="RHD" )
pairwise.fisher.test(wol1$infected, wol1$POP,p.adjust.method = "bon")

#suummarize infection levels
infection <- as.data.frame.matrix(table(wol$POP, wol$infected))
names(infection) <- c("un", "inf")
infection$prop <- infection$inf/(infection$un+infection$inf)
infection$POP <- rownames(infection)
#merge the lat long info for plotting
infection <- full_join(infection, coords, by = "POP")
#calculate standard errors for the proportion
infection <- infection %>% mutate(se=sqrt((prop*(1-prop))/(inf+un)))
#note france won't be drawn
UK <- map_data("world") %>% filter(region=="UK")

#plot the map
p1 <- ggplot() +
  geom_polygon(data = UK, aes(x=long, y = lat, group = group), fill="grey",alpha=0.3) +
  geom_point( data=infection%>%filter(prop>0), aes(x=lon, y=lat, size=prop),col="red", alpha=0.7) + geom_label(data=infection, aes(x=lon, y=lat, label=POP), size=3, hjust = 0, nudge_x = 0.35)+
  theme_void() + ylim(50,59) + xlim(NA,2.2) + coord_map() +scale_size_continuous(name="Proportion of Individuals Infected", breaks=c(0.25,.50,.75), range=c(5,10)) + 
  theme(legend.position = c(0.85, 0.8))
p1


#barplot with map inset
#relevel the order for POP
infection$POP <- factor(infection$POP, c("RHD", "RVS", "OBN", "MLG", "DGC", "BER", "TUL"))
ibp<-ggplot(infection %>% filter(prop>0), aes(x=POP, y=prop)) + geom_bar(stat="identity", alpha=0.5) + theme_bw() +
ylim(c(0,1)) + xlab("Population") + ylab("Proportion of Individuals Infected")+ coord_flip() +
geom_errorbar(aes(x=POP, ymin=prop-se, ymax=prop+se), width=0.1, colour="black", alpha=0.5, size=0.75)+
theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"), 
      plot.margin = unit(c(1,1,1,1), "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())



#plot the map
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
require(rgeos)
library("ggspatial")

UK2 <- ne_countries(scale = "medium", returnclass = "sf", country="united kingdom")


p3 <- ggplot(data = UK2) +
  geom_sf(alpha=0.3) + theme_void() +  coord_sf( xlim=c(-8,1),ylim = c(54, 59), expand = TRUE) +
  geom_label(data=infection, aes(x=lon, y=lat, label=POP), size=2, hjust = 0) +
  theme(plot.background = element_rect( colour = "black", size=1))

library(grid)
ibpm <- ibp+
  annotation_custom(ggplotGrob(p3), xmin = -1.5, xmax = 4, 
                    ymin = 0.67, ymax = 1)
# ggtitle("A") plot.title = element_text(size= 20, hjust=0.05, vjust=-7, margin = margin(b = -0.1, t = -1, l = 2, unit = "cm"))
#boxplot for infection levels of norther populations

north <- wol %>% filter(POP=="BER" | POP=="TUL" | POP== "DGC" | POP=="MLG" | POP=="OBN")


#summarize the data
summ_north <- north%>% group_by(POP, sex) %>% 
  mutate(mean=mean(percentage_classified_wolcbachia), stde=sd(percentage_classified_wolcbachia)/sqrt(n())) %>%
  distinct(POP,.keep_all=TRUE)

#relevel the order for POP
north$POP <- factor(north$POP, c("OBN", "MLG", "DGC", "BER", "TUL"))


p2 <- ggplot(north, aes(y=percentage_classified_wolcbachia, x=POP, group=sex, colour=sex)) +
      geom_jitter(size=3, alpha=0.5, width = 0.25) +  coord_flip() + theme_bw() + scale_colour_manual(values=c("black", "grey80"))+
      xlab("Population") + ylab("Percentage of Classified Reads Mapped to Wolbachia") + geom_point(data=summ_north, aes(y=mean, x=POP, fill=sex),colour="black",size=6, pch=21, alpha=0.5,position=position_dodge(width=0.5))+
      #geom_errorbar(data=summ_north,aes(ymin = mean - stde, ymax = mean + stde),width = 0.25, size = 0.5,position=position_dodge(width=0.5) )+
      guides(colour=FALSE) + scale_y_continuous(breaks=c(0,10,20,30,40,50))+scale_fill_manual(values=c("black", "grey80"))+
      theme(legend.position=c(0.8,0.2),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"), 
            plot.margin = unit(c(1,1,1,1), "lines"),
            legend.title = element_text( size = 14),
            legend.text = element_text( size = 12),
            legend.key.size = unit(3,"line"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())

  
p2

#+ ggtitle("B") plot.title = element_text(size= 20, hjust=-0.15, vjust=7, margin = margin( l = 2, unit = "cm"))

pdf(file = "Fig4.pdf", width = 14, height = 8)
ggarrange(ibpm,p2, labels = c("A", "B"), font.label = list(size = 20, color = "black", face = "bold", family = NULL))
dev.off()


#testing custome annotations
#solution from: https://stackoverflow.com/questions/10525957/how-to-draw-lines-outside-of-plot-area-in-ggplot2

Text1 = textGrob("Outer\n Hebrides")

p3 <- ggplot(north, aes(y=percentage_classified_wolcbachia, x=POP, colour=sex)) +
  geom_jitter(size=2, alpha=0.5, width = 0.25) +  coord_flip() + theme_bw() +
  xlab("Population") + ylab("Pecentages of classified reads mapped to Wolbachia") + geom_point(data=summ_north, aes(y=mean, x=POP, fill=sex),colour="black",size=5, pch=21, alpha=0.5,position=position_dodge(width=0.5))+
  geom_errorbar(data=summ_north,aes(ymin = mean - stde, ymax = mean + stde),width = 0.25, size = 0.5,position=position_dodge(width=0.5) )+
  guides(colour=FALSE) + scale_y_continuous(breaks=c(0,10,20,30,40,50))+
  theme(legend.position=c(0.8,0.2),axis.text=element_text(size=14, face="bold"),axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(1,5,1,1), "lines"),
        legend.title = element_text( size = 14),
        legend.text = element_text( size = 12),
        legend.key.size = unit(3,"line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


p4 <- p3 +  annotation_custom(Text1,  xmin = 4, xmax = 5, ymin = 63, ymax = 64) +
      annotation_custom(grob = linesGrob(), xmin = 3.75, xmax =5.25, ymin = 58, ymax = 58)

#gets rid of the clipping
gt <- ggplot_gtable(ggplot_build(p4))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
x <- grid.draw(gt)

dev.off()


#supplemental figures histogramm abova
p5 <- ggplot(wolbachia, aes(x=log(percentage_classified_wolcbachia))) + geom_density(fill="grey80", alpha=0.5) + theme_bw() +
  xlab(expression("Log"[2] *"(percentage of classified reads mapped Wolbachia)")) +geom_vline(xintercept = 0, linetype="dotted") +
  theme(axis.text=element_text(size=14, face="bold"),axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(1,5,1,1), "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )

p6 <- ggplot(wolbachia, aes(x=log10(percentage_classified_wolcbachia), y=log10(n_reads))) + geom_point() + theme_bw()+
  xlab(expression("Log"[10] *"(percentage of classified reads mapped to Wolbachia)")) + 
  ylab(expression("Log"[10] *"(Total Reads)")) + 
  theme(axis.text=element_text(size=14, face="bold"),axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(1,5,1,1), "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )

cor.test(log10(wolbachia$percentage_classified_wolcbachia), log10(wolbachia$n_reads), method="spearman")

pdf(file = "FigS6AB.pdf", width = 8, height = 14)
ggarrange(p6,p5, labels = c("A", "B"), ncol=1, font.label = list(size = 20, color = "black", face = "bold", family = NULL))
dev.off()
wolbachia[log(wolbachia$percentage_classified_wolcbachia)>0,]

#----------------------------------------------------------------------------------------------------
#supplemental table stats
wol1<- subset(wol,  POP=="BER" | POP=="TUL" |POP=="DGC" |POP=="MLG" |POP=="OBN" |POP=="RHD" )
library(reporttools)
pairwise.fisher.test(wol1$infected, wol1$POP,p.adjust.method = "bon")

#-----------------------------------------------------------------------------------------
#correlation of percent wolbachia with average covrage of stack
wolradstats <- read.csv("/media/data_disk/PROJECTS/Saad/CommonBlue/scripts/FigScripts/TableS6_RADstats_Wolbachia_short.csv")
names(wolbachia)[1] <- "Sample"
corr <- inner_join(wolradstats, wolbachia,"Sample")
ggplot(corr, aes(x=log(Average.Coverage), y=log(percentage_classified_wolcbachia))) + geom_point()+ theme_bw()
cor.test(corr$Average.Coverage, corr$percentage_classified_wolcbachia, method="spearman")
