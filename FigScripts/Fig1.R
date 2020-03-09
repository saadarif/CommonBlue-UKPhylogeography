#Fig1 Plot of sampling localaities in the UK
#2019/01/28 SA
library(tidyverse)
library(ggplot2)
require(maps)
library('ggsn') #for scale bar and symbol


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
#convert rownames to column
coords <- tibble::rownames_to_column(coords, "POP")

#read in data for samples
samples <- read.table("/media/data_disk/PROJECTS/Saad/CommonBlue/INFO/popmap_sex_copy.tsv", colClasses = c("NULL", "factor", "factor"))
samples <- as.data.frame(table(samples$V2))
colnames(samples) <- c("POP", "Samples")

#combine the data
sampling1 <- inner_join(coords, samples, by="POP")

#filter out data for the FRN site as it won't be plotted
sampling <- sampling1 %>% filter(POP != "FRN")

#plot on GB map with freq corresponding to size
range(sampling$lat)
range(sampling$lon)


#Get the region
UK <- map_data("world") %>% filter(region=="UK")

#plot the map
p1 <- ggplot() +
  geom_polygon(data = UK, aes(x=long, y = lat, group = group), fill="grey",alpha=0.3) +
  geom_point( data=sampling, aes(x=lon, y=lat, size=Samples), alpha=0.75) + geom_label(data=sampling, aes(x=lon, y=lat, label=POP), size=3, hjust = 0, nudge_x = 0.2)+
  theme_void() + ylim(50,59) + xlim(NA,2.2) + coord_map() + scale_size_continuous(name="Sample Size", breaks=c(6,10,16))+
  theme(legend.position = c(0.85, 0.8)) 

library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
require(rgeos)
library("ggspatial")

Fr <- ne_countries(scale = "medium", returnclass = "sf", continent="Europe")


p3 <- ggplot(data = Fr) +
  geom_sf(alpha=0.3) + theme_void() +  coord_sf( xlim=c(-5,10),ylim = c(41, 53), expand = TRUE) +
  geom_point( data=sampling1 %>% filter(POP=="FRN"), aes(x=lon, y=lat, size=Samples), alpha=0.6) + geom_label(data=sampling1 %>% filter(POP=="FRN"), aes(x=lon, y=lat, label=POP), size=3, hjust = 0, nudge_x = 0.2)+
  theme(legend.position ="none" )
p3

pdf(file = "Fig1.pdf", width = 8, height = 11)
p1
dev.off()



#-----------------------------------------------------------------------------------
#alternate maps with scale bar and north
p2 <- ggplot() +
  geom_polygon(data = UK, aes(x=long, y = lat, group = group), fill="grey",alpha=0.3) +
  geom_point( data=sampling, aes(x=lon, y=lat, size=Samples), alpha=0.75) + geom_label(data=sampling, aes(x=lon, y=lat, label=POP), size=3, hjust = 0, nudge_x = 0.2)+
  theme_void() + ylim(50,59) + xlim(NA,2.2) + coord_map() +scale_size_continuous(name="No. Specimens", breaks=c(6,10,16)) + 
  scalebar(UK, dist = 100, dist_unit = "km",transform = T, model = "WGS84", location = "bottomright",anchor=c(x=2.2, y=50.25), st.size=3, height=0.01, st.dist = 0.01)+
  theme(legend.position = c(0.85, 0.8),legend.title = element_text( size = 14),legend.text = element_text( size = 12) ) 
north2(p2, x = 0.74, y = 0.1, scale = 0.03, symbol = 12)


#plot.background = element_rect(fill = "lightskyblue", color = NA), 
#panel.background = element_rect(fill = "lightskyblue", color = NA), 
#legend.background = element_rect(fill = "lightskyblue", color = NA)

#alternative map with sf
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
require(rgeos)
library("ggspatial")

UK2 <- ne_countries(scale = "medium", returnclass = "sf", country="united kingdom")


p3 <- ggplot(data = UK2) +
  geom_sf(alpha=0.3) + theme_void() +  coord_sf( xlim=c(-8,3),ylim = c(50, 59), expand = TRUE) +
  geom_point( data=sampling, aes(x=lon, y=lat, size=Samples), alpha=0.6) + geom_label(data=sampling, aes(x=lon, y=lat, label=POP), size=3, hjust = 0, nudge_x = 0.2)+
  scale_size_continuous(name="No. Specimens", breaks=c(6,10,16)) + theme(legend.position = c(0.85, 0.8),legend.title = element_text( size = 14),legend.text = element_text( size = 12) )+
  annotation_scale(location = "br", width_hint = 0.3) +
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.79, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_orienteering, height = unit(.7, "cm"), width = unit(.7, "cm"))
p3

#Prepare supplemental table
tab <- read.table("/media/data_disk/PROJECTS/Saad/CommonBlue/INFO/popmap_sex_copy.tsv", colClasses = c("NULL", "factor", "factor"))
sexes <- as.data.frame(table(tab$V2, tab$V3))
#this gives data in long format covnert to wide
df <- spread(sexes, Var2, Freq)
names(df) <- c("POP", "f", "m")
df2 <- inner_join(sampling1, df, by="POP")
write.csv(df2, "TableS1", sep="\t")