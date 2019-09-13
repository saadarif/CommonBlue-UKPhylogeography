#plotting infection rates on a map
library(maptools)
library(rgdal)
library(ggplot2)
library(scatterpie)
library(plyr)
library("ggmap")
library("ggrepel")

#read in infection data
infection <- read.delim("/media/data_disk/PROJECTS/Saad/CommonBlue/Wolbachia/Wolbachia_infection.csv", header=T)

#calculate average infection rates for each pop and sex
infection$sex <-  as.factor(sapply(strsplit(as.character(infection$X.sample), '_'), function(x){paste0(x[2])}))

ave_data <-aggregate(percentage_classified_wolcbachia ~ POP + sex, data = infection, FUN= "mean" )
ave_data$uninfected <- 100-ave_data$percentage_classified_wolcbachia

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
RNL=c(-0.281520, 53.254429) #should be Chamber's Farm wood CFW , lincolnshire 10
RHD=c(-1.477433, 54.713907) #Raisby hill durham 13
BMD=c(-1.125439, 51.782563) #Oxford 16


coords <- rbind(MLG, MDC, BER, TUL, DGC, RVS, OBN, BWD, FRN, PCP, ETB, MMS, RNL, RHD, BMD)
colnames(coords) <- c("lat", "lon")
#get in same order as the ave data
ord_coord <- as.data.frame(coords[order(rownames(coords)),])
#check
for (i in 1:dim(ave_data)[1])
{
  ave_data$lon[i] <- ord_coord$lat[rownames(ord_coord)==ave_data$POP[i]]
  ave_data$lat[i] <- ord_coord$lon[rownames(ord_coord)==ave_data$POP[i]]
}

names(ave_data)[3] <- "Infection"
#get sample sizes
sample_size <- count(infection, vars = c("POP", "sex"))
for (i in 1:dim(ave_data)[1])
{
  ave_data$sample[i] <- sample_size$freq[sample_size$POP==ave_data$POP[i] & sample_size$sex==ave_data$sex[i]]
}
#delete the french sample, no point showing on map
#ave_data_s<- subset(ave_data, POP!="FRN")

#get the map
UK <- map_data ("world", regions="UK")
ave_data_s_f <- subset(ave_data_s, sex=="f")
ave_data_s_m <- subset(ave_data_s, sex=="m")
#mapplot1 <- ggplot(UK) + 
#  geom_map(data = UK, map=UK,aes( map_id=region), col = "white", fill = "gray50") + theme_bw()

p <- ggplot(UK, aes(long, lat)) +
    geom_map(map=UK, aes(map_id=region), fill=NA, color="black") +
   coord_quickmap()

#mapplot1<- ggplot(data = UK) + 
#geom_polygon() #+
#coord_map()
p + geom_scatterpie(aes(x=lon, y=lat, group=sex, r=radius/30), 
                data = ave_data_s_f, cols =colnames(ave_data_s_f[,c(3:4)]),color=NA )  

#another way
map <- get_map(location = c(-2.65, 54.6), zoom = 6, maptype = "terrain-background", color="bw")
#location = c(lon = -2.5, lat =  54
p <- ggmap(map) +
geom_point(data = ave_data_s, mapping = aes(x = lon, y = lat,
                                      col = Infection, size = sample, shape=sex),  position=position_jitter(w = 0.45, h = 0.0)) + scale_colour_distiller(palette = "RdPu", direction = 1)  
p

#p + geom_point(data = ave_data_s_f, mapping = aes(x = lon, y = lat,
   #                                               col = Infection, size = sample), position=position_jitter(w = 0.5, h = 0.0)) + scale_colour_distiller(palette = "RdPu", direction = 1)
#geom_text(aes(label=ave_data_s_m$sex),hjust=0, vjust=0)


#one more try
p <- ggplot(UK, aes(long, lat)) +
  geom_map(map=UK, aes(map_id=region), col = "white", fill = "gray90")  + theme_bw() + coord_quickmap()
p + geom_point(data = ave_data_s, mapping = aes(x = lon, y = lat,
                                                col = Infection, size = sample, shape=sex),  position=position_jitter(w = 0.6, h = 0.0)) + scale_colour_distiller(palette = "YlOrRd", direction = 1)  
