#missingness pca

#read in the missing data matrix from radpainter put

setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m5_M4_n4/populations.r50.p15_moh_0.65/RadPainter/miss75/")
missmat =read.delim("fineRADpainter.lociFilt.samples25%missFilt_reordered_missingnessMatrix.out")
misspca <- prcomp(missmat[,2:dim(missmat)[2]], center = T)
#get plate info
lane1 <- read.delim("/media/data_disk/PROJECTS/Saad/CommonBlue/INFO/barcodes_plate1", header=F)
lane2 <- read.delim("/media/data_disk/PROJECTS/Saad/CommonBlue/INFO/barcodes_plate2", header=F)

#combine
lane1$lane <- "1"
lane2$lane <- "2"
lanes <- rbind(lane1,lane2)
#extract only retained individuals
laneInfo <- lanes[lanes$V2 %in% rownames(misspca$rotation),2:3]
#add popinfo
laneInfo$pops <- as.factor(sapply(strsplit(as.character(laneInfo$V2), '_'), function(x){paste0(x[1])}))
#add sex
laneInfo$sex <- as.factor(sapply(strsplit(as.character(laneInfo$V2), '_'), function(x){paste0(x[2])}))
#colour the groups if necessayr
col <- funky(15)
plot(misspca$rotation[,1], misspca$rotation[,2], pch=c(16, 17)[as.numeric(laneInfo$lane)], col=transp(col,.8), 
     xlab=paste0("PC1 (", round(misspca$sdev[1]^2/sum(misspca$sdev^2) *100,2), "%)"), ylab=paste0("PC2 (",  round(misspca$sdev[2]^2/sum(misspca$sdev^2 )*100,2),"%)"), ylim=c(-0.4,0.4))
legend("bottomright", legend = c("Lane 1","Lane 2") , 
       pch=c(16, 17), bty="n" , cex = 1, horiz = FALSE, inset = c(0.03, 0.05))

#barplots of missing data
col <- funky(2)
mind <- read.delim("ineRADpainter.lociFilt.samples25%missFilt_reordered_missingness.out")
barplot(as.vector(t(mind)), ylim=c(0,0.5), ylab= "Percent Missig Data", col=laneInfo$lane)

#per pop
library(ggplot2)
stopifnot(names(mind) == laneInfo$V2)
laneInfo$miss <- as.vector(t(mind))
ggplot(laneInfo, aes(x = pops, y = miss, fill = lane)) + geom_boxplot() + theme_bw()
boxplot(as.vector(t(mind))~ laneInfo$pops , col=laneInfo$lane)
