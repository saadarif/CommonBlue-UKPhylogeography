#!/usr/bin/env Rscript

#read in files from the terminal, full paths seperates by space
args = commandArgs(trailingOnly=T)

# test that both arguements are present
if (length(args)!=2) {
  stop("Please provide full paths to both kmer frequency files (seperated by a space)", call.=FALSE)
} 

# turn off warnings for quiet output
options(warn=-1)

#use ggplot2
library("ggplot2")

#read in the files
plate1=args[1]
plate2=args[2]

#read data with hardcodes paths, not used at the moment
plate1_kmers = read.table(plate1,header=T)
plate2_kmers = read.table(plate2,header=T)

#combine the df's
plate1_kmers$group <- "plate1"
plate2_kmers$group <- "plate2"
both<-rbind(plate1_kmers,plate2_kmers)
both$group <- as.factor(both$group)


#make the plots
p <- ggplot(both, aes(x = KmerFrequency, y = Count, group=group))+ geom_line(aes(colour=group))+scale_x_continuous(limits = c(0, 4000), name="K-mer Frequency")+scale_y_continuous(trans = 'log10', limits=c(100,10000000), name=expression("Log"[10]*" K-mer Count"))
p <- p + theme_bw()
p <-p + theme(text = element_text(size=20),legend.title = element_blank(),legend.position = c(.8,.75),legend.text=element_text(size=20, face = "bold"))
p <- p + coord_fixed(ratio = 0.25)
#initialize and save plots as pdf
pdf("/media/data_disk/PROJECTS/Saad/CommonBlue/plots/kmer_dists_plots.pdf")
print(p)
dev.off()
