#Figures for determining suitable parameters
#from stacks testing
# Feb 2020 saad

library("ggplot2")
library(grid)
library(gridExtra)
library(ggpubr)
library(tidyverse)
#directory for saving
setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/scripts/FigScripts/")


#There is a folder for stacks 2-12 in each folder there is a folder
#called "results" with a file called n_snps_per_locus.tsv which has the relevant data
DIR<- dir("/media/data_disk/PROJECTS/Saad/CommonBlue/tests.denovo", pattern="stacks*", full.names = T)
#append results to end of each file
res_dir <- paste0(DIR, "/results")
#get the full paths to each file
out_files <- sapply(res_dir, dir, pattern="\\.tsv$", full.names=TRUE)
#read all files into a list


all <- lapply(out_files, read.csv, sep="\t")

# the final data should have a column for m (M=n=2) a column for all loci and a column for only polymorphic loci
m_analysis=data.frame(m=numeric(), status=character(), n=numeric())
#populate the data frame
for (i in all){
  Mn_2 = subset(i, M==2 & n==2)
  df_r1 <- data.frame(m=Mn_2$m[1], status="total", n=sum(Mn_2$n_loci[1:length(Mn_2$n_loci)]))
  df_r2 <- data.frame(m=Mn_2$m[1], status="polymorphic", n=sum(Mn_2$n_loci[2:length(Mn_2$n_loci)]))
  m_analysis <- rbind(m_analysis, df_r1, df_r2)
}

p1 <- ggplot(m_analysis, aes(x=m, y=n, group=status, color=status)) + geom_point(size=2)+ theme_bw() +
  xlab("m") + ylab("No. of loci") + scale_colour_discrete(name="",labels=c("Putative Rad Loci", "Polymorphic Rad Loci"))+
  scale_x_continuous(breaks=1:12)+
  theme(legend.position=c(0.75,0.8),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(1,5,1,1), "lines"),
        legend.title = element_text( size = 14),
        legend.text = element_text( size = 12),
        legend.key.size = unit(3,"line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )

#same for M=n
m3 <- read.csv("/media/data_disk/PROJECTS/Saad/CommonBlue/tests.denovo/stacks.m3/results/n_snps_per_locus.tsv", sep="\t")
Mn_analysis=data.frame(Mn=numeric(), status=character(), n=numeric())
#populate the data frame
for (i in 1:8) {
  m3_M <- subset(m3, M==i)
  df_r1 <- data.frame(Mn=m3_M$M[1], status="total", n=sum(m3_M$n_loci[1:length(m3_M$n_loci)]))
  df_r2 <- data.frame(Mn=m3_M$M[1], status="polymorphic", n=sum(m3_M$n_loci[2:length(m3_M$n_loci)]))
  Mn_analysis <- rbind(Mn_analysis, df_r1, df_r2)
}

p2 <- ggplot(Mn_analysis, aes(x=Mn, y=n, group=status, color=status)) + geom_point()+ theme_bw() +
  xlab("M=n") + ylab("No. of loci") + scale_colour_discrete(name="",labels=c("Putative Rad Loci", "Polymorphic Rad Loci"))+
  scale_x_continuous(breaks=1:8)+
  theme(legend.position=c(0.75,0.4),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(1,5,1,1), "lines"),
        legend.title = element_text( size = 14),
        legend.text = element_text( size = 12),
        legend.key.size = unit(3,"line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )

pdf(file = "FigS1AB.pdf", width = 8, height = 14)
ggarrange(p1,p2, labels = c("A", "B"), ncol=1, font.label = list(size = 20, color = "black", face = "bold", family = NULL))
dev.off()

d <- read.csv("/media/data_disk/PROJECTS/Saad/CommonBlue/tests.denovo/stacks.m4/results/n_snps_per_locus.tsv", sep="\t")
m_val <- 4

d = subset(d, M==n & m==m_val)
# Make sure the table is ordered by number of snps.
d = d[order(d$n_snps),]

Mn_values = sort(unique(d$M))

# Write the counts in a matrix.
m = matrix(NA, nrow=length(Mn_values), ncol=max(d$n_snps)+1)
for(i in 1:nrow(d)) {
  m[d$M[i],d$n_snps[i]+1] = d$n_loci[i] # [n_snps+1] as column 1 is for loci with 0 SNPs
}

# Truncate the distributions.
max_n_snps = 10
m[,max_n_snps+2] = rowSums(m[,(max_n_snps+2):ncol(m)], na.rm=T)
m = m[,1:(max_n_snps+2)]
m = m / rowSums(m, na.rm=T)

# Draw the barplot.
pdf('./FigS2.pdf')

col = rev(heat.colors(length(Mn_values)))

x <- barplot(m,
        beside=T, col=col, las=1,
        names.arg=c(0:max_n_snps, paste('>', max_n_snps, sep='')),
        xlab='Number of SNPs',
        ylab='Percentage of loci',
        ylim=c(0,.2),
        main='Distributions of the number of SNPs per locus\nfor a range of M=n values at m=4'
)
legend(90, 0.17, legend=c( Mn_values), fill=c(NA, col))
text(93, 0.18, labels="M=n", cex=1.5)

null=dev.off()

