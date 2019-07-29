#!/usr/bin/env Rscript

#Plots PCA plot of various stacks/populations threshholds
library("adegenet")
library(vcfR)

top="/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/"
setwd(top)
stacks_dir <- dir()
folder_temp="populations.r50.p15_moh_0.65"
#different r cuttoffs
r_vals=c("50", "60","70","75")
vcf_temp="populations.snps.filter2.0.5.recode.vcf"
#different missigness cutoffs
missing_vals=c("0.5", "0.25")

#file for plotting the pcas
pdf('../plots/AllPCAs.pdf')

for (s_dir in stacks_dir){
  m_dir=paste0(top,"/", s_dir)
  #get stacks paramter m value
  m = unlist(strsplit(s_dir, '_'))[2]
  for (r in r_vals){
    folder <- gsub("50", r , folder_temp)
    vcf_dir=paste0(m_dir, "/",folder)
    #read in vcf files with 50% and 25% filters
    for (miss in missing_vals) {
    vcf <- gsub("0.5", miss , vcf_temp)
    vcf_file <- paste0(vcf_dir,"/", vcf)
    print(vcf_file)
    #plot the pcas
    VCF <- read.vcfR(vcf_file)
    z <- vcfR2genlight(VCF)
    #get pops
    popnames <- z@ind.names 
    pops <- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[1])}))
    pop(z) = pops
    #impute missigng values
    w <- tab(z, freq = TRUE, NA.method = "mean")
    
    #perform PCA
    w.pca <- dudi.pca(w, scannf = F, scale=F, nf=3)
    
    #PCA plot
    col <- funky(15)
    s.class(w.pca$li, pop(z) ,xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
            cstar=0, cpoint=2, clabel=.75, grid=FALSE, sub=paste0(m," ", "r:", r, " ", "missing: ", miss))
    #add scree to plot
    #add.scatter.eig(w.pca$eig[1:10], xax=1, yax=2)
    }
    
  }
}

null=dev.off()