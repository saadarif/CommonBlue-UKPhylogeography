#Making SFS from VCF

#USe the vcf2sfs.r

source("/media/data_disk/PROJECTS/Saad/CommonBlue/scripts/RScripts/vcf2sfs.r")

#generate a pop file with clusters from PC but remove all intermeidate pops
library(vcfR)
library(stringi)    
stri_replace_all_fixed("a'b c", pattern = c("'", " "), replacement = c("", "_"), vectorize_all = FALSE)
vcf <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/BayeScan/p15r50miss25/noFRNRVS/populations.snps.filter2.0.25.recode.singlesnp_neut.vcf")
z <- vcfR2genlight(vcf)
popnames <- z@ind.names
pops <- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[1])}))
pops<- str_replace_all(pops, "RNL", "CFW" )
#replace BER and TUL with 
pops<-stri_replace_all_fixed(pops, pattern ="BER", replacement = "OHS")
pops<- stri_replace_all_fixed(pops, "TUL", replacement="OHS")
#replace MLS, OBN and DGC with HGS
pops<-stri_replace_all_fixed(pops, "MLG", "HGS")
pops <-stri_replace_all_fixed(pops, "OBN", "HGS")
pops <- stri_replace_all_fixed(pops, "DGC", "HGS")
#replace BMD, ETB, BWD, MMS, PCP, CFW, MDC as SBI
pops<- str_replace_all(pops, "BMD", "SBI")
pops<- str_replace_all(pops, "ETB", "SBI")
pops<- str_replace_all(pops, "BWD", "SBI")
pops<- str_replace_all(pops, "MMS", "SBI")
pops<- str_replace_all(pops, "PCP", "SBI")
pops<- str_replace_all(pops, "CFW", "SBI")
pops<- str_replace_all(pops, "MDC", "SBI")

#combine and write to file
popmap <- data.frame(ind=z@ind.names, pop=pops)
write.table( popmap, "/media/data_disk/PROJECTS/Saad/CommonBlue/scripts/RScripts/popmapSFS.txt", sep="\t", col.names=FALSE, quote = F, row.names = F)


#read in vcf and popmap file for making the SFS
mygt<-vcf2gt("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/BayeScan/p15r50miss25/noFRNRVS/populations.snps.filter2.0.25.recode_neut.vcf",
            "/media/data_disk/PROJECTS/Saad/CommonBlue/scripts/RScripts/popmapSFS.txt")
mysfs1<-gt2sfs.raw(mygt, "SBI")
plot.sfs(fold.sfs(mysfs1))
mysfs2 <- gt2sfs.raw(mygt, "HGS")
plot.sfs(fold.sfs(mysfs2))
mysfs3 <- gt2sfs.raw(mygt, "OHS")
plot.sfs(fold.sfs(mysfs3))
#2dsfs
mysfs4<-gt2sfs.raw(mygt, c("OHS", "HGS"))
plot.sfs(mysfs4)
mysfs5<-gt2sfs.raw(mygt, c("HGS", "SBI"))
plot.sfs(mysfs5)


