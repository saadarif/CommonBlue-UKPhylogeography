#prepare wolbachia rad loci for making a time tree

system("ln -s ../populations.snps.filter2.0.25.recode.vcf .")
setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/Wolbachia/stacks_m4_M4_n4/populations.r50.p3/")
#convert vcf to phylip
system("python /opt/vcf2phylip/vcf2phylip.py -i  populations.snps.filter2.0.25.recode.vcf -n ")
#move files to current director


#filter invariant sites and save as phylip for raxml
library("phrynomics")
snps <- ReadSNP("populations.snps.filter2.0.25.recode.min4.nexus")
# Remove Invariant/nonbinary sites or loci/sites with too much missing data
snps_i <- RemoveInvariantSites(snps)
#write out phylipc for RAXML
WriteSNP(snps_i, file="populations.snps.filter2.0.25.recode.min4_noIinv.phy", format="phylip")