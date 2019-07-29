#!/usr/bin/env Rscript

#convert a filtered vcf to pylip and run raxml for groups
setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4/populations.r50.p15_moh_0.65/TreeBuilding/")
#convert vcf to phylip
system("python /opt/vcf2phylip/vcf2phylip.py ../populations.snps.filter2.0.25.recode.vcf -n -b -f")
#move files to current director
system("mv ../*min4* .")

#filter invariant sites and save as phylip for raxml
library("phrynomics")
snps <- ReadSNP("populations.snps.filter2.0.25.recode.min4.nexus")
# Remove Invariant/nonbinary sites or loci/sites with too much missing data
snps_i <- RemoveInvariantSites(snps)
#write out phylipc for RAXML
WriteSNP(snps_i, file="populations.snps.filter2.0.25.recode.min4_noIsnps.phy", format="phylip")

#run raxml
system("mkdir -p raxml_out")
system("raxmlHPC-PTHREADS-SSE3 -f a -T 2 -m ASC_GTRCAT --asc-corr=lewis -N 1000 -x 12345 -p 54321 -n m4r50n5miss25icarus-tree \
       -w /media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4/populations.r50.p15_moh_0.65/TreeBuilding/raxml_out \
       -s /media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4/populations.r50.p15_moh_0.65/TreeBuilding/populations.snps.filter2.0.25.recode.min4_noIsnps.phy")