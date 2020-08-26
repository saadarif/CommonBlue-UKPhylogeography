#!/usr/bin/env Rscript

#convert a filtered vcf to pylip and run raxml for groups, use cvf with single snps
system("ln -s ../populations.snps.filter2.0.25.recode.vcf .")
setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/raxml/p15r50miss25neutral/")
#convert vcf to phylip
system("python /opt/vcf2phylip/vcf2phylip.py -i  populations.snps.filter2.0.25.recode.vcf -n ")
#move files to current director


#filter invariant sites and save as phylip for raxml
library("phrynomics")
snps <- ReadSNP("populations.snps.filter2.0.25.recode.neutral.min4.nexus")
# Remove Invariant/nonbinary sites or loci/sites with too much missing data
snps_i <- RemoveInvariantSites(snps)
snps_b <- RemoveNonBinary(snps_i)
#write out phylipc for RAXML
WriteSNP(snps_b, file="populations.snps.filter2.0.25.recode.min4_BInoIinv.phy", format="phylip")

#run raxml
system("mkdir -p raxml_out")
system("raxmlHPC-PTHREADS-AVX -f a -T 10 -m ASC_GTRCAT -V --asc-corr=lewis -N autoMRE -x 12345 -p 54321 -n m4r50n5miss25icarus-tree \
       -w /Storage/PROJECTS/Saad/PS001_CommonBlue/TreeBuilding/raxml_out \
       -s /Storage/PROJECTS/Saad/PS001_CommonBlue/TreeBuilding/populations.snps.filter2.0.25.recode.singlesnp_min4_noIinv.phy")

#1using raxml-ng outgroup is for drawing, following tuorial from https://github.com/amkozlov/raxml-ng/wiki/Tutorial#tree-inference
#--outgroup FRN_m_M05,FRN_f_M02,FRN_m_M03,FRN_f_M04,FRN_m_M06,FRN_m_M07
system("raxml-ng --parse --msa populations.snps.filter2.0.25.recode.singlesnp_min4_noIinv.phy --model GTR+ASC_LEWIS --prefix T2")
#2infer trees --outgroup FRN_m_M05,FRN_f_M02,FRN_m_M03,FRN_f_M04,FRN_m_M06,FRN_m_M07
system("raxml-ng --msa populations.snps.filter2.0.25.recode.neutral.min4.phy --model GTR+G+ASC_LEWIS  --prefix withInv --threads 8 --seed 2 --tree pars{25},rand{25}")
#3check differences in topology
system("raxml-ng --rfdist --tree T4.raxml.mlTrees --prefix RF")
#4bootstrapping
system("raxml-ng --bootstrap --msa populations.snps.filter2.0.25.recode.singlesnp_min4_noIinv.phy --model GTR+ASC_LEWIS --prefix T7 --seed 2 --threads 2")
#4a test bootstrap convergence at a more stringent level
system("raxml-ng --bsconverge --bs-trees T7.raxml.bootstraps --prefix T9 --seed 2 --threads 2 --bs-cutoff 0.01")
#5 left at defualt bs cuttoff (3%) generate bs support stree
system("raxml-ng --support --tree T4.raxml.bestTree --bs-trees T7.raxml.bootstraps --prefix T13 --threads 2")
#generate support trees with tbe
system("raxml-ng --support --tree T4.raxml.bestTree --bs-trees T7.raxml.bootstraps --prefix T14_tbeSup --threads 2 --bs-metric tbe")
#generatre majortiy rule consense trees from bootstreps
system("raxml-ng --consense -tree T7.raxml.bootstraps --prefix T15 --threads 2")

system("raxml-ng --all -msa ../populations.snps.filter2.0.25.recode.singlesnp_min4_noIinv.phy --model GTR+ASC_LEWIS \
        --tree pars{10} --outgroup FRN_m_M05,FRN_f_M02,FRN_m_M03,FRN_f_M04,FRN_m_M06,FRN_m_M07 --bsconverge ")


#building a tree using concataneated variants used for Radpainter
#50% missing data (more inds)
system("python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < /media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4/populations.r50.p15_moh_0.65/RadPainter/maxSnps5/miss50/haps_radpainter.fineRADpainter.lociFilt.samples50%missFilt.txt > /media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4/populations.r50.p15_moh_0.65/TreeBuilding/radpainter_50missfilt5snps.txt")
