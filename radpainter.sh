#!/bin/bash

#perform rad chromosome painting from a stacks VCF file
#requries instalattion of fineRADStructure, please change that path accordingly


#inputs
#$1 = input vcf file

#convert the snps file to haplotype file
/opt/fineRADstructure/RADpainter hapsFromVCF -H p $1 > haps_radpainter

#LD thining
Rscript /opt/fineRADstructure/sampleLD.R -s 1 -n 500 haps_radpainter haps_radpainter_reordered


#Calculate the co-ancestry matrix
/opt/fineRADstructure/RADpainter paint haps_radpainter_reordered
#Assign individuals to populations 
/opt/fineRADstructure/finestructure -x 100000 -y 100000 -z 1000 haps_radpainter_reordered_chunks.out haps_radpainter_reordered_chunks.mcmc.xml
#tree building
/opt/fineRADstructure/finestructure -m T -x 10000 haps_radpainter_reordered_chunks.out haps_radpainter_reordered_chunks.mcmc.xml haps_radpainter_reordered_chunks.mcmcTree.xml

#still need to plot results
