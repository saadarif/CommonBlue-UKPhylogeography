#!/bin/bash

#script that further filters vcf from step 5 to only include unlinked loci 
#(first snp of each loci) for fastStructure analysis

#get base level directory
#top=$(readlink -f $(dirname $0)/..)


#start with the vcf file filtered with missingness and bad loci
#workking wth p15 (all pops) and r50
#vcf_dir=$top/stacks.denovo/stacks_m4_M4_n4/populations.r50.p15_moh_0.65/
#in_dir=$top/stacks.denovo/stacks_m4_M4_n4/populations.r50.p15_moh_0.65/

#cd $in_dir
#mkdir -p fastStructure

#the input vcf file is also filtered for loci with > 5 snps
#and individuals with more than 25% missing data remove
#start with the previously filtered vcf and keep only the first snp from each loci
#awk '!seen[$1]++' populations.snps.filter2.0.25.recode.vcf > fastStructure/populations.snps.filter2.0.25.recode.singlesnp.vcf 


#conver the vcf to faststructrue using PGD spider and previously created spid file
#java -Xmx1024m -Xms512m -jar /opt/PGDSpider2/PGDSpider2-cli.jar \
#-inputfile fastStructure/populations.snps.filter2.0.25.recode.singlesnp.vcf -inputformat vcf \
#-outputfile fastStructure/populations.snps.filter2.0.25.recode.p.snps.fast.str -outputformat STRUCTURE \
#-spid stacks_vcf_to_fastStr.spid 

#K_values="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15"

#mkdir -p  fastStructure/fastStructureOut

#irt run faststructure with a simple prior then try logistic prior if require to look for subtle structure. 
#for K in $K_values; do
#    echo "now executing fastStructre for K = $K ..."
#    python /opt/fastStructure/structure.py -K $K \
#    --input=fastStructure/populations.snps.filter2.0.25.recode.p.snps.fast \
#    --output=fastStructure/fastStructureOut/simple  --full --prior=simple --format=str
#done

#choosing model complexity
#echo  "now running model complexity script..."
#python /opt/fastStructure/chooseK.py --input=fastStructure/fastStructureOut/simple \
#      > fastStructure/fastStructureOut/ model_test_simple.out

#bssed on results of the simple prior run with logsitic prior for subtle strucute or k=2-4
in_dir=/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4/populations.r50.p15_moh_0.65/
cd $in_dir
for K in {2..4}; do
    echo "now executing fastStructre for K = $K ..."
    python /opt/fastStructure/structure.py -K $K \
    --input=fastStructure/populations.snps.filter2.0.25.recode.p.snps.fast \
    --output=fastStructure/fastStructureOut/logistic  --full --prior=logistic --format=str
done

#choosing model complexity
echo  "now running model complexity script..."
python /opt/fastStructure/chooseK.py --input=fastStructure/fastStructureOut/logistic \
      > fastStructure/fastStructureOut/model_test_logistic.out
