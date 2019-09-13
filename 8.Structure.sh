#!/bin/bash

#Structure threader script

#get base level directory
top=$(readlink -f $(dirname $0)/..)

#the structure input is composed of the same vcf file use for faststructure (scritp 6.fastStructure
#i.e. loci with > 5 snps filtered, only first snp from other loci and indivisuals with missingess > 25% removed

#change to appropriate dir

w_d=$top/stacks.denovo/stacks_m4_M4_n4/populations.r50.p15_moh_0.65/Structure
cd $w_d
#convert the vcf to structrue using PGD spider and previously created spid file
#java -Xmx1024m -Xms512m -jar /opt/PGDSpider2/PGDSpider2-cli.jar \
#-inputfile ../fastStructure/populations.snps.filter2.0.25.recode.singlesnp.vcf -inputformat vcf \
#-outputfile populations.snps.filter2.0.25.recode.singlesnp.str -outputformat STRUCTURE \
#-spid ../stacks_vcf_to_Str.spid

#manually curate the mainparams and extraparams file in the same dir
#manually curate indfile for plotting

#run structure threader for k=6 (amin clusters in PCA) with 5 replicates each
mkdir -p StructureResults_K1_8
structure_threader run -K 8 -R 5 -i populations.snps.filter2.0.25.recode.singlesnp.str -o StructureResults_K1_8 \
                        --ind indfile -t 4 -st /usr/local/share/Structure2.3.4_console/structure

