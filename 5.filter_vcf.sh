#!/bin/bash

#script that further filters vcf output from stacks populations
#for sites with too many snps (>5), individuals with too much missing data (>20%)
#and any additional contaminants found from centrifuge analysis or elsewhere

#get base level directory
top=$(readlink -f $(dirname $0)/..)

#directory where vcf output is
vcf_dir=$top/stacks.denovo/stacks_m5_M4_n4_new/populations.r50.p15_moh_0.65

#remove loci with more x amount of snps
bad=5
cut -f 1  $vcf_dir/populations.snps.vcf | sort | uniq -c| sort |\
 mawk -v a="$bad" '{if($1 > a) print $2}' > $vcf_dir/pos.tmp

while read line; do 
	grep -w "^$line" $vcf_dir/populations.snps.vcf | cut -f1,2 
 done < $vcf_dir/pos.tmp > $vcf_dir/badloci

#filter using vcf tools
vcftools --vcf $vcf_dir/populations.snps.vcf --exclude-positions $vcf_dir/badloci --recode \
--recode-INFO-all --out $vcf_dir/populations.snps.filter1

#next filter based on missingness in idividuals
#we use indiviuals missing more than 25% data or missing more than 50%
miss="0.25 0.5"
vcftools --vcf $vcf_dir/populations.snps.filter1.recode.vcf --missing-indv --out $vcf_dir/out
#this shoud output a file called  out.imiss

for b in $miss; do
	mawk -v a="$b" '{if($5 > a) print $1}' $vcf_dir/out.imiss  > $vcf_dir/lowDP.indv
        vcftools --vcf $vcf_dir/populations.snps.filter1.recode.vcf --remove $vcf_dir/lowDP.indv \
        --recode --recode-INFO-all --out $vcf_dir/populations.snps.filter2.${b}
done
