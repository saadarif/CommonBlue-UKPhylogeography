#!/bin/bash
top=$(readlink -f $(dirname $0)/..)
echo $top
cd $top/tests.denovo/control

#iterate over the m_values
m_values="3 4 5 6 7 8 9 10" 

mkdir -p results
cd results

echo "tallying the numbers"
echo -e '#par_set\tM\tn\tm\tshared_loci\tremoved_loci\tallele_mis\tn_snps' > CON_error_rates.tsv
#allele error rate = allele mismatch/shared loci
#snp error rate = n_snps/shared_loci
for m in $m_values ;do
	M_values="1 2 3 4 5 6 7 8"
	for M in $M_values ;do
	n=$M
	#for shared_loci and removed loci (sum equals total loci) and n_variant sites
	log_file=../stacks.m${m}/stacks.M${M}/populations.p2/populations.log
	#for allele mistmatch
	haplotype_file=../stacks.m${m}/stacks.M${M}/populations.p2/populations.haplotypes.tsv
	#vcf file for snp_mismatches
	vcf_file=../stacks.m${m}/stacks.M${M}/populations.p2/populations.snps.vcf
	
	#get shared and removed loci
	n_shared=$(cat $log_file | grep "Kept" | awk '{print $2}')
	n_removed=$(cat $log_file | grep "Removed" | awk '{print $2}')
	
	#get allele mistmatch from haplotype file
	allele_mismatches=$(awk '$3!=$4  {print $0}' $haplotype_file | wc -l)
	#remove 1 count for the header
	let allele_mismatches=allele_mismatches-1
	#ignore header files in vcf and extract and save genotypes for each rep
	tail -n +16 $vcf_file | cut -f10 | cut -d : -f1 > geno1.tmp
	tail -n +16 $vcf_file | cut -f11 | cut -d : -f1 > geno2.tmp
	n_snps=$(paste geno1.tmp geno2.tmp | awk '$1!=$2 {print$0}' |wc -l)
	#add any variant sites that were filtered for presumably not sharing snps
	filtered=$(cat $log_file | grep "Kept" | awk '{print $8}')
	let n_snps=n_snps+filtered
	#put all variables in a file
	echo -e "$n_shared\t$n_removed\t$allele_mismatches\t$n_snps" > variables.tmp
	# Cat the content of this file, prefixing each line with information on this
        # parameter combination.
        line_prefix="M$M-n$n-m$m\t$M\t$n\t$m\t"
        sed -r "s/^/$line_prefix/" variables.tmp >> CON_error_rates.tsv
	#remove the .tmp files
	rm *.tmp
	done
done

