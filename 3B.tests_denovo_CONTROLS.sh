#!/bin/bash

#uses the  control samples from floragenex
#to determine parameter values that may estimating locus, allele and snp error rates 
#between the two control samples which should be the same template

#set path to base dir
top=$(readlink -f $(dirname $0)/..)
#set popmap file
popmap=$top/INFO/popmap_controls.tsv
# there are only two samples and each plate is listed as a different population
#set reads dir
reads_dir=$top/cleanData/controls

#move to  directory for tests for control samples
cd $top/tests.denovo/control
#iterate over m values
m_values="3 4 5 6 7 8 9 10" 
#make subdirectories for different values of m
for m in $m_values; do
        mkdir -p stacks.m$m
done

for m in $m_values ;do
	 #iterate over M
        M_values="1 2 3 4 5 6 7 8"
        #create subdirectories for different values of  M
        for M in $M_values; do
                 mkdir -p stacks.m$m/stacks.M$M
        done

        for M in $M_values ;do
                n=$M
                echo "Running Stacks for m=$m M=$M, n=$n..."
                out_dir=stacks.m$m/stacks.M$M
                log_file=$out_dir/denovo_map.oe
                denovo_map.pl -T 2  --samples $reads_dir --popmap $popmap -o $out_dir  -M $M -n $n -m $m  &> $log_file
        #end of M/N values loop
        done

	# Check that all runs have completed.
        echo "Checking that all denovo_map runs have completed..."
        ls stacks.m$m/stacks.M*/denovo_map.oe | wc
        wc -l stacks.m$m/stacks.M*/denovo_map.oe
        grep -iE '\b(err|e:|warn|w:|fail|abort)' stacks.m$m/stacks.M*/denovo_map.oe stacks.m$m/stacks.M*/denovo_map.log
        grep -L 'denovo_map\.pl is done' stacks.m$m/stacks.M*/denovo_map.log

	# Run populations with '-p 2' (loci present on samples from both plates)
        for M in $M_values ;do
                stacks_dir=stacks.m$m/stacks.M$M
                out_dir=$stacks_dir/populations.p2
                mkdir -p $out_dir
                log_file=$out_dir/populations.oe
                populations -P $stacks_dir -O $out_dir -M $popmap -p 2 --fasta_samples_raw --vcf --fasta_samples &> $log_file
        done

#close m loop
done

#collect statistics on  #shared loci  #total loci #allele_mismatches
# #total_alleles #snpsbetweenreps #total variantsites
# to calculate  total shared loci, allele error rate and snp error rate
echo $top
cd $top/tests.denovo/control

#These error rate calculations are based on Mastretta-Yanes et al. 2014 Mol Ecol Resources

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


#Calculate the error rates and plot them for varyin m(at M=n=3) and varying M (at m=n=3)
echo "Calculating error rates and plotting the results ..."
Rscript $top/scripts/RScripts/5_plot_control_errorRates.R
