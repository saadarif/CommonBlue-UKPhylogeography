#!/bin/bash
#modified from https://bitbucket.org/rochette/rad-seq-genotyping-demo/raw
top=$(readlink -f $(dirname $0)/..)

#using 24 rep samples in popmap.test_samples.tsv

# STEP 15-A-i: Chose parameter combinations to survey.
m_values="3 4 5 6"

# STEP 15-A-ii: Change directory.
cd $top/tests.denovo

# STEP 15-A-iv: Run denovo_map on the subset of samples.
popmap=../INFO/popmap.test_samples.tsv
for m in $m_values;do
	
	#create all subdirectories for M/N
	M_values="1 2 3 4 5 6 7 8"
	for M in $M_values; do
		mkdir stacks.m$m.M$M
	done
	
	for M in $M_values ;do
		n=$M
		m=$m
		echo "Running Stacks for m=$m M=$M, n=$n..."
		reads_dir=../cleanData
		out_dir=stacks.m$m.M$M
		log_file=$out_dir/denovo_map.oe
		denovo_map.pl -T 6 --samples $reads_dir -O $popmap -o $out_dir -b 1 -S -M $M -n $n -m $m &> $log_file
	#end of M/N values loop
	done
#end of m values loop
done

# STEP 15-A-v: Check that all runs have completed.
echo "Checking that all denovo_map runs have completed..."
ls stacks.m*/denovo_map.oe | wc
wc -l stacks.m*/denovo_map.oe
grep -iE '\b(err|e:|warn|w:|fail|abort)' stacks.m*/denovo_map.oe stacks.m*/denovo_map.log
grep -L 'denovo_map\.pl is done' stacks.m*/denovo_map.log

# STEP 15-A-vi: Check coverages (in file `stacks.M1/denovo_map.log`).

# STEP 15-A-vii: Run populations with '-r 0.80' (loci present in 80% of samples)
for m in $m_values ;do	
	for M in $M_values ;do
		stacks_dir=stacks.m$m.M$M
		out_dir=$stacks_dir/populations.r80
		mkdir $out_dir
		log_file=$out_dir/populations.oe
		populations -P $stacks_dir -O $out_dir -r 0.80 &> $log_file
	done

# STEP 15-A-viii: Compare the results obtained with different parameters
mkdir results
cd results

# Extract the SNPs-per-locus distributions. These distributions
# are reported in the log of populations.
echo "Tallying the numbers..."
echo -e '#par_set\tM\tn\tm\tn_snps\tn_loci' > n_snps_per_locus.tsv
for m in $m_values ;do
	for M in $M_values ;do
		n=$M
		log_file=../stacks.m$m.M$M/populations.r80/batch_1.populations.log

		# Extract the numbers for this parameter combination.
		sed -n '/^#n_snps\tn_loci/,/^[^0-9]/ p' $log_file | grep -E '^[0-9]' > $log_file.snps_per_loc

		# Cat the content of this file, prefixing each line with information on this
		# parameter combination.
		line_prefix="M$M-n$n-m$m\t$M\t$n\t$m\t"
		sed -r "s/^/$line_prefix/" $log_file.snps_per_loc >> n_snps_per_locus.tsv
	#close M loop
	done
#close m loop
done
# Plot the results with R.
echo "Plotting the number of loci..."
Rscript $top/scripts/RScripts/4.plot_n_loci.R
echo "Plotting the distribution of the number of SNPs..."
Rscript $top/scripts/RScripts/4.plot_n_snps_per_locus.R
