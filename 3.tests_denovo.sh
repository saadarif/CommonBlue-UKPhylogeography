#!/bin/bash
#modified from https://bitbucket.org/rochette/rad-seq-genotyping-demo/raw
top=$(readlink -f $(dirname $0)/..)

#using 24 rep samples in popmap.test_samples.tsv

# STEP 15-A-i: Chose parameter combinations to survey.
M_values="1 2 3 4 5 6 7 8 9"

# STEP 15-A-ii: Change directory.
cd $top/tests.denovo

#create the subdirectories
for M in $M_values; do
     mkdir -p stacks.M$M
done
#this the only option for denovo_map.pl in stacks 2.0
m=3
# STEP 15-A-iv: Run denovo_map on the subset of samples.
popmap=../INFO/popmap.test_samples.tsv
	
for M in $M_values ;do
	n=$M
	echo "Running Stacks for m=$m M=$M, n=$n..."
	reads_dir=../cleanData
	out_dir=stacks.M$M
	log_file=$out_dir/denovo_map.oe
	denovo_map.pl -T 6  --samples $reads_dir --popmap $popmap -o $out_dir  -M $M -n $n  &> $log_file
	#end of M/N values loop
	done
#end of m values loop
done

# STEP 15-A-v: Check that all runs have completed.
echo "Checking that all denovo_map runs have completed..."
ls stacks.M*/denovo_map.oe | wc
wc -l stacks.M*/denovo_map.oe
grep -iE '\b(err|e:|warn|w:|fail|abort)' stacks.M*/denovo_map.oe stacks.M*/denovo_map.log
grep -L 'denovo_map\.pl is done' stacks.M*/denovo_map.log

# STEP 15-A-vi: Check coverages (in file `stacks.M1/denovo_map.log`).

# STEP 15-A-vii: Run populations with '-r 0.80' (loci present in 80% of samples)
for M in $M_values ;do
	stacks_dir=stacks.M$M
	out_dir=$stacks_dir/populations.r80
	mkdir -p $out_dir
	log_file=$out_dir/populations.oe
	populations -P $stacks_dir -O $out_dir -r 0.80 &> $log_file
	done

# STEP 15-A-viii: Compare the results obtained with different parameters
mkdir -p results
cd results

# Extract the SNPs-per-locus distributions. These distributions
# are reported in the log of populations.
echo "Tallying the numbers..."
echo -e '#par_set\tM\tn\tm\tn_snps\tn_loci' > n_snps_per_locus.tsv
for M in $M_values ;do
	n=$M
	log_file=../stacks.M$M/populations.r80/populations.log.distribs

	# Extract the numbers for this parameter combination.
	stacks-dist-extract  $log_file snps_per_loc_postfilters | grep -E '^[0-9]' > $log_file.snps_per_loc

	# Cat the content of this file, prefixing each line with information on this
	# parameter combination.
	line_prefix="M$M-n$n-m$m\t$M\t$n\t$m\t"
	sed -r "s/^/$line_prefix/" $log_file.snps_per_loc >> n_snps_per_locus.tsv
	#close M loop
done

# Plot the results with R.
echo "Plotting the number of loci..."
Rscript $top/scripts/RScripts/4.plot_n_loci.R
echo "Plotting the distribution of the number of SNPs..."
Rscript $top/scripts/RScripts/4.plot_n_snps_per_locus.R
