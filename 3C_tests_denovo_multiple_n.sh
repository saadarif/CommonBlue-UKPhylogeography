#!/bin/bash
#modified from https://bitbucket.org/rochette/rad-seq-genotyping-demo/raw
top=$(readlink -f $(dirname $0)/..)

#using 24 rep samples in popmap.test_samples.tsv

#Final parameters of testing
#for values of m=M=3:4 and n=n+1,n-1
#m=M=n=3:4 was already done previously
#try different values of m
m_values="3 4"

# STEP 15-A-iv: Run denovo_map on the subset of samples.
popmap=$top/INFO/popmap.test_samples.tsv
	
for m in $m_values; do
	#move to base dir
	cd $top/tests.denovo
	#iterate over M
	M=$m
	let nminus=m-1
	let nplus=m+1
	n_values="$nminus $nplus"

	for n in $n_values ;do
		#create subdirectories for different values of  M
                mkdir -p stacks.m$m/stacks.M$M.n$n
		echo "Running Stacks for m=$m M=$M, n=$n..."
		reads_dir=../cleanData
		out_dir=stacks.m$m/stacks.M$M.n$n
		log_file=$out_dir/denovo_map.oe
		denovo_map.pl -T 6  --samples $reads_dir --popmap $popmap -o $out_dir  -M $M -n $n -m $m  &> $log_file
	#end of M/N values loop
	done

	# STEP 15-A-v: Check that all runs have completed.
	echo "Checking that all denovo_map runs have completed..."
	ls stacks.m$m/stacks.M$M.n*/denovo_map.oe | wc
	wc -l stacks.m$m/stacks.M$M.n*/denovo_map.oe
	grep -iE '\b(err|e:|warn|w:|fail|abort)' stacks.m$m/stacks.M$M.n*/denovo_map.oe stacks.m$m/stacks.M$M.n*/denovo_map.log
	grep -L 'denovo_map\.pl is done' stacks.m$m/stacks.M$M.n*/denovo_map.log

	# STEP 15-A-vi: Check coverages (in file `stacks.M1/denovo_map.log`).

	# STEP 15-A-vii: Run populations with '-r 0.80' (loci present in 80% of samples)
	for n in $n_values ;do
		stacks_dir=stacks.m$m/stacks.M$M.n$n
		out_dir=$stacks_dir/populations.r80
		mkdir -p $out_dir
		log_file=$out_dir/populations.oe
		populations -P $stacks_dir -O $out_dir -r 0.80 &> $log_file
	done

	# STEP 15-A-viii: Compare the results obtained with different parameters
	mkdir -p stacks.m$m/results_var_n
	cd stacks.m$m/results_var_n
	
	#FIXME
	# Extract the SNPs-per-locus distributions. These distributions
	# are reported in the log of populations.
	#result should include the previous run of m=M=n
	n_values="$nminus $m $nplus"
	echo "Tallying the numbers..."
	echo -e '#par_set\tM\tn\tm\tn_snps\tn_loci' > n_snps_per_locus.tsv
	for n in $n_values ;do
		log_file=../stacks.M{3,4}*/populations.r80/populations.log.distribs

		# Extract the numbers for this parameter combination.
		stacks-dist-extract  $log_file snps_per_loc_postfilters | grep -E '^[0-9]' > $log_file.snps_per_loc

		# Cat the content of this file, prefixing each line with information on this
		# parameter combination.
		line_prefix="M$M-n$n-m$m\t$M\t$n\t$m\t"
		sed -r "s/^/$line_prefix/" $log_file.snps_per_loc >> n_snps_per_locus.tsv
		#close M loop
	done
	#plot results manually
 	#The current plotting scripts assume M=n
	# Plot the results with R.
	#echo "Plotting the number of loci..."
	#Rscript $top/scripts/RScripts/4.plot_n_loci.R $m
	#echo "Plotting the distribution of the number of SNPs..."
	#Rscript $top/scripts/RScripts/4.plot_n_snps_per_locus.R $m
	
#close the m loop
done
