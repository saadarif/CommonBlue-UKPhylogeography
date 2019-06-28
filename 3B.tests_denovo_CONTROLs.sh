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
                populations -P $stacks_dir -O $out_dir -p 2 &> $log_file
        done

#close m loop
done
