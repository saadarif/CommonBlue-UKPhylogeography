#!/bin/bash

#extract per sample locus coverage after running denovo_map.pl

#assumes Pop is the first prefix seperated by "_"

#define the directory prefix of multiple runs of denovo_map.pl
dir=""

#extract the coverage (at the end of ustacks)
sed -n '/^Depths of Coverage for Processed Samples\:/,/^$/ p' $dir/denovo_map.log |sed '/^$/d' | \
	sed '1d'| awk -F ":" '{print $2}' | sed -e 's/^[ \t]*//' | sed -e 's/x//' > cov_temp
#extract sample names
sed -n '/^Depths of Coverage for Processed Samples\:/,/^$/ p' denovo_map.log |sed '/^$/d' | \
	 sed '1d' | awk -F ":" '{print $1}' > samples_temp
#extract pop names as prefix  before first "_"
sed -n '/^Depths of Coverage for Processed Samples\:/,/^$/ p' denovo_map.log |sed '/^$/d' | \
	 sed '1d' | awk -F ":" '{print $1}' | awk -F "_" '{print $1}' > pop_temp
#combine all the files and then remove them
paste samples_temp pop_temp cov_temp >> coverage_stats
#add a line prefix to this file and then append to a file with coverage from all runs
