#!/bin/bash
top=$(readlink -f $(dirname $0)/..)

# STEP 16, STEP 17-A-i: Change directory.
cd $top/stacks.denovo

# STEP 17-A-ii: Run ustacks on every sample, use m=3 M=n=4 based on r80 rule from tests denovo
# M=n=4 aslo results in the lowest allele and snp error rates of the control samples.
#However in the controls m=4 recovered the most shared loci (marignally) even though total loci was less than m3
m=5
M=4
n=4
index=1

mkdir -p stacks_m${m}_M${M}_n${n}
cd stacks_m${m}_M${M}_n${n}
#NOTE the following samples have been excluded from popmap.tsv
#due to <500,000 reads
#OBN_m_114 
#RNL_m_225
#RNL_m_233


for sample in $(cut -f1 ../../INFO/popmap.tsv) ;do
	echo "Running ustacks on $sample.."
	fq_file=../../cleanData/$sample.fq.gz
	logfile=$sample.ustacks.oe
	ustacks -f $fq_file -i $index -o ./ -M $M -m $m -p 4 &> $logfile
	index=$(( $index + 1 ))
done

# STEP 17-A-iii: Check that all ustacks runs have completed.
echo "Checking that all ustacks runs have completed..."
ls *.ustacks.oe | wc
wc -l *.ustacks.oe
grep -iE '\b(err|e:|warn|w:|fail|abort)' *.ustacks.oe
grep -L 'ustacks is done' *.ustacks.oe

# STEP 17-A-iv: Extract the sample coverages from the `*.ustacks.oe` files
# and check the values are consistent with previously obtained coverages.

# STEP 17-A-v: Pick the samples to include in the catalog.
# For each population, we pick the 10 samples with the highest coverage excluding any
#samples with abnormally high reads or abnornmally low reads
#for pops with less than 10 add all individuals
awk '$3=="y" {print $0}' ../../INFO/n_reads_per_sample.tsv | cut -f1,4 > ../../INFO/popmap.catalog.tsv

#n_per_pop=10
#excluded_samples=sj_1484.07
#for pop in $(cut -f2 ../info/popmap.tsv | sort -u) ;do
#	sort -k2,2nr ../info/n_reads_per_sample.tsv \
#		| grep -v $excluded_samples \
#		| grep "^$pop" \
#		| head -n $n_per_pop \
#		| cut -f1 \
#		| sed -r 's/([A-Za-z]+)(.+)/\1\2\t\1/'

# STEP 17-A-vi: Run cstacks
echo "Running cstacks..."
cstacks -P ./ -M ../../INFO/popmap.catalog.tsv -n $n -p 4 &> cstacks.oe

# STEP 17-A-vii: Run sstacks on every sample
sstacks -P ./ -M ../../INFO/popmap.tsv  -p 4 &> sstacks.oe

# (Check that all sstacks runs have completed.)
# Run tsv2bam to transpose the data so it is stored by locus, instead of by sample. We will include

echo "Running tsvbam ..."
tsv2bam -P ./ -M ../../INFO/popmap.tsv  -t 4

#Run gstacks: build a paired-end contig from the metapopulation data (if paired-reads provided),
# align reads per sample, call variant sites in the population, genotypes in each individual.
# NOT SURE IF THIS IS REQUIRED EITHER
echo "Running gstacks ..."
gstacks  -P ./ -M ../../INFO/popmap.tsv  -t 4

# STEP 21: Filter genotypes with populations.

#filter at r 0.8
mkdir -p populations.r75.p15_moh_0.65
min_samples=0.75
min_pops=15
min_maf=0.05
max_obs_het=0.65
populations -P ./ -M ../../INFO/popmap.tsv -r $min_samples -p $min_pops --min-maf=$min_maf --max-obs-het=$max_obs_het \
	   -O populations.r75.p15_moh_0.65 --fstats \
	   --vcf --genepop --radpainter -t 4 &> populations.r75.p15_moh_0.65/populations.oe
#filter at r 0.7
mkdir -p populations.r70.p15_moh_0.65
min_samples=0.70
populations -P ./ -M ../../INFO/popmap.tsv -r $min_samples -p $min_pops --min-maf=$min_maf --max-obs-het=$max_obs_het \
           -O populations.r70.p15_moh_0.65 --fstats --vcf --genepop --radpainter -t 4 &> populations.r70.p15_moh_0.65/populations.oe
#filter at r 0.6
mkdir -p populations.r60.p15_moh_0.65
min_samples=0.60
populations -P ./ -M ../../INFO/popmap.tsv -r $min_samples -p $min_pops --min-maf=$min_maf --max-obs-het=$max_obs_het \
           -O populations.r60.p15_moh_0.65 --fstats --vcf --genepop --radpainter -t 4 &> populations.r60.p15_moh_0.65/populations.oe
#filter at 0.5
mkdir -p populations.r50.p15_moh_0.65
min_samples=0.50
populations -P ./ -M ../../INFO/popmap.tsv -r $min_samples -p $min_pops --min-maf=$min_maf --max-obs-het=$max_obs_het \
           -O populations.r50.p15_moh_0.65 --fstats --vcf --genepop --radpainter -t 4 &> populations.r50.p15_moh_0.65/populations.oe

#adittional pop filters
mkdir -p populations.r50.p11_moh_0.65
populations -P ./ -M ../../INFO/popmap.tsv -r 0.5 -p 11 --min-maf=0.05 --max-obs-het=0.65 \
           -O populations.r50.p11_moh_0.65 --fstats --vcf --phylip --phylip-var -t 4 &> populations.r50.p11_moh_0.65/populations.oe

mkdir -p populations.r50.p8_moh_0.65
populations -P ./ -M ../../INFO/popmap.tsv -r 0.5 -p 8 --min-maf=0.05 --max-obs-het=0.65 \
           -O populations.r50.p8_moh_0.65 --fstats --vcf --phylip --phylip-var -t 4 &> populations.r50.p8_moh_0.65/populations.oe
