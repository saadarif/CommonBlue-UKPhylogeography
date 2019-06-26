#!/bin/bash
#This script is modified from Rochettet and Catchen 2017 Nat Prot
#https://bitbucket.org/rochette/rad-seq-genotyping-demo/src/default/demo_scripts

#set project root directory
top=$(readlink -f $(dirname $0)/..)

# STEP 9: Check the per-sample coverages.

# Extract the number of reads from the log of process_radtags.
cd $top/INFO

echo -e '#sample\tn_reads' > n_reads_per_sample.tsv
for plate in plate1 plate2 ;do
	# Retrieve the part of the log between 'Barcode...' and the next empty line,
	# then discard the first and last lines and keep the 2nd and 6th columns.
	sed -n '/^Barcode\tFilename\t/,/^$/ p' ../cleanData/$plate/process_radtags.rawdata.log \
		| sed '1 d; $ d' \
		| cut -f2,6 \
		>> n_reads_per_sample.tsv
done

# Plot these numbers.
echo "Plotting per-sample coverage..."
Rscript ../scripts/RScripts/2.plot_n_reads_per_sample.R


#The following samples were removed for having reads <500,00
#RNL_m_225, RNL_m_232, OBN_m114
#The folloing samples were removed for having reads >6,000,000
#BER_m_028, BER_m_029, DGC_M_100, MMS_m_208
