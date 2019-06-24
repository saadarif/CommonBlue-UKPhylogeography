#!/bin/bash

#visualzie erorr profiles of cleaned reads using kmer_filter

plate1_reads=/media/data_disk/PROJECTS/Saad/CommonBlue/cleanData/plate1
plate2_reads=/media/data_disk/PROJECTS/Saad/CommonBlue/cleanData/plate2

#concatenate all reads into single file for generating distribution
kmer_filter -f $plate1_reads/plate1_all.fq.gz -i 'gzfastq' --k_dist > $plate1_reads/kmer_dist_plate1

kmer_filter -f $plate2_reads/plate2_all.fq.gz -i 'gzfastq' --k_dist > $plate2_reads/kmer_dist_plate2

