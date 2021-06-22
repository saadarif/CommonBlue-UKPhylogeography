#!/usr/bin/env Rscript
#from https://bitbucket.org/rochette/rad-seq-genotyping-demo/raw/c22c0c66b8fd88b47647b697bec2e7e5f3cbfb80/demo_scripts/R_scripts

d = read.delim('./n_reads_per_sample.tsv')

pdf('/media/data_disk/PROJECTS/Saad/CommonBlue/plots/n_reads_per_sample.pdf', height=12, width=12)
dotchart(sort(setNames(d$n_reads, d$X.sample)), xlim=c(0, max(d$n_reads)))
hist(d$n_reads, nclass=20)
