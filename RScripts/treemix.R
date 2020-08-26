#plotting and analysing treemix results
#May2020

library(RColorBrewer)
library(R.utils)
source("/opt/treemixv1.13/src/plotting_funcs.R")

setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/Treemix/p15r50miss25SingleSNP_neutral/")

#choose prefix to plot
stem <- "p15r50miss25SingleSNP_neutral.treemix.in.gz.1"
plot_tree(stem, scale=T)

#check model residuals
#order populations by lat

plot_resid(stem)
