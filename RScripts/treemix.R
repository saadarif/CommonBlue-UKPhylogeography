#plotting and analysing treemix results
#May2020

library(RColorBrewer)
library(R.utils)
library(OptM) #for testinng optimal migration events

source("/opt/treemixv1.13/src/plotting_funcs.R")

setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/Treemix/p15r50miss25SingleSNP/mOptimal/")

#Test optimal M
folder <- getwd()
optimalM <- optM(folder)
plot_optM(optimalM)
setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/Treemix/p15r50miss25SingleSNP_neutral////")


#choose prefix to plot
stem <- "p15r50miss25SingleSNP_neutral.treemix.in.gz.0"
plot_tree(stem, scale=T, mbar=F, ybar=0.9)

#check model residuals
#order populations by lat

plot_resid(stem)

setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/Treemix/p15r50miss25SingleSNP_noOut//")


#choose prefix to plot
stem <- "p15r50miss25SingleSNP.treemix.in.gz.1"
plot_tree(stem, scale=T, mbar=F, disp=0.001,ybar=0.1)
#flip=c(c(413,412))

#check model residuals
#order populations by lat

plot_resid(stem, "poporder")
