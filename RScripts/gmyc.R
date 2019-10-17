#GYMC for delineatign species from co1 data
#install.packages(c("ape", "paran", "rncl"))
#install.packages("splits", repos = "http://R-Forge.R-project.org")

#read in the tree
setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/CO1_BOLD/BEAST/bMT_CCS_sc/run1/")
library(ape)
coal_tr <- read.nexus("constant_coalescent_tree.nex")


#gmyc analaysis
library(splits)
coal_gmyc <- gmyc(coal_tr)

summary(coal_gmyc)

coal_support <- gmyc.support(coal_gmyc)        # estimate support
is.na(coal_support[coal_support == 0]) <- TRUE # only show values for affected nodes
plot(coal_tr, cex=.6, no.margin=TRUE)          # plot the tree
nodelabels(ifelse(coal_support>0.9, round(coal_support,2), "NA"), cex=.5)     # plot the support values on the tree

spec.list(coal_gmyc)