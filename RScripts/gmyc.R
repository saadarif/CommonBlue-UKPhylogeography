#GYMC for delineatign species from co1 data
#install.packages(c("ape", "paran", "rncl"))
#install.packages("splits", repos = "http://R-Forge.R-project.org")

#read in the tree
setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/CO1_BOLD/BEAST/bMT_CCS_sc_nosplit/")
library(ape)
coal_tr <- read.nexus("run1_run2_combined_mcc")


#gmyc analaysis
library(splits)
coal_gmyc <- gmyc(coal_tr, method = "single", interval = c(0, 10), quiet =
                    FALSE)

summary(coal_gmyc)

coal_support <- gmyc.support(coal_gmyc)        # estimate support
is.na(coal_support[coal_support == 0]) <- TRUE # only show values for affected nodes
plot(coal_tr, cex=.6, no.margin=TRUE)          # plot the tree
nodelabels(round(coal_support, 2), cex=.7)    # plot the support values on the tree

spec.list(coal_gmyc)
