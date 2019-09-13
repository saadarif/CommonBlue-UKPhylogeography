library(gghybrid)

#make sure the infile below only has either spaces or tabs between fields. The conversion from PGDspider does not guarantee this
gg_snps <- read.data("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4/populations.r50.p15_moh_0.65/fastStructure/gghybrid_in"
                  , MISSINGVAL=-9, MARKERNAME = 0, nprecol = 2, precol.headers = 0, NUMLOCI=2699)

#Data preparation and filtering. Here I'm filtering out loci that have a minor allele
#frequency greater than 0.1 in both parental reference sets. There are also options for
#filtering by difference in parental allele frequencies, and for number of allele copies
#in each parental reference set (variable among loci due to missing data).

#The function uses objects produced by 'read.data'#

prepdata=data.prep(data=gg_snps$data, loci=gg_snps$loci, alleles=gg_snps$alleles, S0=c("BMD", "ETB" ,"MMS", "RNL","PCP","BWD","MDC"), S1=c("TUL", "BER"), precols=gg_snps$precols, max.S.MAF = 0.1, return.genotype.table=T, return.locus.table=T)
#'return.genotype.table=T' makes an optional table of genotypes, where for each locus
#an individual's genotype (assuming diploidy) will be 0 (two copies of the allele with
#relatively higher frequency in the 'S0' parental set), 1 (heterozygote), 2 (two copies
#of the designated 'S1' allele). This table isn't needed in downstream functions, but
#could be useful e.g. for estimating parental linkage disequilibria (associations of
#alleles from the same parent species).

#'return.locus.table=T' is also optional and not needed downstream. It's just a table
#with one row per marker, giving some information on parental allele frequencies, sample
#sizes etc.

###

#Next, run hybrid index estimation#

#This function uses objects produced by both the previous functions#

hindlabel= esth(data.prep.object = prepdata$data.prep,
                read.data.precols = gg_snps$precols,
                include.Source = TRUE,	#Set to TRUE if you want hybrid indices for the parental reference individuals#
                plot.ind = c("BER_m_027","OBN_m_119","MLG_m_012","RHD_m_241",
                             "FRN_m_M07","BMD_f_149"),
                plot.col = c("blue","green","yellow","purple","pink","red"),
                nitt=10000,burnin=5000)

#The plots ('plot.ind' and 'plot.col' options) are optional. They plot accepted posterior hybrid 
#index values in real time, in this case for 5 randomly chosen individuals.

#Plot a subset of the estimated hybrid indices (the resulting object 'abc' is useful for making a legend)#

setkey(hindlabel$hi,POPID)	#function from data.table, for rapid sorting and subsetting#

#
abc = plot_h(data=hindlabel$hi[c("MLG", "MDC", "BER", "TUL", "DGC", "RVS", "OBN", "BWD", "FRN", "PCP", "ETB", "MMS", "RNL", "RHD", "BMD")],#All of POPIDs#
             test.subject=hindlabel$test.subject,
             mean.h.by="POPID",			#Calculate the mean hybrid index for each value of the "POPID" column#
             sort.by=c("mean_h","POPID","h_posterior_mode"),	#Order test subjects along the x axis by the mean hybrid 
             #index calculated above and also by individual hybrid index#
             col.group="POPID",
             group.sep="POPID",
             fill.source=TRUE,
             basic.lines=FALSE,
             source.col=c("blue","red"),
             source.limits=c("blue","red"),
             cex=1,pch=16,
             cex.lab=1.5,cex.main=1.5,ylim=c(0,1))
#

#Reshape the plot window as you want#

#Add a legend using the 'plot_h' object 'abc'#

setkey(abc,rn)		#Order data by row number#
legend("topleft",	#Place the legend in the top left of the figure#
       abc[,POPID], 		#Name of the field by which data point colours are grouped#
       bg="white",			#Background colour#
       text.col=c("black"), #Text colour#
       pch=22, 				#Text size#
       col=abc[,col.Dark2], #Name of the field containing colour information#
       pt.bg=abc[,col.Dark2],	#Name of the field containing colour information#
       ncol=2,				#Number of columns for the group names#
       cex=1, pt.cex=1)

###
