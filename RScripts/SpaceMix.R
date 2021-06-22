#SpaceMix

#Dec 2020

#Space Mix analysis
#https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005703#pgen-1005703-g008

#VCF file was converted to input using plus commits
#https://github.com/gbradburd/SpaceMix/pull/4/commits/afcbb99782ee7d8c241080c6c654a6afd7ec1035

library(vcfR)
library(dartR)
library(tidyverse)

#read in vcf to get coordinates
MLG=c(-5.834874, 56.991962)# check field notes 16
MDC=c(-3.843251, 53.295675) #12
BER=c(-7.213534, 57.713854)#16
TUL=c(-7.025334, 58.185530)# check field notes 16
DGC=c(-4.016982, 57.878096)#14
RVS=c(-2.596053, 56.023971)#6
OBN=c(-5.482266, 56.433569)#14
BWD=c(-3.898247, 50.579755) #devon 7
FRN=c(1.982669, 44.157335) #Will's site in France 6
PCP=c(-4.308992, 51.674219) #Pembrey county park 13
ETB=c( 0.243747, 50.769187) #eastbourne 14
MMS=c(1.255096, 52.168021) #martin's meadows in suffolk 14
CFW=c(-0.281520, 53.254429) #should be Chamber's Farm wood CFW , lincolnshire 10
RHD=c(-1.477433, 54.713907) #Raisby hill durham 13
BMD=c(-1.125439, 51.782563) #Oxford 16


coords <- rbind(MLG, MDC, BER, TUL, DGC, RVS, OBN, BWD, FRN, PCP, ETB, MMS, CFW, RHD, BMD)
colnames(coords) <- c("lon", "lat")
coords <- as.data.frame(coords)
#convert rownames to column
coords <- tibble::rownames_to_column(coords, "POP")


VCF <- read.vcfR("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/populations.snps.filter2.0.25.recode.neutral.vcf")
z <- vcfR2genlight(VCF)
popnames <- z@ind.names
pops <- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[1])}))
pops<- str_replace_all(pops, "RNL", "CFW")
pop(z) = pops
ploidy(z) <- 2
latlong = data.frame( pop=z@pop, lat=numeric(length=length(z@pop)), lon=numeric(length=length(z@pop)))

for (i in 1:dim(latlong)[1]){
  print(i)
  latlong$lat[i] <- coords[coords$POP==latlong[i,1],2]
  latlong$lon[i] <- coords[coords$POP==latlong[i,1],3]
  
}

colnames(latlong) <- c("pop","lat", "lon")
#add lat long to gen light object
z@other$latlong <- latlong[,2:3]
#should now be able to use the dartR IBD function
#add sex information as well
popnames <- z@ind.names 
sex<- as.factor(sapply(strsplit(popnames, '_'), function(x){paste0(x[2])}))
ind <- cbind( z@ind.names, pops, sex, latlong[,2:3])
ind$`z@ind.names` <- as.character(ind$`z@ind.names` )
colnames(ind)[1] <- "ind"
z@other$ind.metrics <- ind

setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/SpaceMix/")

library(SpaceMix)
vignette("spacemix_vignette")
# read in the data and combine
counts <- as.matrix(read.table("populations.snps.filter2.0.25.recode.singlesnp_neut.vcf.counts"))
size <- as.matrix(read.table("populations.snps.filter2.0.25.recode.singlesnp_neut.vcf.size"))

rownames(counts) <- NULL
colnames(counts) <- NULL
rownames(size) <- NULL
colnames(size) <- NULL
sm_data <- list(allele.counts=as.matrix(counts), sample.size=as.matrix(size), population.coordinates=as.matrix(cbind(z@other$latlong[,2],z@other$latlong[,1])))
str(sm_data)

#Running the Spacemix analysis
# Data option: allele counts and sample sizes
# Fast Model option: estimating geogenetic locations
# Long Model option: estimating geogenetic locations and 
#                    admixture source locations
# Spatial priors: default variance,
#                   observed geographic sampling locations
run.spacemix.analysis(n.fast.reps = 5,
                      fast.MCMC.ngen = 5e6,
                      fast.model.option = "source_and_target",
                      long.model.option = "source_and_target",
                      data.type = "counts",
                      sample.frequencies=NULL,
                      mean.sample.sizes=NULL,
                      counts = sm_data$allele.counts,
                      sample.sizes = sm_data$sample.size,
                      sample.covariance=NULL,
                      target.spatial.prior.scale=NULL,
                      source.spatial.prior.scale=NULL,
                      spatial.prior.X.coordinates = sm_data$population.coordinates[,1],
                      spatial.prior.Y.coordinates = sm_data$population.coordinates[,2],
                      round.earth = T,
                      long.run.initial.parameters=NULL,
                      k = nrow(sm_data$allele.counts),
                      loci = ncol(sm_data$sample.size),
                      ngen = 1e8,
                      printfreq = 1e5,
                      samplefreq = 1e5,
                      mixing.diagn.freq = 50,
                      savefreq = 1e7,
                      directory="SingleSnpNeut_mod4_run1",
                      prefix = "SingleSnpNeut_mod4_run1")



#the following run is for an analysis not estimating admixture
#Running the Spacemix analysis
# Data option: allele counts and sample sizes
# Fast Model option: estimating geogenetic locations
# Long Model option: estimating geogenetic locations 
# Spatial priors: default variance,
#                   observed geographic sampling locations
setwd("/media/data_disk/PROJECTS/Saad/CommonBlue/stacks.denovo/stacks_m4_M4_n4_new/populations.r50.p15_moh_0.65/SpaceMix/")

run.spacemix.analysis(n.fast.reps = 5,
                      fast.MCMC.ngen = 5e6,
                      fast.model.option = "target",
                      long.model.option = "target",
                      data.type = "counts",
                      sample.frequencies=NULL,
                      mean.sample.sizes=NULL,
                      counts = sm_data$allele.counts,
                      sample.sizes = sm_data$sample.size,
                      sample.covariance=NULL,
                      target.spatial.prior.scale=NULL,
                      source.spatial.prior.scale=NULL,
                      spatial.prior.X.coordinates = sm_data$population.coordinates[,1],
                      spatial.prior.Y.coordinates = sm_data$population.coordinates[,2],
                      round.earth = T,
                      long.run.initial.parameters=NULL,
                      k = nrow(sm_data$allele.counts),
                      loci = ncol(sm_data$sample.size),
                      ngen = 1e8,
                      printfreq = 1e5,
                      samplefreq = 1e5,
                      mixing.diagn.freq = 50,
                      savefreq = 1e7,
                      directory="./SingleSnpNeut_mod2_run1",
                      prefix = "SingleSnpNeut_mod2_run1")
#--------------------------------------------------------------------------------------
#trace plots
load("SingleSnpNeut_mod2_run1/SingleSnpNeut_mod2_run1_LongRun/SingleSnpNeut_mod2_run1_space_MCMC_output1.Robj")
plot(Prob,xlab="MCMC iterations",ylab="value",
     main="Posterior probability trace plot",type='l')
# Trace plots of alpha parameters of the spatial covariance function
matplot(t(nugget),type='l',
        xlab="MCMC iterations",ylab="Parameter value",
        main="Trace plot of nugget parameters")



#----------------------------------------------
#joint marginals
# Joint marginal plot of a0 and a1
#   colored by where in the MCMC these 
#   parameters took their values
plot(a0,a1,xlab="a0",ylab="a1",
     main="Joint marginal of a0 and a1",pch=20,
     col=adjustcolor(rainbow(1000,start=4/6,end=6/6),0.3))
legend(x="bottomright",pch=19,cex=0.8,
       col=rainbow(1000,start=4/6,end=6/6)[c(1,500,1000)],
       legend=c("Sampled MCMC iteration 1",
                "Sampled MCMC iteration 500",
                "Sampled MCMC iteration 1000"))

#----------------------------------------------------------------------------------------
#acceptance rate plots
plot(accept_rates$a0_accept_rate,
     xlab="MCMC iterations",ylab="Acceptance rate",
     main="Acceptance rate of a0",type='l',
     ylim=c(0.35,0.6))
abline(h=0.44,col="gray",lty=2)
# Acceptance rates of nugget parameters over the 
#   course of the MCMC analysis
matplot(t(accept_rates$nugget_accept_rate),
        xlab="MCMC iterations",ylab="Acceptance rate",
        main="Acceptance rates of nuggets",type='l',
        ylim=c(0.3,0.7))
abline(h=0.44,col="gray",lty=2)

#-----------------------------------------------------------------------------
#Evaluating model adequacy
# first, load the standardized (mean-centered and normalized)
#   allele frequency data object.  This object, which is the 
#   "MCN.frequencies.list" (Mean Centered and Normalized) is 
#   saved in the Long Run directory, and is generated if the 
#   user has specified either allele count or allele frequeny 
#   data. 
#   Note that it is not generated if the user has specified the 
#   sample covariance.
load("SingleSnpNeut_mod4_run1/SingleSnpNeut_mod4_run1_LongRun/SingleSnpNeut_mod4_run1_MCN.frequencies.list.Robj")
# Now, calculate the sample covariance from the mean centered 
#   and normalized sample allele frequencies.
sample.covariance <- cov(t(MCN.frequencies.list$mean.centered.normalized.sample.frequencies),
                         use="pairwise.complete.obs")


# Create a matrix that will perform a mean-centering 
#   on the parametric covariance matrix
# Then, mean-center the parametric ovariance matrix.
k <- nrow(MCN.frequencies.list$mean.centered.normalized.sample.frequencies)
MC.matrix <- diag(k) - matrix(1/last.params$inv.mean.sample.sizes / 
                                (sum(1/last.params$inv.mean.sample.sizes)),
                              nrow=k,ncol=k,byrow=TRUE)

MC.parametric.covariance <- (MC.matrix) %*%     
 last.params$admixed.covariance %*% 
  t(MC.matrix)
# Finally, compare the sample covariance to the parametric
#   covariance.  Ideally, there will be a very tight correspondence 
#   between the data and the model.  If there is not, it may 
#   be an indication either that the MCMC has not converged on 
#   the stationary distribution or that the process that generated 
#   the data is only poorly approximated by SpaceMix's model.

# The sample and parametric covariances can be plotted 
#   against each other (if model fit is good they should 
#   fall on the x=y red line)
index.matrix <- upper.tri(sample.covariance,diag=TRUE)
plot(sample.covariance[index.matrix], 
     MC.parametric.covariance[index.matrix],
     col=adjustcolor("black",0.3),pch=20,
     xlab="Sample covariance",
     ylab="Parametric covariance",
     main="Model adequacy:\n matrix comparison")
abline(0,1,col="red")

# Or the patterns of decay of covariance with 
#   geographic distance can be compared between 
#   the data and the model.
plot(last.params$D[1:k,1:k][index.matrix], 
     sample.covariance[index.matrix],
     pch=19,col="black",
     xlab="geogenetic distance",
     ylab="covariance",
     main="Model adequacy:\n IBD patterns")
points(last.params$D[1:k,1:k][index.matrix], 
       MC.parametric.covariance[index.matrix],col="red",pch=20)
legend(x="topright",pch=19,col=c(1,2),
       legend=c("observed","model estimate"))

#----------------------------------------------------------
#plotting the geogenetic map

#generate a vector of sample colors
sample.colors <- rainbow(n=148,start=4/6,end=6/6)[as.numeric(cut(sm_data$population.coordinates[,1],148))]
# And now generate a sample map list using a 95% 
#   credible interval on parameter estimates without 
#   `burning' (i.e., discarding) any sampled iterations
#   of the MCMC.
example.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "SingleSnpNeut_mod2_run1/SingleSnpNeut_mod2_run1_LongRun/SingleSnpNeut_mod2_run1_space_MCMC_output1.Robj",
                                                    geographic.locations = sm_data$population.coordinates,
                                                    name.vector = z@ind.names,
                                                    color.vector = sample.colors,
                                                    quantile=.95,
                                                    burnin=0)

# Now we generate a map of the output showing sample names 
#   at the locations of the maximum a posteriori (MAP) 
#   geogenetic location parameter estimates
make.spacemix.map(spacemix.map.list = example.spacemix.map.list,
                  text=TRUE,
                  ellipses=FALSE,
                  source.option=FALSE)
