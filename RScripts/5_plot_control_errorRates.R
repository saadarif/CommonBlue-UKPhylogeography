#!/usr/bin/env Rscript

#error rate calculation were inspired by Mastretta-Yanes et al 2014 Mol Ecol resources

d = read.delim('./CON_error_rates.tsv')

#add total loci here tl = shared = removed
d$total_loci=d$shared_loci+d$removed_loci

#add allele error rate where ae = allele mismatches/total shared loci
d$allele_error <- d$allele_mis/d$shared_loci

#add snp error rate where se = n_snps/total shared loci
d$snp_error <- d$n_snps/d$shared_loci

#plot the shared and total loci for different values of m while M=N=2 , default M=2
d_1 = subset(d, M==2)

# Make sure the table is ordered
d_1 = d_1[order(d_1$m),]

#save the output
pdf('./control_errorRates.pdf')
# Allele, Snp error rates Number of loci as m increases
# ==========

plot(NULL,
     xlim=range(d_1$m),
     ylim=range(c(0, d_1$total_loci)),
     xlab='m , M=n=2',
     ylab='Number of loci',
     main='Number of loci as m increases',
     xaxt='n',
     las=2
)
abline(h=0:20*1000, lty='dotted', col='grey50')
axis(1, at=c(1:10))
legend('bottomright', c('Total loci', 'shared loci'), lty=c('solid', 'dashed'))

# Total number of loci.
lines(d_1$m, d_1$total_loci)
points(d_1$m, d_1$total_loci, cex=0.5)

# Number of shared loci.
lines(d_1$m, d_1$shared_loci, lty='dashed')
points(d_1$m, d_1$shared_loci, cex=0.5)

#--------------------------------------------------
#allele and snp error rates
plot(NULL,
     xlim=range(d_1$m),
     ylim=c(0, 0.15),
     xlab='m , M=n=2',
     ylab='Error rates',
     main='SNP and Allele error rates as m increases',
     xaxt='n',
     las=2
)
abline(h=c(0.025, 0.05, 0.075, 0.1, 0.125), lty='dotted', col='grey50')
axis(1, at=c(1:10))
legend('topright', c('Allele Error rate', 'SNP Error rate'), lty=c('solid', 'dashed'))

#allele error rates.
lines(d_1$m, d_1$allele_error)
points(d_1$m, d_1$allele_error, cex=0.5)

#SNP error rate
lines(d_1$m, d_1$snp_error, lty='dashed')
points(d_1$m, d_1$snp_error, cex=0.5)


#========================================================================
#Allele, Snp error rates Number of loci as M increases
#plot the shared and total loci for different values of M while m=N=2
d_1 = subset(d, m==3)

# Make sure the table is ordered
d_1 = d_1[order(d_1$M),]

plot(NULL,
     xlim=range(d_1$M),
     ylim=range(c(0, d_1$total_loci)),
     xlab='M , m=n=3',
     ylab='Number of loci',
     main='Number of loci as M increases',
     xaxt='n',
     las=2
)
abline(h=0:20*1000, lty='dotted', col='grey50')
axis(1, at=c(1:8))
legend('bottomright', c('Total loci', 'shared loci'), lty=c('solid', 'dashed'))

# Total number of loci.
lines(d_1$M, d_1$total_loci)
points(d_1$M, d_1$total_loci, cex=0.5)

# Number of shared loci.
lines(d_1$M, d_1$shared_loci, lty='dashed')
points(d_1$M, d_1$shared_loci, cex=0.5)

#-----------------------------------------------------
#allele and snp error rates
plot(NULL,
     xlim=range(d_1$M),
     ylim=c(0, 0.15),
     xlab='M , m=n=3',
     ylab='Error rates',
     main='SNP and Allele error rates as M increases',
     xaxt='n',
     las=2
)
abline(h=c(0.025, 0.05, 0.075, 0.1, 0.125), lty='dotted', col='grey50')
axis(1, at=c(1:10))
legend('topright', c('Allele Error rate', 'SNP Error rate'), lty=c('solid', 'dashed'))

#allele error rates.
lines(d_1$M, d_1$allele_error)
points(d_1$M, d_1$allele_error, cex=0.5)

#SNP error rate
lines(d_1$M, d_1$snp_error, lty='dashed')
points(d_1$M, d_1$snp_error, cex=0.5)


null=dev.off()
