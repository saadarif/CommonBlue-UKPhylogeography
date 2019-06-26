#!/usr/bin/env Rscript

d = read.delim('./n_snps_per_locus.tsv')
# Keep only M==n, and loop through values of m:3
d = subset(d, M==n & m==3)
# Make sure the table is ordered by number of snps.
d = d[order(d$n_snps),]

Mn_values = sort(unique(d$M))

# Write the counts in a matrix.
m = matrix(NA, nrow=length(Mn_values), ncol=max(d$n_snps)+1)
for(i in 1:nrow(d)) {
	m[d$M[i],d$n_snps[i]+1] = d$n_loci[i] # [n_snps+1] as column 1 is for loci with 0 SNPs
}

# Truncate the distributions.
max_n_snps = 10
m[,max_n_snps+2] = rowSums(m[,(max_n_snps+2):ncol(m)], na.rm=T)
m = m[,1:(max_n_snps+2)]
m = m / rowSums(m, na.rm=T)

# Draw the barplot.
pdf('./n_snps_per_locus.pdf')

col = rev(heat.colors(length(Mn_values)))

barplot(m,
	beside=T, col=col, las=1,
	names.arg=c(0:max_n_snps, paste('>', max_n_snps, sep='')),
	xlab='Number of SNPs',
	ylab='Percentage of loci',
	main='Distributions of the number of SNPs per locus\nfor a range of M==n values'
	)
legend('topright', legend=c('M==n', Mn_values), fill=c(NA, col))

null=dev.off()
