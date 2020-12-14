#!/bin/bash

# a script to extract fasta sequences from
# a list of STACK loci (numbers only)

#$1 = a list of stack loci with numbers only
#$2 = the fasta file to extract sequences from

#add the "CLocus_" term in front of all loci numbers for search
pattern="CLocus_"
sed "s/^/$pattern/g" $1 > loci.txt

#use grep to extract the sequence from the new file 
#-A1 assumes only a single line of sequence after name
grep -w -A 1 -f  loci.txt $2 --no-group-separator > out_loci.fa
