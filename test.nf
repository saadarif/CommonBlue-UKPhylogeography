#!/usr/bin/env nextflow

files = Channel.fromPath( "/media/data_disk/PROJECTS/Saad/CommonBlue/cleanData/*/kmer_*" )
params.outdir = "/media/data_disk/PROJECTS/Saad/CommonBlue/scripts"

process wordcount {
   publishDir params.outdir, mode: 'copy', pattern: '*.txt'
    input:
    file x from files

    output:
    file '*file.txt' 

    """
     wc -l $x > ${x}_file.txt 
    """
}


