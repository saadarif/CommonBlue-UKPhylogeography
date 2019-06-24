#!/usr/bin/env nextflow

files = Channel.fromPath( "/media/data_disk/PROJECTS/Saad/CommonBlue/cleanData/*/kmer_*" )
params.outdir = "/media/data_dist/PROJECTS/Saad/CommonBlue/scripts"
process wordcount {
    publishDir params.outdir, saveAs:{ filename -> "foo_$filename" }
    input:
    file x from files

    output:
    file '*.txt'

    """
    < wc -l $x 
    """
}

result.println { it.trim() }

