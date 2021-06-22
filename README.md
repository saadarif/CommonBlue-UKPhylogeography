## BASH/R scripts and processed data files for analysis from:

>The dirty north: Evidence for multiple colonisations and *Wolbachia* infections shaping the genetic structure of the widespread butterfly *Polyommatus icarus* in the British Isles.
>Saad Arif, Michael Gerth, William G. Hone-Millard, Maria D. S. Nunes, Leonardo Dapporto, Timothy G. Shreeve.
>*bioRxiv*, 2020.09.03.267203; doi: https://doi.org/10.1101/2020.09.03.2672

1. BASH Scripts in the top folder were used for calling and filtering genotypes in both *P. icarus*
and the Wolbachia that infect them. The Stacks genotyping scripts in the top folder are largely based on scripts from Rochette and Catchen (Nat. Prot., 2017)[https://www.nature.com/articles/nprot.2017.123]
Paths would require some modfication to use these scripts.

2. The folder **Data** has some of the processed downstream VCF files (ddRADSeq, Wolbachia) and fasta files for the mtDNA used for most of the analysis in the above citation (see details in the text). It also contains barcodes for demultiplexing the 
associated illumina data archived in the NCBI SRA (SRR11238035 & SRR11238036). The mtDNA (fas) contains sequences for all individuals with >= 634 bp that was used for the haplotype network in Figure 1.

3. The folder **RScripts** contain several R scripts using the VCF files above for various analyses including those in the main text and supplementary info. There are also scripts for analyses not presented in the paper. Location of sources files will need to be changed in these scripts.

Please excuse the typos and if there any questions/concerns, please contact me at sarif@brookes.ac.uk


