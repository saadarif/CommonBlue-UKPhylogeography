#!/bin/bash

#stored processed reads
clean_dir=/media/data_disk/PROJECTS/Saad/CommonBlue/cleanData
process=process_radtags
#barcode files
bc_plate1=/media/data_disk/PROJECTS/Saad/CommonBlue/INFO/barcodes_plate1
bc_plate2=/media/data_disk/PROJECTS/Saad/CommonBlue/INFO/barcodes_plate2
#plate fastqfiles
plate1fq=/media/data_disk/PROJECTS/Saad/CommonBlue/rawdata/Picarus20180603.fastq.gz
plate2fq=/media/data_disk/PROJECTS/Saad/CommonBlue/rawdata/Picarus20190702.fastq.gz

#directory for runlogs
runlogs=/media/data_disk/PROJECTS/Saad/CommonBlue/runlogs

#Clean data
#plate1 barcodes are truncated to remove past1 cutsite
$process -f $plate1fq --inline_null -c -q -r \
	-e pstI  -b $bc_plate1  -D  --retain_header\
	-o $clean_dir &> $runlogs/process_radtags.plate1.oe

#plate2 uses barcode and cutsite as the barcode
$process -f $plate2fq --inline_null -c -q -r \
          -b $bc_plate2  -D  --retain_header\
        -o $clean_dir &> $runlogs/process_radtags.plate2.oe

