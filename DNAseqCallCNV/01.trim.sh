#!/bin/bash

outdir=$1 # the output direcotry for trimmed reads
fq1=$2 # the r1
fq2=$3 # the r2

trim_galore --fastqc --paired --quality 20 --phred33 --stringency 3 --gzip --length 36  --output_dir $outdir $fq1 $fq2
