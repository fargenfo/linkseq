#!/bin/bash

module load fastqc-0.11.7

fastq_path=$1
out_dir=$2
tmp_dir=$3

adapters=/opt/fastqc-0.11.7/Configuration/adapter_list.txt
contaminants=/opt/fastqc-0.11.7/Configuration/contaminant_list.txt
limits=/opt/fastqc-0.11.7/Configuration/limits.txt

fastqs=`ls $fastq_path/*`
fastqc --outdir $out_dir --dir $tmp_dir -a $adapters --contaminants $contaminants -l $limits $fastqs
