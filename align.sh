#!/bin/bash

module load longranger-2.2.2

sample_id=$1
fastq_path=$2
threads=$3

basedir=/mnt/fargen/experiments/joint_call

# Reference sequence.
reference=$basedir/resources/reference_10x_genomics/refdata-GRCh38-2.1.0

longranger align --id=$sample_id \
   --reference=$reference \
   --fastqs=$fastq_path \
   --localcores=$threads \
   --uiport=3003

