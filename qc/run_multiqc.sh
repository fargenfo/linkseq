#!/bin/bash

# TODO:
# Results from VariantEval can be supplied here, if they are available.

fastqc_path=$1
bamqc_path=$2
recal_path=$3
snpeff_path=$4
out=$?

multiqc $fastqc_path $bamqc_path -o $out

