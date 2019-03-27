#!/bin/bash

# TODO:
# Is --skip-duplicated necessary?
# Is supplying targets necessary? (qualimap doesn't accept my BED file).

bam=$1
outdir=$2

basedir=/mnt/fargen/experiments/joint_call

targets=$basedir/resources/sureselect_human_all_exon_v6_utr_grch38/S07604624_Padded.bed

qualimap=$basedir/software/qualimap_v2.2.1/qualimap

$qualimap bamqc \
    -gd HUMAN \
    -bam $bam \
    -outdir $outdir \
    --skip-duplicated \
    --collect-overlap-pairs \
    -nt 1 \
    --java-mem-size=10G
#    -gff $targets \

