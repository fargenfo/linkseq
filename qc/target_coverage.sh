#!/bin/bash

# NOTE: WIP

bam=$1
bed=$2
out=$3

bedtools coverage -a $bed -b $bam -mean > $out
