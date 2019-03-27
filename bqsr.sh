#!/bin/bash

bam=$1
out=$2

# TODO:
# Indels for known sites resource, yay or nay?

basedir=/mnt/fargen/experiments/joint_call
gatk=$basedir/software/gatk-4.1.0.0/gatk

# Folder containing resources.
resources_dir=$basedir/resources

# Resources. Reference sequence and target regions.
reference=$resources_dir/reference_10x_genomics/refdata-GRCh38-2.1.0/fasta/genome.fa
targets=$resources_dir/sureselect_human_all_exon_v6_utr_grch38/S07604624_Padded.bed
dbsnp=$resources_dir/gatk_bundle/Homo_sapiens_assembly38.dbsnp138/Homo_sapiens_assembly38.dbsnp138.vcf
mills=$resources_dir/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

$gatk BaseRecalibrator \
    -I $bam \
    -R $reference \
    --known-sites $dbsnp \
    --known-sites $mills \
    -O recal_data.table \
    --tmp-dir=tmp \
    --java-options "-Xmx10g -Xms10g"

$gatk ApplyBQSR \
    -R $reference \
    -I $bam \
    --bqsr-recal-file recal_data.table \
    -L $targets \
    -O $out \
    --tmp-dir=tmp \
    --java-options "-Xmx10g -Xms10g"

$gatk AnalyzeCovariates \
    -bqsr recal_data.table \
    -plots AnalyzeCovariates.pdf
