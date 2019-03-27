#!/bin/bash

# TODO:
# There's a ton of information on the GATK website about HaplotypeCaller and its usage in variant discovery. Read it.
# Call variants present in dbsnp?

bam=$1
out=$2

basedir=/mnt/fargen/experiments/joint_call
gatk=$basedir/software/gatk-4.1.0.0/gatk

# Folder containing resources.
resources_dir=$basedir/resources

# Resources. Reference sequence and target regions.
reference=$resources_dir/reference_10x_genomics/refdata-GRCh38-2.1.0/fasta/genome.fa
targets=$resources_dir/sureselect_human_all_exon_v6_utr_grch38/S07604624_Padded.bed
dbsnp=$resources_dir/gatk_bundle/Homo_sapiens_assembly38.dbsnp138/Homo_sapiens_assembly38.dbsnp138.vcf

# NOTE:
# stand_call_conf is 0 by default in GVCF mode.

$gatk HaplotypeCaller  \
    -I $bam \
    -O $out \
    -R $reference \
    -L $targets \
    --dbsnp $dbsnp \
    -ERC GVCF \
    --create-output-variant-index \
    --annotation MappingQualityRankSumTest \
    --annotation QualByDepth \
    --annotation ReadPosRankSumTest \
    --annotation RMSMappingQuality \
    --annotation FisherStrand \
    --annotation Coverage \
    --verbosity INFO \
    --tmp-dir=tmp \
    --java-options "-Xmx10g -Xms10g"

