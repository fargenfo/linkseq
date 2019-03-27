#!/bin/bash

vcf=$1

basedir=/mnt/fargen/experiments/joint_call
gatk=$basedir/software/gatk-4.1.0.0/gatk

# Folder containing resources.
resources_dir=$basedir/resources

# Resources. Reference sequence and target regions.
reference=$resources_dir/reference_10x_genomics/refdata-GRCh38-2.1.0/fasta/genome.fa
dbsnp=$resources_dir/dbsnp_Homo_sapiens_assembly38/Homo_sapiens_assembly38.dbsnp138.vcf

# Validate variants sice we used a non-GAKT tool.
$gatk ValidateVariants \
    -V $vcf \
    -R $reference \
    --dbsnp $dbsnp

