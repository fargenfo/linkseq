#!/bin/bash

vcf=$1
out=$2

# TODO:
# Maybe add more annotations in VariantAnnotator.

module load tabix-0.2.5

basedir=/mnt/fargen/experiments/joint_call
snpEff=$basedir/software/snpEff/snpEff.jar
snpSift=$basedir/software/snpEff/SnpSift.jar
gatk=$basedir/software/gatk-4.1.0.0/gatk

# Folder containing resources.
resources_dir=$basedir/resources

# Resources. Reference sequence and target regions.
reference=$resources_dir/reference_10x_genomics/refdata-GRCh38-2.1.0/fasta/genome.fa
dbsnp=$resources_dir/gatk_bundle/Homo_sapiens_assembly38.dbsnp138/Homo_sapiens_assembly38.dbsnp138.vcf


java -jar $snpEff \
    -i vcf \
    -o vcf \
    -csvStats snpEff_stats.csv \
    hg38 \
    -v \
    $vcf > snpeff.vcf

# Validate variants sice we used a non-GAKT tool.
$gatk ValidateVariants \
    -V snpeff.vcf \
    -R $reference \
    --dbsnp $dbsnp

# VariantAnnotator is still in beta (as of 20th of March 2019).
$gatk VariantAnnotator \
    -R $reference \
    -V snpeff.vcf \
    -O $out


