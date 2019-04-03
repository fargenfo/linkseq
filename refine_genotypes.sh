#!/bin/bash

# TODO:
# Consider whether the chosen supporting dataset is suitable.
# Supply pedigree information (does this use more than just trios?)
# Consider which filters to apply.

vcf=$1
out=$2

basedir=/mnt/fargen/experiments/joint_call
gatk=$basedir/software/gatk-4.1.0.0/gatk

# Folder containing resources.
resources_dir=$basedir/resources

# Resources. Reference sequence and target regions.
reference=$resources_dir/reference_10x_genomics/refdata-GRCh38-2.1.0/fasta/genome.fa
dbsnp=$resources_dir/gatk_bundle/Homo_sapiens_assembly38.dbsnp138/Homo_sapiens_assembly38.dbsnp138.vcf
targets=$resources_dir/sureselect_human_all_exon_v6_utr_grch38/S07604624_Padded.bed

# Adds fields PP and PG, and updates GQ and GT (genotype) fields.
$gatk CalculateGenotypePosteriors \
    -V $vcf \
    -O posteriors.vcf \
    -supporting $kGphase3

$gatk VariantFiltration \
    -R $reference \
    -V posteriors.vcf \
    -O $out \
    --genotype-filter-name "GQ20" \
    --genotype-filter-expression "GQ<=20"

