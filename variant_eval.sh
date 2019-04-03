#!/bin/bash

# TODO:
# -comp and -kown-name can be used to specify a different set of known variants.
# If I want to know the novel concordant rate, I must supply another known variants file with -comp.
# The "Novelty" stratification module is default.
# Multithreading.
# Consider which evaluation and stratification modules to use. The "Filter" and "Sample" stratification modules seem like a good idea.

vcf=$1

basedir=/mnt/fargen/experiments/joint_call
gatk=$basedir/software/gatk-4.1.0.0/gatk

# Folder containing resources.
resources_dir=$basedir/resources

# Resources. Reference sequence and target regions.
reference=$resources_dir/reference_10x_genomics/refdata-GRCh38-2.1.0/fasta/genome.fa
dbsnp=$resources_dir/gatk_bundle/Homo_sapiens_assembly38.dbsnp138/Homo_sapiens_assembly38.dbsnp138.vcf
targets=$resources_dir/sureselect_human_all_exon_v6_utr_grch38/S07604624_Padded.bed

# VariantEval fails if the output file doesn't already exist. NOTE: this should be fixed in a newer version of GATK, as of the 19th of February 2019.
echo -n > variant_eval.table

$gatk VariantEval \
    -R $reference \
    --eval $vcf \
    --output variant_eval.table \
    --dbsnp $dbsnp \
    -L $targets \
    -no-ev \
    --eval-module TiTvVariantEvaluator \
    --eval-module CountVariants \
    --eval-module CompOverlap \
    --eval-module ValidationReport \
    --stratification-module Filter \


