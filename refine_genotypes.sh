#!/bin/bash

# TODO:
# Consider whether to use a supporting dataset. I commented out the "-supporting" argument,
# because it biases the data toward the population the supporting dataset is based on.
# Supply pedigree information (does this use more than just trios?)
# Consider which filters, if any, to apply.

vcf=$1
out=$2

basedir=/mnt/fargen/experiments/joint_call
gatk=$basedir/software/gatk-4.1.0.0/gatk

# Folder containing resources.
resources_dir=$basedir/resources

# Resources. Reference sequence and target regions.
kGphase3=$resources_dir/gatk_bundle/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf

# Adds fields PP and PG, and updates GQ and GT (genotype) fields.
$gatk CalculateGenotypePosteriors \
    -V $vcf \
    -O $out
#    -supporting $kGphase3

#$gatk VariantFiltration \
#    -R $reference \
#    -V posteriors.vcf \
#    -O $out \
#    --genotype-filter-name "GQ20" \
#    --genotype-filter-expression "GQ<=20"

