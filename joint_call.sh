#!/bin/bash

# TODO:
# Figure out some way of passing a list of VCFs to GenomicsDBImport.

n_threads=$1

basedir=/mnt/fargen/experiments/joint_call
gatk=$basedir/software/gatk-4.1.0.0/gatk

# Folder containing resources.
resources_dir=$basedir/resources

# Resources. Reference sequence, target regions, and a filter for CNV calls.
reference=$resources_dir/reference_10x_genomics/refdata-GRCh38-2.1.0/fasta/genome.fa
targets=$resources_dir/sureselect_human_all_exon_v6_utr_grch38/S07604624_Padded.bed
dbsnp=$resources_dir/gatk_bundle/Homo_sapiens_assembly38.dbsnp138/Homo_sapiens_assembly38.dbsnp138.vcf

export TILEDB_DISABLE_FILE_LOCKING=1
$gatk GenomicsDBImport \
    -V $basedir/data/gvcf/FN000009.g.vcf.gz \
    -V $basedir/data/gvcf/FN000010.g.vcf.gz \
    -V $basedir/data/gvcf/FN000011.g.vcf.gz \
    -L $targets \
    --genomicsdb-workspace-path data/genomicsdb/run2 \
    --merge-input-intervals \
    --tmp-dir=tmp \
    --java-options "-Xmx50g -Xms50g"

$gatk GenotypeGVCFs \
    -V gendb://data/genomicsdb/run1 \
    -R $reference \
    -O variants.vcf \
    --tmp-dir=tmp \
    --java-options "-Xmx50g -Xms50g"
