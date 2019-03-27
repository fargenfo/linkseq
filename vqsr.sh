#!/bin/bash

# TODO:
# Increasing --max-gaussians may work for larger sample sizes. For two samples, --max-gaussians=4 failed. exoseq uses 4.
# Filtering with VQSR based on DP is not recommended for exome data. Don't know if I'm currently doing this.

vcf=$1
out=$2

basedir=/mnt/fargen/experiments/joint_call
gatk=$basedir/software/gatk-4.1.0.0/gatk
gatk3=$basedir/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar

# Folder containing resources.
resources_dir=$basedir/resources

# Resources. Reference sequence and target regions.
reference=$resources_dir/reference_10x_genomics/refdata-GRCh38-2.1.0/fasta/genome.fa
dbsnp=$resources_dir/gatk_bundle/Homo_sapiens_assembly38.dbsnp138/Homo_sapiens_assembly38.dbsnp138.vcf
mills=$resources_dir/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
kG=$resources_dir/gatk_bundle/1000G_phase1.snps.high_confidence.hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
omni=$resources_dir/gatk_bundle/1000G_omni2.5.hg38/1000G_omni2.5.hg38.vcf.gz
hapmap=$resources_dir/gatk_bundle/hapmap_3.3.hg38.vcf.gz/hapmap_3.3.hg38.vcf.gz

snps=snps.vcf
snps_recal=snps_recal.vcf
tranches=snps.tranches
recal=snps.recal
plots=snps.plots.R

$gatk SelectVariants \
    -R $reference \
    -V $vcf \
    -O $snps \
    --select-type-to-include SNP

$gatk VariantRecalibrator \
    -R $reference \
    -V $snps \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
    -resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $kG \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    --max-gaussians 4 \
    -O $recal \
    --tranches-file $tranches \
    --rscript-file $plots

$gatk ApplyVQSR \
    -R $reference \
    -V $snps \
    -O $snps_recal \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file $tranches \
    --recal-file $recal \
    -mode SNP

indels=indels.vcf
indels_recal=indels_recal.vcf
tranches=indels.tranches
recal=indels.recal
plots=indels.plots.R

$gatk SelectVariants \
    -R $reference \
    -V $vcf \
    -O $indels \
    --select-type-to-include INDEL \
    --select-type-to-include MIXED \
    --select-type-to-include MNP \
    --select-type-to-include SYMBOLIC \
    --select-type-to-include NO_VARIATION

$gatk VariantRecalibrator \
    -R $reference \
    -V $indels \
    -resource:mills,known=false,training=true,truth=true,prior=12.0 $mills \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode INDEL \
    --max-gaussians 4 \
    -O $recal \
    --tranches-file $tranches \
    --rscript-file $plots

$gatk ApplyVQSR \
    -R $reference \
    -V $indels \
    -O $indels_recal \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file $tranches \
    --recal-file $recal \
    -mode INDEL

# FIXME: Merge variants some other way. CombineVariants may not be included in GATK4 in the near future.
java -jar $gatk3 -T CombineVariants \
    -R $reference \
    -V:snps $snps_recal \
    -V:indels $indels_recal \
    -o recal.vcf \
    --genotypemergeoption PRIORITIZE \
    --rod_priority_list snps,indels

