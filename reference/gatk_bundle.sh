#!/bin/bash
#
# Download resources from GATK Resource Bundle. Links below can be found at the link below, although the available
# resources may vary from time to time. These resources were present on the 22-03-2019.
# https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0
#
# More information about the GATK Resource Bundle at https://software.broadinstitute.org/gatk/download/bundle

mkdir 1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38
cd 1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx
cd ..

mkdir 1000G_omni2.5.hg38
cd 1000G_omni2.5.hg38
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi
cd ..

mkdir 1000G_phase1.snps.high_confidence.hg38
cd 1000G_phase1.snps.high_confidence.hg38
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
cd ..

mkdir Axiom_Exome_Plus.genotypes.all_populations.poly.hg38
cd Axiom_Exome_Plus.genotypes.all_populations.poly.hg38
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi
cd ..

mkdir Homo_sapiens_assembly38.dbsnp138
cd Homo_sapiens_assembly38.dbsnp138
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
cd ..

mkdir Homo_sapiens_assembly38
cd Homo_sapiens_assembly38
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
cd ..

mkdir Homo_sapiens_assembly38.known_indels
cd Homo_sapiens_assembly38.known_indels
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
cd ..

mkdir Mills_and_1000G_gold_standard.indels.hg38
cd Mills_and_1000G_gold_standard.indels.hg38
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
cd ..

mkdir hapmap_3.3.hg38.vcf.gz
cd hapmap_3.3.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi
cd ..

mkdir wgs_calling_regions.hg38
cd wgs_calling_regions.hg38
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list
cd ..

