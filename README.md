
# LinkSeq -- GATK best-practices pipeline adapted to linked-reads [WIP]

**Note:** This pipeline is a work in progress.

## Table of Contents  
* [Overview](#overview)
* [Workflow](#workflow)
* [Align reads and call variants with `single_sample.nf`](#align-reads-and-call-variants-with-single_samplenf)
* [Joint genotyping with `multi_sample.nf`](#joint-genotyping-with-multi_samplenf)
* [Phase variants with `phase.nf`](#phase-variants-with-phasenf)
* [Reference resources](#reference-resources)
	* [Barcode whitelist](#barcode-whitelist)
	* [GATK resources](#gatk-resources)
	* [Exome sequencing targets](#exome-sequencing-targets)

## Overview

This pipeline aligns [linked-reads from 10x Genomics](https://www.10xgenomics.com/linked-reads/) and calls variants with GAKT. Specifically, germline short variant discovery (SNPs and indels) is performed according to [GATK best-practices](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145).

The pipeline is written in [Nextflow](https://www.nextflow.io/) and contains two sub-pipelines. The steps of these pipelines are summarized below.

* `main.nf `:
	* Align reads with [EMA](https://github.com/arshajii/ema/), and recalibrate BAM (BQSR)
	* Call variants with GATK's `HaplotypeCaller`, yielding a GVCF (which can be used in `joint_genotyping.nf`)
	* Genotype GVCF with GATK's `GenotypeGVCFs`, yielding a single-sample VCF
	* Annotate variant effect with `SnpEff` and filter variants
	* Phase VCF with [HapCUT2](https://github.com/vibansal/HapCUT2)
	* Attach phasing from VCF to BAM using [WhatsHap](https://whatshap.readthedocs.io/en/latest/index.html)
	* QC of variants with GATK's `VariantEval`
	* QC of BAM with [Qualimap](http://qualimap.bioinfo.cipf.es/)
	* QC report using [MultiQC](multiqc.info/)
* `joint_genotyping.nf`: joint genotyping and variant annotation and fitering.
	* Joint genotyping of many GVCFs from `main.nf` using GATK's `GenotypeGVCFs`
	* Variant filtering using VQSR from GATK
	* Refine genotypes using GATK's `CalculateGenotypePosteriors`
	* Annotate variant effect with `SnpEff`
	* QC of variants with GATK's `VariantEval`
	* QC report using MultiQC
	* **NOTE:** as of yet, no phasing in this pipeline

## Workflow

* Basecall and demultiplex raw sequencing data and trim reads with `linkseq-demux` (https://github.com/olavurmortensen/linkseq-demux).
* For each sample, align reads, call variants, and phase VCF and BAM with `main.nf`.
* Perform joint genotyping of all samples with `joint_genotyping.nf` (no phasing).

## Align reads and call variants with `single_sample.nf`

**TODO**

## Joint genotyping with `multi_sample.nf`

**TODO**

## Phase variants with `phase.nf`

**TODO**

## Reference resources

### Barcode whitelist

The whitelist contains barcodes in the 10X Genomics GemCode technology, and can be obtained in the LongRanger software bundle here:

longranger-2.y.z/longranger-cs/2.y.z/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt

LongRanger can be obtained on the 10X Genomics website:

https://support.10xgenomics.com/genome-exome/software/overview/welcome

It's easier to download the file from this link though:
```
http://cb.csail.mit.edu/cb/ema/data/4M-with-alts-february-2016.txt
```

### GATK resources

To run this pipeline you need some reference databases from the [GATK resource bundle](https://software.broadinstitute.org/gatk/download/bundle), as well as an exome targets file.

The `reference/gatk_bundle.sh` script downloads all the resources needed from the GATK resource bundle; note that if the GATK dev team change any of these resources, this script may fail. We downloaded the resources from the site below on the 22/03-2019:

> GATK resource bundle
> https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0

### Exome sequencing targets

Our particular sequencing experiment uses the Agilent SureSelect Human All Exon V6 UTR kit to capture the exome. `reference/sureselect_human_all_exon_v6_utr_grch38` contains the target BED file we use and some details are in the README.


