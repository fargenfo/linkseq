# LinkSeq -- GATK best-practices pipeline adapted to linked-reads

[![Docker build](https://img.shields.io/badge/Docker%20build-Available-informational)](https://hub.docker.com/repository/docker/olavurmortensen/linkseq)

## Table of Contents  
* [Overview](#overview)
* [Workflow](#workflow)
* [Setup](#setup)
* [Usage example](#usage-example)
	* [Run LinkSeq](#run-linkseq)
	* [A note on debugging](#a-note-on-debugging)
* [Reference resources](#reference-resources)
	* [Barcode whitelist](#barcode-whitelist)
	* [GATK resources](#gatk-resources)
	* [Exome sequencing targets](#exome-sequencing-targets)
	* [SnpEff data](#snpeff-data)

## Overview

This pipeline aligns [linked-reads from 10x Genomics](https://www.10xgenomics.com/linked-reads/) and calls variants with GAKT. Specifically, germline short variant discovery (SNPs and indels) is performed according to [GATK best-practices](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145).

The pipeline is written in [Nextflow](https://www.nextflow.io/). The main steps of the pipeline are summarized below.

* Align reads with [EMA](https://github.com/arshajii/ema/), and recalibrate BAM (BQSR)
* Call variants with GATK's `HaplotypeCaller`, yielding a GVCF (which can be used in joint genotyping)
* Genotype GVCF with GATK's `GenotypeGVCFs`, yielding a single-sample VCF
* Annotate variant effect with `SnpEff` and filter variants
* Phase VCF with [HapCUT2](https://github.com/vibansal/HapCUT2)
* Attach phasing from VCF to BAM using [WhatsHap](https://whatshap.readthedocs.io/en/latest/index.html)
* QC of variants with GATK's `VariantEval`
* QC of BAM with [Qualimap](http://qualimap.bioinfo.cipf.es/)
* QC report using [MultiQC](multiqc.info/)

<img src="https://github.com/olavurmortensen/linkseq/blob/master/linkseq_overview.png" width=500>

## Workflow

LinkSeq has a few sister pipelines. This section describes how they fit together.

* [LinkSeq Demux](https://github.com/olavurmortensen/linkseq-demux): Basecall and demultiplex raw sequencing data and trim reads.
* LinkSeq: For each sample, align reads, call variants, and phase VCF and BAM.
* [olavurmortensen/gatk-joint-genotyping](https://github.com/olavurmortensen/gatk-joint-genotyping): Perform joint genotyping of all samples.
* [LinkSeq Phase [WIP]](https://github.com/olavurmortensen/linkseq-phase): Phase multi-sample joint genotyped VCF (work in progress).

## Setup

You have two options: (1) use conda to install dependencies, as described below, or (2) use the [olavurmortensen/linkseq](https://hub.docker.com/repository/docker/olavurmortensen/linkseq) Docker image.

Install dependencies with `conda` using the `environment.yml` file:

```
conda env create -f environment.yml
```

Activate the environment (check the name of the environment in the `environment.yml` file, it should be `linkseq`):

```
conda activate linkseq
```

Pull this project with `nextflow`:

```
nextflow pull https://github.com/olavurmortensen/linkseq
```

## Usage example

We will run `linkseq` on the "tiny-bcl" example dataset from 10x Genomics. Before running this example, run the example from the `linkseq-demux` (https://github.com/olavurmortensen/linkseq-demux) pipeline to basecall, demultiplex and trim the raw sequences.

For information about reference data used in pipeline, see the *Reference resources* section below.

### Run LinkSeq

The easiest way to run the pipeline is by defining a configuration file for nextflow. The configuration file below defines the input files to the pipeline as well as some runtime settings.

```
params {
    sample = "Sample1"
    
    // Read 1 and 2 FASTQ files.
    fastq_r1 = "tiny-fastq/Sample1/fastqs/Sample1_L005_R1*fastq.gz"
    fastq_r2 = "tiny-fastq/Sample1/fastqs/Sample1_L005_R2*fastq.gz"
    
    // Reference data.
    whitelist = "resources/whitelist/4M-with-alts-february-2016.txt"
    reference = "resources/gatk_bundle/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta"
    dbsnp = "resources/gatk_bundle/Homo_sapiens_assembly38.dbsnp138/Homo_sapiens_assembly38.dbsnp138.vcf"
    targets = "resources/sureselect_human_all_exon_v6_utr_grch38/S07604624_Padded_modified.bed"
    snpeff_datadir = "resources/snpeff_data"
    
    // How many bins to use in EMA alignment.
    bcbins = 10
    
    outdir = "outs"
}

// Resources available to each process.
process {
    executor = 'local'
    memory = '10GB'
    cpus = 1 
}

// Total resources avilable to pipeline.
executor {
    name = 'local'
    cpus = 10
    memory = '100GB'
    queueSize = 100 
}

// Capture exit codes from upstream processes when piping.
process.shell = ['/bin/bash', '-euo', 'pipefail']
```

Note that the number of "bins" in EMA alignment is set to 10 only because this is a tiny example dataset. For whole-genome data, for example, the EMA developers recommend 500 bins.

We can run the pipeline with the following command. We use `-with-trace` to get a log file with process progress. If the pipeline fails and we re-run it, `-resume` means it will continue from where it left off.

```bash
nextflow run olavurmortensen/linkseq -resume -with-trace
```

When the pipeline has completed, we can run `tree -L 3 outs/` to get an overview of the outputs, which we can see below. There is one folder for each sample, with aligned reads (`bam`), one with variants (`vcf`), and with GVCFs for joint genotyping (`gvcf`). The `multiqc_logs` folder contains various QC reports, these reports are primarily to create a MultiQC report (combined QC report).

```bash
outs/
├── Sample1
│   ├── bam
│   │   ├── Sample1.bam
│   │   ├── Sample1.bam.bai
│   │   └── bx_stats.csv
│   ├── gvcf
│   │   ├── gvcf.g.vcf
│   │   └── gvcf.g.vcf.idx
│   └── vcf
│       ├── Sample1.vcf.gz
│       ├── Sample1.vcf.gz.tbi
│       └── phasing
├── multiqc
│   ├── multiqc_data
│   │   ├── multiqc.log
│   │   ├── multiqc_data.json
│   │   ├── multiqc_fastqc.txt
│   │   ├── multiqc_gatk_varianteval.txt
│   │   ├── multiqc_general_stats.txt
│   │   ├── multiqc_qualimap_bamqc_genome_results.txt
│   │   ├── multiqc_snpeff.txt
│   │   └── multiqc_sources.txt
│   └── multiqc_report.html
└── multiqc_logs
    ├── AnalyzeCovariates
    │   └── Sample1
    ├── SnpEff
    │   └── Sample1
    ├── VariantEval
    │   └── Sample1
    ├── WhatsHap
    │   └── Sample1
    ├── bqsr_after
    │   └── Sample1
    ├── bqsr_before
    │   └── Sample1
    ├── bx_stats
    ├── fastqc
    │   └── Sample1_fastqc.zip
    └── qualimap
        └── Sample1
```

### A note on debugging

If you need to debug the pipeline, because it failed, it can be useful to inspect the pipeline "trace":

```bash
less trace.txt
```

And the files in the `work/` directory:

```bash
tree work/ | less
```

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

### Interval BED file

You should always use an interval BED file when running this pipeline. If you're running whole-genome sequencing, you might want to use WGS calling regions such as those from the GATK Resource Bundle:

> GATK Resource Bundle
>
> https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle

Specifically this file:

> WGS calling regions hg38
>
> https://storage.cloud.google.com/genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list

Our particular sequencing experiment uses the Agilent SureSelect Human All Exon V6 UTR kit to capture the exome. The folder `reference/sureselect_human_all_exon_v6_utr_grch38` contains an **example BED file**, and a `README` that explains where this file comes from. If `LinkSeq` is unwilling to accept your BED file, this example may help debug the problem.

### SnpEff data

Use the `reference/snpeff_data.sh` script to download SnpEff database for `hg38`. In `LinkSeq`, the path to this database is used together with the `dataDir` argument in SnpEff.

The handle of this database can also be found by running `snpEff databases | grep hg38`.

This is particularly useful when running the pipeline in the Docker container, as it saves a lot of time. If we do not supply `-dataDir` in SnpEff, it will download the database every time it is run. This download may even sporadically fail, crashing the whole pipeline.
