
# LinkSeq -- GATK best-practices pipeline adapted to linked-reads [WIP]

**Note:** This pipeline is a work in progress.

## Table of Contents  
* [Overview](https://github.com/olavurmortensen/linkseq#overview)
* [Workflow](https://github.com/olavurmortensen/linkseq#workflow)
* [Basecalling/demultiplexing and trimming with `demux`](https://github.com/olavurmortensen/linkseq#basecallingdemultiplexing-and-trimming-with-demux)
	* [Trimming](https://github.com/olavurmortensen/linkseq#trimming)
	* [Setup](https://github.com/olavurmortensen/linkseq#setup)
	* [Running on tiny-bcl](https://github.com/olavurmortensen/linkseq#running-on-tiny-bcl)
		* [Running demux pipeline](https://github.com/olavurmortensen/linkseq#run-demux-pipeline)
		* [Output](https://github.com/olavurmortensen/linkseq#output)
* [Align reads and call variants with `single_sample.nf` and `multi_sample.nf`](https://github.com/olavurmortensen/linkseq#align-reads-and-call-variants-with-single_samplenf-and-multi_samplenf)
* [Phase variants with `phase.nf`](https://github.com/olavurmortensen/linkseq#phase-variants-with-phasenf)
* [Reference resources](https://github.com/olavurmortensen/linkseq#reference-resources)
	* [Barcode whitelist](https://github.com/olavurmortensen/linkseq#barcode-whitelist)
	* [GATK resources](https://github.com/olavurmortensen/linkseq#gatk-resources)
	* [Exome sequencing targets](https://github.com/olavurmortensen/linkseq#exome-sequencing-targets)

## Overview

This pipeline basecalls/demultiplexes and aligns [linked-reads from 10x Genomics](https://www.10xgenomics.com/linked-reads/) and calls variants with GAKT. Specifically, germline short variant discovery (SNPs and indels) is performed according to [GATK best-practices](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145).

The pipeline is written in [Nextflow](https://www.nextflow.io/) and contains seven sub-pipelines:

* `demux.nf`: basecall/demultiplex raw BCL data with [bcl2fastq](https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf)
* `single_sample.nf `: align reads with [EMA](https://github.com/arshajii/ema/), recalibrate BAM and call variants with GATK
* `multi_sample.nf`: joint genotyping and variant annotation and fitering.
* Alignment only:
	* `ema_align.nf`: align reads with [EMA](https://github.com/arshajii/ema/)
	* `lr_align.nf`: align reads with [longranger ALIGN](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/other-pipelines)
	* `bwa_align.nf`: align reads with [BWA](http://bio-bwa.sourceforge.net/)
* `phase.nf`: phase variants with [HapCUT2](https://github.com/vibansal/HapCUT2)

## Workflow

* Basecall and demultiplex raw sequencing data with `demux.nf`
* Process all your samples individualy with `single_sample.nf`
* Process your samples jointly with `multi_sample.nf`

## Basecalling/demultiplexing and trimming with `demux`

This Nextflow pipeline basecalls and demultiplexes linked-reads from 10x Genomics. To run this pipeline, the 8-base sample indexes are needed, corresponding to the 10x Genomics indexes (e.g. `SI-GA-A1`).

This pipeline makes some assumptions about the input data. For example, it makes the assumption that it is paired-end sequencing, and therefore uses `--use-bases-mask=Y*,I*,Y*` in `bcl2fastq`, and assumes that the read lengths (and index length) is found in `RunInfo.xml`.

### Trimming

There are four trimming steps in the `demux.nf` pipeline, each of which is listed below.

* `trim_adapters` trims adapter sequences from 3' end using [BBtools/BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/).
* `bctrim` trims linked-read barcodes from read 2. The `trimR2bc.py` is written by Elisabet Thomsen and can be found in the [bin folder](https://github.com/olavurmortensen/linkseq/blob/master/bin/trimR2bc.py) of this project. When the insert size is small and read 1 and 2 overlap, read 2 may be contaminated by the barcode attached to read 1.
* `polyG_trim` trims poly-G tails from reads using [fastp](https://github.com/OpenGene/fastp).
* `quality_trim_read1`/`quality_trim_read2` trims low quality bases from read 1 and 2 respectively.

### Setup

Install dependencies with `conda` using the [conda_envs/demux.yml file](https://github.com/olavurmortensen/linkseq/blob/master/conda_envs/demux.yml):

```
conda env create -f demux.yml
```

Activate the environment (check the name of the environment, it should be `linkseq-demux`):

```
conda activate linkseq-demux
```

Pull this project with `nextflow`:

```
nextflow pull https://github.com/olavurmortensen/linkseq
```

### Running on tiny-bcl

Here's how to run this pipeline on the "tiny-bcl" example dataset from 10x Genomics. First of all, download the tiny-bcl tar file and the Illumina Experiment Manager sample sheet: tiny-bcl-samplesheet-2.1.0.csv.

> Download the tiny-bcl data:
> https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/mkfastq#example_data

Next, edit the `[Data]` part of the samplesheet from the following:

```
Lane,Sample_ID,index,Sample_Project
5,Sample1,SI-GA-C5,tiny_bcl
```

To the following:

```
Lane,Sample_ID,index
5,Sample1,CGACTTGA
5,Sample1,TACAGACT
5,Sample1,ATTGCGTG
5,Sample1,GCGTACAC
```

Using that the index `SI-GA-C5` corresponds to the four octamers `CGACTTGA,TACAGACT,ATTGCGTG,GCGTACAC` (10X Genomics have a tool for this on [their website](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/bcl2fastq-direct)). Notice that we also removed the `Sample_Project` column.

Also add adapter sequences to the `[Settings]` part of the samplesheet:
```
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```

#### Run demux pipeline

Use the `example_configs/demux.config` to get an idea how to define input parameters. For the `tiny-bcl` dataset, you should set `rundir` to `tiny-bcl-2.0.0` and `samplesheet` to `tiny-bcl-samplesheet-2.1.0.csv`. Additionally, you need the barcode whitelist, see the [barcode whitelist](https://github.com/olavurmortensen/linkseq#barcode-whitelist) section on how to obtain this list.

Change the memory and CPU specifications in the configuration (under `process` and `executor`) to suit your needs before continuing.

When you've made the config file and activated the `linkseq-demux` environment, run the pipeline like this:

```
nextflow run olavurmortensen/linkseq/demux.nf -c [your config]
```

#### Output

The output from the pipeline, run on the `tiny-bcl` data, is shown below. The compressed FASTQ data is in `outs/Sample1/fastqs`. The pipeline runs `FastQC` for quality control and the reports are in `outs/Sample1/fastqc`. There are various logs, from the basecalling itself via `Bcl2Fastq`, from the various trimming steps, from read synchronization, and from `FastQC`.

```
$ tree outs/
outs/
├── bcl2fastq.log
└── Sample1
    ├── fastqc
    │   ├── fastqc.log
    │   ├── Sample1_L005_R1_fastqc.html
    │   ├── Sample1_L005_R2_fastqc.html
    │   └── zips
    │       ├── Sample1_L005_R1_fastqc.zip
    │       └── Sample1_L005_R2_fastqc.zip
    ├── fastqs
    │   ├── Sample1_L005_R1.fastq.gz
    │   └── Sample1_L005_R2.fastq.gz
    └── logs
        ├── adapter_trim
        │   └── L005.log
        ├── bctrim
        │   └── L005.log
        ├── polyG_trim
        │   └── L005.log
        ├── quality_trim
        │   ├── L005_R1.log
        │   └── L005_R2.log
        └── sync_reads
            └── L005.log
```


## Align reads and call variants with `single_sample.nf` and `multi_sample.nf`

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


