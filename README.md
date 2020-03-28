
# LinkSeq -- GATK best-practices pipeline adapted to linked-reads [WIP]

**Note:** This pipeline is a work in progress.

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

**TODO:** an explaination of the trimming steps.

### Setup

Pull this project with `nextflow`:

```
nextflow pull https://github.com/olavurmortensen/linkseq
```

Install dependencies with `conda` using the [conda_envs/demux.yml file](https://github.com/olavurmortensen/linkseq/blob/master/conda_envs/demux.yml):

```
conda env create -f demux.yml
```

Activate the environment (check the name of the environment, it should be `linkseq-demux`):

```
conda activate linkseq-demux
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

#### Run demux pipeline

Use the `example_config/demux.config` to get an idea how to define input parameters. For the `tiny-bcl` dataset, you should set `rundir` to `tiny-bcl-2.0.0` and `samplesheet` to `tiny-bcl-samplesheet-2.1.0.csv`. Additionally, you need the barcode whitelist, which you can get in the [reference/whitelist subfolder of this project](https://github.com/olavurmortensen/linkseq/tree/master/reference/whitelist).

When you've made the config file and activated the `linkseq-demux` environment, run the pipeline like this:

```
nextflow run olavurmortensen/linkseq/demux.nf -c [your config]
```

#### Output



## Reference resources

To run this pipeline you need some reference databases from the [GATK resource bundle](https://software.broadinstitute.org/gatk/download/bundle), as well as an exome targets file.

The `reference/gatk_bundle.sh` script downloads all the resources needed from the GATK resource bundle; note that if the GATK dev team change any of these resources, this script may fail. We downloaded the resources from the site below on the 22/03-2019:

> GATK resource bundle
> https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0

Our particular sequencing experiment uses the Agilent SureSelect Human All Exon V6 UTR kit to capture the exome. `reference/sureselect_human_all_exon_v6_utr_grch38` contains the target BED file we use and some details are in the README.

## Align and call variants

**TODO**

