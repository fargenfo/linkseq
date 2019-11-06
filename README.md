
# exolink -- GATK best-practices pipeline adapted to linked-reads [WIP]

This pipeline aligns [linked-reads from 10x Genomics](https://www.10xgenomics.com/linked-reads/) and calls variants with GAKT. Specifically, germline short variant discovery (SNPs and indels) is performed according to [GATK best-practices](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145).

The pipeline is written in [Nextflow](https://www.nextflow.io/) and contains six sub-pipelines:

* `ema_align.nf`: align reads with [EMA](https://github.com/arshajii/ema/)
* `lr_align.nf`: align reads with [longranger ALIGN](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/other-pipelines)
* `bwa_align.nf`: align reads with [BWA](http://bio-bwa.sourceforge.net/)
* `single_sample.nf `: recalibrate BAM and call variants with GATK
* `multi_sample.nf`: joint genotyping and variant annotation and fitering.
* `phase.nf`: phase variants with [HapCUT2](https://github.com/vibansal/HapCUT2)

## Workflow

* Basecall and demultiplex raw sequencing data with our [demuxlink](https://github.com/olavurmortensen/demuxlink) pipline
* Process all your samples individualy with first `ema_align.nf` and `single_sample.nf`
* Process your samples jointly with `multi_sample.nf`

