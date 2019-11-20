#!/usr/bin/env nextflow
/*
Author: Ólavur Mortensen <olavur@fargen.fo>
*/

// Input parameters.
params.fastq_path = null
params.reference = null
params.targets = null
params.outdir = null
params.help = false

// TODO: make help string
// Help message
helpMessage = """
Parameters:
--outdir            Desired path/name of folder to store output in.
""".stripIndent()

// Show help when needed
if (params.help){
    log.info helpMessage
        exit 0
}

// Make sure necessary input parameters are assigned.
assert params.fastq_path != null, 'Input parameter "fastq_path" cannot be unasigned.'
assert params.reference != null, 'Input parameter "reference" cannot be unasigned.'
assert params.targets != null, 'Input parameter "targets" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

println "P I P E L I N E     I P U T S    "
println "================================="
println "fastq_path         : ${params.fastq_path}"
println "targets            : ${params.targets}"
println "outdir             : ${params.outdir}"

// Get file handlers for input files.
reference = file(params.reference)
targets = file(params.targets)
outdir = file(params.outdir)
fastq_path = file(params.fastq_path)

// A single FASTQ file to get the sample name from and to construct the readgroup from.
fastq_sample_ch = Channel.fromPath(params.fastq_path + '/*L001*R1*.gz')

// Get the sample name from one of the FASTQ files.
process get_samplename {
    input:
    file fastq from fastq_sample_ch

    output:
    stdout sample_ch

    script:
    """
    get_samplenames.py $fastq
    """
}

// Align FASTQ reads to reference with LongRanger ALIGN command.
// https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/other-pipelines
process align_reads {
    input:
    val sample from sample_ch

    output:
    set sample, file("$sample/outs/possorted_bam.bam"), file("$sample/outs/possorted_bam.bam.bai") into aligned_bam_qc_ch

    script:
    """
    longranger align --id=$sample \
        --reference=$reference \
        --fastqs=$fastq_path \
        --localcores=${task.cpus} \
        --localmem=${task.memory}
    """
}

// Run Qualimap for QC metrics of aligned and recalibrated BAM.
process qualimap_analysis {
    publishDir "$outdir/bam/aligned/lr/$sample/qc", mode: 'copy'

    input:
    set sample, file(bam), file(bai) from aligned_bam_qc_ch

    output:
    set sample, file("qualimap_results") into qualimap_results_ch

    script:
    """
    # This first line adds two columns to our BED file with target regions, as QualiMap expects these.
    # The fifth and sixth column are respectively just "0" and ".", which has no information about the
    # regions.
    awk 'BEGIN{OFS="\\t"}{ if(NR > 2) { print \$1,\$2,\$3,\$4,0,"." } }' $targets > 'targets_6_fields.bed'

    # FIXME:
    # remove this when done testing.
    # and uncomment awk command above
    #echo 'track name="dummy" description="dummy BED" color=0,0,128 db=hg38' > 'targets_6_fields.bed'
    #echo 'chr17    1  83257441 allchr17   0   .' > 'targets_6_fields.bed'
    # FIXME

    # Make sure QualiMap doesn't attemt to open a display server.
    unset DISPLAY
    # Run QualiMap.
    qualimap bamqc \
        -gd HUMAN \
        -bam $bam \
        -gff 'targets_6_fields.bed' \
        -outdir "qualimap_results" \
        --skip-duplicated \
        --collect-overlap-pairs \
        -nt ${task.cpus} \
        --java-mem-size=${task.memory.toGiga()}G
    """
}

