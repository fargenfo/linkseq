#!/usr/bin/env nextflow
/*
Author: Ólavur Mortensen <olavur@fargen.fo>
*/

// Input parameters.
params.fastq_paths = null
params.reference = null
params.targets = null
params.threads = null
params.mem = null
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
assert params.fastq_paths != null, 'Input parameter "fastq_paths" cannot be unasigned.'
assert params.reference != null, 'Input parameter "reference" cannot be unasigned.'
assert params.targets != null, 'Input parameter "targets" cannot be unasigned.'
assert params.threads != null, 'Input parameter "threads" cannot be unasigned.'
assert params.mem != null, 'Input parameter "mem" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

println "P I P E L I N E     I P U T S    "
println "================================="
println "fastq_paths         : ${params.fastq_paths}"
println "targets            : ${params.targets}"
println "threads            : ${params.threads}"
println "mem                : ${params.mem}"
println "outdir             : ${params.outdir}"

// TODO: get rid of all this "reference_fa" business. At some point I thought that for longranger I couldn't
// supply the path directly to the FASTA, but had to supply the path to the "reference folder". I think longranger
// can do both.

// Get file handlers for input files.
reference = file(params.reference)  // Directory of 10x reference.
reference_fa = file(params.reference + '/fasta/genome.fa')  // Reference fasta file.


// Turn the file with FASTQ paths into a channel with [sample, path] tuples.
fastq_paths_ch = Channel.fromPath(params.fastq_paths)
fastq_paths_ch
    .splitCsv(header: true)
    .map { it -> tuple(it.sample, it.fastq_path) }
    .into { fastq_align_ch; fastq_print_ch; fastq_qc_ch }

println("Processing data:\nSample\tFASTQ path")
fastq_print_ch.subscribe { println(it[0] + "\t" + it[1]) }

// Align FASTQ reads to reference with LongRanger ALIGN command.
// https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/other-pipelines
process align_reads {
    memory = "${params.mem}GB"
    cpus = params.threads

    input:
    set sample, fastq_path from fastq_align_ch

    output:
    set sample, file("$sample/outs/possorted_bam.bam"), file("$sample/outs/possorted_bam.bam.bai") into aligned_bam_prepare_ch, aligned_bam_apply_ch

    script:
    """
    longranger align --id=$sample \
        --reference=$reference \
        --fastqs=$fastq_path \
        --localcores=${params.threads} \
        --localmem=${params.mem}
    """
}

// Run Qualimap for QC metrics of aligned and recalibrated BAM.
process qualimap_analysis {
    memory = "${params.mem}GB"
    cpus = params.threads

    publishDir "${params.outdir}/bamqc", mode: 'copy',
        saveAs: {filename -> "$sample"}

    input:
    set sample, file(bam), file(bai) from recalibrated_bam_qualimap_ch

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
        -nt ${params.threads} \
        --java-mem-size=${params.mem}G
    """
}

