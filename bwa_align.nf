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
println "reference          : ${params.reference}"
println "targets            : ${params.targets}"
println "outdir             : ${params.outdir}"

// Get file handlers for input files.
reference = file(params.reference)
targets = file(params.targets)
outdir = file(params.outdir)

// Get FASTQ files, assuming multiple lanes and paired reads, and that the files are gzip compressed.
fastq_ch = Channel.fromFile(params.fastq_path + '*L*R*.gz')

// A single FASTQ file to get the sample name from and to construct the readgroup from.
Channel.fromPath(params.fastq_path + '/*L001*R1*.gz').into { fastq_sample_ch; fastq_rg_ch }

// Construct a readgroup from the filename of one of the input FASTQ files.
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

// Construct a readgroup from the sequence identifier in one of the input FASTQ files.
process get_readgroup {
    input:
    file fastq from fastq_rg_ch

    output:
    stdout readgroup_ch

    script:
    """
    get_readgroups.py $fastq
    """
}

// Create a new channel with (sample, readgroup, fastq path) tuples for each lane and read.
sample_rg_fastq_ch = sample_ch.combine( readgroup_ch.combine(fastq_ch) )

// Align the reads with BWA.
process align {
    input:
    set sample, rg, file(fastq) from sample_rg_fastq_ch

    output:
    set sample, file("aligned.bam") into aligned_merge_bam_ch

    script:
    """
    bwa mem -p -t ${task.cpus} -M -R '$rg' $reference $fastq | \
        samtools view -b -o aligned.bam
    """
}


// Merge BAMs. All BAMs have the same readgroup, so the RG and PG headers wil be combined.
process merge_bams {
    input:
    set sample, file(bams) from aligned_bam_merge_ch.collect()

    output:
    set sample, file("merged.bam") into merged_bam_sort_ch

    script:
    bam_list = (bams as List).join(' ')
    """
    samtools merge -@ ${task.cpus} -O bam -l 0 -c -p "merged.bam" $bam_list
    """
}


// Coordinate sort BAM.
process sort_bam {
    input:
    set sample, file(bam) from merged_bam_sort_ch

    output:
    set sample, file("sorted.bam") into sorted_bam_markdup_ch

    script:
    """
    samtools sort -@ ${task.cpus} -O bam -l 0 -m 4G -o "sorted.bam" $bam
    """
}

// Mark duplicates in BAM.
process mark_dup {
    input:
    set sample, file(bam) from sorted_bam_markdup_ch

    output:
    set sample, file("marked_dup.bam") into marked_bam_index_ch

    script:
    """
    gatk MarkDuplicates -I $bam -O "marked_dup.bam" -M "marked_dup_metrics.txt"
    """
}

// Index the BAM.
process index_bam {
    publishDir "$outdir/bam/aligned/bwa/$sample", mode: 'copy', pattern: '*.bam',
        saveAs: { filename -> "${sample}.bam" }
    publishDir "$outdir/bam/aligned/bwa/$sample", mode: 'copy', pattern: '*.bam.bai',
        saveAs: { filename -> "${sample}.bam.bai" }

    input:
    set sample, file(bam) from marked_bam_index_ch

    output:
    set sample, file("marked_dup.bam"), file("indexed.bam.bai") into indexed_bam_qc_ch

    script:
    """
    gatk BuildBamIndex -I $bam -O "indexed.bam.bai"
    """
}

// Run Qualimap for QC metrics of aligned and recalibrated BAM.
process qualimap_analysis {
    publishDir "$outdir/bam/aligned/bwa/$sample/qc", mode: 'copy'

    input:
    set sample, file(bam), file(bai) from indexed_bam_qc_ch

    output:
    file "qualimap_results" into qualimap_results_ch

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
        --java-mem-size="${task.memory.toGiga()}G"
    """
}









