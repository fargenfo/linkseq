#!/usr/bin/env nextflow
/*
Author: Ólavur Mortensen <olavur@fargen.fo>
*/

// Input parameters.
params.bam_paths = null
params.reference = null
params.dbsnp = null
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
assert params.bam_paths != null, 'Input parameter "bam_paths" cannot be unasigned.'
assert params.reference != null, 'Input parameter "reference" cannot be unasigned.'
assert params.dbsnp != null, 'Input parameter "dbsnp" cannot be unasigned.'
assert params.targets != null, 'Input parameter "targets" cannot be unasigned.'
assert params.threads != null, 'Input parameter "threads" cannot be unasigned.'
assert params.mem != null, 'Input parameter "mem" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

println "P I P E L I N E     I P U T S    "
println "================================="
println "bam_paths          : ${params.bam_paths}"
println "reference          : ${params.reference}"
println "dbsnp              : ${params.dbsnp}"
println "targets            : ${params.targets}"
println "threads            : ${params.threads}"
println "mem                : ${params.mem}"
println "outdir             : ${params.outdir}"

// Get file handlers for input files.
bam_paths = file(params.bam_paths)
reference = file(params.reference)
dbsnp = file(params.dbsnp)
targets = file(params.targets)


// FIXME
// remove when done testing
//targets = "chr17"
// FIXME


// FIXME:
// This has not been tested.
//
// Turn the file with FASTQ paths into a channel with [sample, path] tuples.
bam_paths_ch = Channel.fromPath(params.bam_paths)
bam_paths_ch
    .splitCsv(header: true)
    .map { it -> tuple(it.sample, it.bam_path, it.bai_path) }
    .into { aligned_bam_prepare_ch; aligned_bam_apply_ch }

println("Processing data:\nSample\tBAM path\tBAI path")
fastq_print_ch.subscribe { println(it[0] + "\t" + it[1] + "\t" + it[2]) }


// FIXME
// remove when done testing
//process make_small_bam {
//    input:
//    set sample, file(bam), file(bai) from aligned_bam_prepare_ch
//
//    output:
//    set sample, file("small.bam"), file("small.bam.bai") into small_bam_ch
//
//    script:
//    """
//    samtools view -b -o "small.bam" -T $reference_fa $bam "chr17"
//    samtools index -b "small.bam"
//    """
//}
//
//aligned_bam_prepare_ch = null
//aligned_bam_apply_ch = null
//small_bam_ch.into { aligned_bam_prepare_ch; aligned_bam_apply_ch }

/*
The next three processes, prepare_bqsr_table, analyze_covariates, and apply_bqsr, deal with base quality score
recalibration, in preparation for GATK best practices.
BQSR: https://software.broadinstitute.org/gatk/documentation/article?id=44
*/

// Generate recalibration table for BQSR.
process prepare_bqsr_table {
    memory = "${params.mem}GB"
    cpus = params.threads

    input:
    set sample, file(bam), file(bai) from aligned_bam_prepare_ch

    output:
    set sample, file('bqsr.table') into bqsr_table_ch, bqsr_table_copy_ch

    script:
    """
    mkdir tmp
    gatk BaseRecalibrator \
            -I $bam \
            -R $reference_fa \
            --known-sites $dbsnp \
            -O 'bqsr.table' \
            --tmp-dir=tmp \
            --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    """
}

// Evaluate BQSR.
process analyze_covariates {
    memory = "${params.mem}GB"
    cpus = params.threads

    publishDir "${params.outdir}/bam/analyze_covariates", mode: 'copy', overwrite: true, saveAs: { filename -> "${sample}_$filename" }

    input:
    set sample, file(bqsr_table) from bqsr_table_ch

    output:
    set sample, file('AnalyzeCovariates.pdf') into bqsr_analysis_ch

    script:
    """
    gatk AnalyzeCovariates \
        -bqsr $bqsr_table \
        -plots 'AnalyzeCovariates.pdf'
    """
}

// Pair the BAM and the BQSR table in one "set" channel.
data_apply_bqsr_ch = aligned_bam_apply_ch.join(bqsr_table_copy_ch)

// Apply recalibration to BAM file.
process apply_bqsr {
    memory = "${params.mem}GB"
    cpus = params.threads

    publishDir "${params.outdir}/bam", mode: 'copy', overwrite: true

    input:
    set sample, file(bam), file(bai), file(bqsr_table) from data_apply_bqsr_ch

    output:
    set sample, file("${sample}.bam"), file("${sample}.bam.bai") into recalibrated_bam_call_ch, recalibrated_bam_qualimap_ch

    script:
    """
    mkdir tmp
    gatk ApplyBQSR \
        -R $reference_fa \
        -I $bam \
        --bqsr-recal-file $bqsr_table \
        -L $targets \
        -O "${sample}.bam" \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    mv "${sample}.bai" "${sample}.bam.bai"
    """
}

// Call variants in sample with HapltypeCaller, yielding a GVCF.
process call_sample {
    memory = "${params.mem}GB"
    cpus = params.threads

    publishDir "${params.outdir}/gvcf", mode: 'copy', overwrite: true

    input:
    set sample, file(bam), file(bai) from recalibrated_bam_call_ch

    output:
    set sample, file("${sample}.g.vcf"), file("${sample}.g.vcf.idx") into gvcf_ch

    script:
    """
    mkdir tmp
    gatk HaplotypeCaller  \
        -I $bam \
        -O "${sample}.g.vcf" \
        -R $reference_fa \
        -L $targets \
        --dbsnp $dbsnp \
        -ERC GVCF \
        --create-output-variant-index \
        --annotation MappingQualityRankSumTest \
        --annotation QualByDepth \
        --annotation ReadPosRankSumTest \
        --annotation RMSMappingQuality \
        --annotation FisherStrand \
        --annotation Coverage \
        --verbosity INFO \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    """
}

/*
Below we perform QC of data.
*/

// Run Qualimap for QC metrics of recalibrated BAM.
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

