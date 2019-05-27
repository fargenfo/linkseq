#!/usr/bin/env nextflow
/*
Author: Ólavur Mortensen <olavur@fargen.fo>
*/

// Input parameters.
params.fastq_paths = null
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
assert params.fastq_paths != null, 'Input parameter "fastq_paths" cannot be unasigned.'
assert params.reference != null, 'Input parameter "reference" cannot be unasigned.'
assert params.dbsnp != null, 'Input parameter "dbsnp" cannot be unasigned.'
assert params.targets != null, 'Input parameter "targets" cannot be unasigned.'
assert params.threads != null, 'Input parameter "threads" cannot be unasigned.'
assert params.mem != null, 'Input parameter "mem" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

println "P I P E L I N E     I P U T S    "
println "================================="
println "fastq_paths         : ${params.fastq_paths}"
println "reference          : ${params.reference}"
println "dbsnp              : ${params.dbsnp}"
println "targets            : ${params.targets}"
println "threads            : ${params.threads}"
println "mem                : ${params.mem}"
println "outdir             : ${params.outdir}"

// Get file handlers for input files.
reference = file(params.reference)  // Directory of 10x reference.
reference_fa = file(params.reference + '/fasta/genome.fa')  // Reference fasta file.
dbsnp = file(params.dbsnp)
targets = file(params.targets)



// FIXME
// remove when done testing
targets = "chr17"
// FIXME




// Turn the file with FASTQ paths into a channel with [sample, path] tuples.
fastq_paths_ch = Channel.fromPath(params.fastq_paths)
fastq_paths_ch
    .splitCsv(header: true)
    .map { it -> tuple(it.sample, it.fastq_path) }
    .into { fastq_align_ch; fastq_print_ch; fastq_qc_ch }

println("Processing data:\nSample\tFASTQ path")
fastq_print_ch.subscribe { println(it[0] + "\t" + it[1]) }

// Channel for the path to the FASTQ directory.
//fastq_paths_ch = Channel.from(params.fastq_path)

// TODO: point directly to fastq folder
// Align FASTQ reads to reference with LongRanger ALIGN command.
// https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/other-pipelines
process align_reads {
    memory = "${params.mem}GB"
    cpus = params.threads

    input:
    //val fastq_path from fastq_align_ch
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

process make_small_bam {
    input:
    set sample, file(bam), file(bai) from aligned_bam_prepare_ch

    output:
    //set file("small.bam"), file("small.bam.bai") into small_bam_genotyping_ch, small_bam_extract_ch, small_bam_link_ch, small_bam_phase_bam_ch
    set sample, file("small.bam"), file("small.bam.bai") into small_bam_ch

    script:
    """
    samtools view -b -o "small.bam" -T $reference_fa $bam "chr17"
    samtools index -b "small.bam"
    """
}

aligned_bam_prepare_ch = null
aligned_bam_apply_ch = null
small_bam_ch.into { aligned_bam_prepare_ch; aligned_bam_apply_ch }

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
    //set sample, file(bqsr_table) from bqsr_table_copy_ch
    //set sample, file(bam), file(bai) from aligned_bam_apply_ch
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

// Prepare input for FastQC. FastQC needs the paths to FASTQ files. Below, we get each path in a list,
// and then we map the list to a string.
fastq_concat_ch = fastq_qc_ch
    .map { tuple(it[0], file(it[1] + "/*.fastq{.gz,}")) }
    .map { tuple(it[0], it[1].join(' ')) }

// Concatenate all the FASTQ files of one sample into one FASTQ file. This way, we get only one FastQC
// report. Also, conveniently, we avoid using the original file names, possibly containing confidiential
// sample names, in the FastQC report.
process concat_fastq {
    input:
    set sample, fastq_list from fastq_concat_ch

    output:
    set sample, file("${sample}.fastq.gz") into fastq_fastqc_ch

    script:
    """
    zcat $fastq_list | bgzip -c > "${sample}.fastq.gz"
    """
}

// Run FastQC for QC metrics of raw data.
// Note that FastQC will allocate 250 MB of memory per thread used. Since FastQC is not a bottleneck of
// this pipeline, it will be run with a single thread.
process fastqc_analysis {
    publishDir "${params.outdir}/fastqc/${sample}", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set sample, fastq from fastq_fastqc_ch

    output:
    set sample, file('*.{zip,html}') into fastqc_report_ch
    set sample, file('.command.out') into fastqc_stdout_ch

    script:
    """
    mkdir tmp
    fastqc -q --dir tmp --outdir . $fastq
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
    #awk 'BEGIN{OFS="\\t"}{ if(NR > 2) { print \$1,\$2,\$3,\$4,0,"." } }' $targets > 'targets_6_fields.bed'

    # FIXME:
    # remove this when done testing.
    # and uncomment awk command above
    echo 'track name="dummy" description="dummy BED" color=0,0,128 db=hg38' > 'targets_6_fields.bed'
    echo 'chr17    1  83257441 allchr17   0   .' >> 'targets_6_fields.bed'
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

