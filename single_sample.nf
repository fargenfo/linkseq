#!/usr/bin/env nextflow
/*
Author: Ólavur Mortensen <olavur@fargen.fo>
*/

// Input parameters.
params.sample = null
params.fastq_path = null
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
assert params.sample != null, 'Input parameter "sample" cannot be unasigned.'
assert params.fastq_path != null, 'Input parameter "fastq_path" cannot be unasigned.'
assert params.reference != null, 'Input parameter "reference" cannot be unasigned.'
assert params.dbsnp != null, 'Input parameter "dbsnp" cannot be unasigned.'
assert params.targets != null, 'Input parameter "targets" cannot be unasigned.'
assert params.threads != null, 'Input parameter "threads" cannot be unasigned.'
assert params.mem != null, 'Input parameter "mem" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

println "P I P E L I N E     I P U T S    "
println "================================="
println "sample             : ${params.sample}"
println "fastq_path         : ${params.fastq_path}"
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

// Channel for the path to the FASTQ directory.
fastq_paths_ch = Channel.from(params.fastq_path)

// Align FASTQ reads to reference with LongRanger ALIGN command.
// https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/other-pipelines
process align_reads {
    memory = "${params.mem}GB"
    cpus = params.threads

    input:
    val fastq_path from fastq_paths_ch

    output:
    file "${params.sample}/outs/possorted_bam.bam" into aligned_bam_prepare_ch, aligned_bam_apply_ch

    script:
    """
    longranger align --id=${params.sample} --sample=${params.sample} \
        --reference=$reference \
        --fastqs=$fastq_path \
        --localcores=${params.threads} \
        --localmem=${params.mem}
    """
}

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
    file bam from aligned_bam_prepare_ch

    output:
    file 'bqsr.table' into bqsr_table_ch, bqsr_table_copy_ch

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

    publishDir "${params.outdir}/bam/analyze_covariates", mode: 'copy', overwrite: true, saveAs: { filename -> "${params.sample}_$filename" }

    input:
    file bqsr_table from bqsr_table_ch

    output:
    file 'AnalyzeCovariates.pdf' into bqsr_analysis_ch

    script:
    """
    gatk AnalyzeCovariates \
        -bqsr $bqsr_table \
        -plots 'AnalyzeCovariates.pdf'
    """
}

// Apply recalibration to BAM file.
process apply_bqsr {
    memory = "${params.mem}GB"
    cpus = params.threads

    publishDir "${params.outdir}/bam", mode: 'copy', overwrite: true, saveAs: { filename -> "${params.sample}_$filename" }

    input:
    file bqsr_table from bqsr_table_copy_ch
    file bam from aligned_bam_apply_ch

    output:
    file 'recalibrated.bam' into recalibrated_bam_call_ch, recalibrated_bam_qualimap_ch
    file 'recalibrated.bai' into recalibrated_idx_ch

    script:
    """
    mkdir tmp
    gatk ApplyBQSR \
        -R $reference_fa \
        -I $bam \
        --bqsr-recal-file $bqsr_table \
        -L $targets \
        -O 'recalibrated.bam' \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    """
}

// Call variants in sample with HapltypeCaller, yielding a GVCF.
process call_sample {
    memory = "${params.mem}GB"
    cpus = params.threads

    publishDir "${params.outdir}/gvcf", mode: 'copy', overwrite: true

    input:
    file bam from recalibrated_bam_call_ch

    output:
    file "${params.sample}.g.vcf" into gvcf_ch

    script:
    """
    mkdir tmp
    gatk HaplotypeCaller  \
        -I $bam \
        -O "${params.sample}.g.vcf" \
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

// Path to FASTQ files. The first '*' matches the Illumina flowcell ID string.
fastq_ch = Channel.fromPath("${params.fastq_path}/*/${params.sample}/*.fastq.gz")

// Run FastQC for QC metrics of raw data.
// Note that FastQC will allocate 250 MB of memory per thread used. Since FastQC is not a bottleneck of
// this pipeline, it will be run with a single thread.
process fastqc_analysis {
    publishDir "${params.outdir}/fastqc/${params.sample}", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    val fastq_list from fastq_ch.toList()

    output:
    file '*.{zip,html}' into fastqc_report_ch
    file '.command.out' into fastqc_stdout_ch

    script:
    fastqs = (fastq_list as List).join(' ')
    """
    mkdir tmp
    fastqc -q --dir tmp --outdir . $fastqs
    """
}

// Run Qualimap for QC metrics of aligned and recalibrated BAM.
process qualimap_analysis {
    memory = "${params.mem}GB"
    cpus = params.threads

    publishDir "${params.outdir}/bamqc", mode: 'copy',
        saveAs: {filename -> "${params.sample}"}

    input:
    file bam from recalibrated_bam_qualimap_ch

    output:
    file "qualimap_results" into qualimap_results_ch

    script:
    """
    awk 'BEGIN{OFS="\\t"}{ if(NR > 2) { print \$1,\$2,\$3,\$4,0,"." } }' $targets > 'targets_6_fields.bed'
    unset DISPLAY
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

