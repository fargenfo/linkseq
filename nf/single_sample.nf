#!/usr/bin/env nextflow

// TODO:
// tmp folders for various processes.

// Input parameters.
//params.fastq_paths = null
params.sample = null
params.fastq_path = null
params.reference = null
params.dbsnp = null
params.targets = null
params.threads = null
params.mem = null
params.outdir = null
params.help = false

// TODO:
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
//assert params.fastq_paths != null, 'Input parameter "fastq_paths" cannot be unasigned.'
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
//println "fastq_paths        : ${params.fastq_paths}"
println "fastq_path         : ${params.fastq_path}"
println "reference          : ${params.reference}"
println "dbsnp              : ${params.dbsnp}"
println "targets            : ${params.targets}"
println "threads            : ${params.threads}"
println "mem                : ${params.mem}"
println "outdir             : ${params.outdir}"

// Get file handlers for input files.
//fastq_paths = file(params.fastq_paths)
reference = file(params.reference)
dbsnp = file(params.dbsnp)
targets = file(params.targets)

// TODO: take in multiple paths from multiple samples.
// Turn the file with FASTQ paths into a channel with [sample, path] tuples.
//fastq_paths_ch = Channel.fromPath(fastq_paths)
//fastq_paths_ch = fastq_paths_ch.splitCsv(header: true).map { it -> [it.sample, it.fastq_path] }

// NOTE: only while developing the pipeline for running a single sample.
fastq_paths_ch = Channel.from(params.fastq_path)

process align_reads {
    //echo true

    input:
    val fastq_path from fastq_paths_ch

    output:
    file "$params.sample/outs/possorted_bam.bam" into aligned_bam_ch, aligned_bam_copy_ch

    script:
    """
    echo longranger align --id=$params.sample --sample=$params.sample \
        --reference=$reference \
        --fastqs=$fastq_path \
        --localcores=$params.threads \
        --localmem=$params.mem
    #    --uiport=3003
    mkdir -p '$params.sample/outs'
    touch '$params.sample/outs/possorted_bam.bam'
    """
}


process prepare_bqsr_table {
    //echo true

    input:
    file bam from aligned_bam_ch

    output:
    file 'bqsr.table' into bqsr_table_ch, bqsr_table_copy_ch

    script:
    """
    echo gatk BaseRecalibrator \
            -I $bam \
            -R $reference \
            --known-sites $dbsnp \
            -O 'bqsr.table' \
            --tmp-dir=tmp \
            --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    touch 'bqsr.table'
    """

}

process analyze_covariates {
    publishDir "${params.outdir}/bam/analyze_covariates", mode: 'copy', overwrite: true, saveAs: { filename -> "${params.sample}_$filename" }

    input:
    file bqsr_table from bqsr_table_ch

    output:
    file 'AnalyzeCovariates.pdf' into bqsr_analysis_ch

    script:
    """
    echo gatk AnalyzeCovariates \
        -bqsr $bqsr_table \
        -plots 'AnalyzeCovariates.pdf'
    touch 'AnalyzeCovariates.pdf'
    """
}

// TODO:
// Is it necessary to output the .bai file in order to save it with publishDir?
process apply_bqsr {
    publishDir "${params.outdir}/bam", mode: 'copy', overwrite: true, saveAs: { filename -> "${params.sample}_$filename" }

    input:
    file bqsr_table from bqsr_table_copy_ch
    file bam from aligned_bam_copy_ch

    output:
    file 'recalibrated.bam' into recalibrated_bam_ch, recalibrated_bam_copy_ch
    file 'recalibrated.bam.bai' into recalibrated_idx_ch

    script:
    """
    echo gatk ApplyBQSR \
        -R $reference \
        -I $bam \
        --bqsr-recal-file $bqsr_table \
        -L $targets \
        -O 'recalibrated.bam' \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    touch 'recalibrated.bam'
    touch 'recalibrated.bam.bai'
    """
}

process call_sample {
    input:
    file bam from recalibrated_bam_ch

    output:
    file 'haplotypecalled.gvcf' into gvcf_ch

    script:
    """
    echo gatk HaplotypeCaller  \
        -I $bam \
        -O 'haplotypecalled.gvcf' \
        -R $reference \
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
    touch 'haplotypecalled.gvcf'
    """
}

// Path to FASTQ files. The first '*' matches the Illumina flowcell ID string.
fastq_ch = Channel.fromPath("$params.fastq_path/*/$params.sample/*.fastq.gz")

// TODO:
// Because my fastqc installation is not installed at /etc/fastqc, it fails to find
// adapters, contaminants, and limits at /etc/fastqc/Configuration. The pipeline should
// assume that these files are in the correct place.
process fastqc_analysis {
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    val fastq_list from fastq_ch.toList()

    output:
    file '*.{zip,html}' into fastqc_report_ch
    file '.command.out' into fastqc_stdout_ch

    script:
    fastqs = (fastq_list as List).join(' ')
    """
    #fastqc -q --dir tmp --outdir fastqc_report
    mkdir tmp
    fastqc -q --dir tmp --outdir . -a $adapters --contaminants $contaminants -l $limits $fastqs
    """
}


process qualimap_analysis {
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    file bam from recalibrated_bam_copy_ch

    output:
    file 'qualimap_results' into qualimap_results_ch

    script:
    """
    mkdir qualimap_results
    qualimap bamqc \
        -gd HUMAN \
        -bam $bam \
        -gff $targets \
        -outdir 'qualimap_results' \
        --skip-duplicated \
        --collect-overlap-pairs \
        -nt $params.threads \
        --java-mem-size=${params.mem}G
    """

