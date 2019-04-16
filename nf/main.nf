#!/usr/bin/env nextflow



// Input parameters.
params.fastq_paths = null
params.reference = null
params.threads = null
params.mem = null
params.outdir = null
params.help = false

// TODO:
helpMessage = """ 
Parameters:
--fastq_paths       A file with FASTQ paths.
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
assert params.threads != null, 'Input parameter "threads" cannot be unasigned.'
assert params.mem != null, 'Input parameter "mem" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

println "P I P E L I N E     I P U T S    "
println "================================="
println "fastq_paths        : ${params.fastq_paths}"
println "reference          : ${params.reference}"
println "threads            : ${params.threads}"
println "mem                : ${params.mem}"
println "outdir             : ${params.outdir}"

// Get file handlers for input files.
fastq_paths_file = file(params.fastq_paths)

// Turn the file with FASTQ paths into a channel with [sample, path] tuples.
fastq_paths_ch = Channel.fromPath(fastq_paths_file)
fastq_paths_ch = fastq_paths_ch.splitCsv(header: true).map { it -> [it.sample, it.fastq_path] }

process align_reads {
    tag "$sample"

    input:
    set sample, fastq_path from fastq_paths_ch

    output:
    file 'fastq_path.txt' into temp_ch

    script:
    """
    longranger align --id=$sample \
        --reference=$reference \
        --fastqs=$fastq_path \
        --localcores=$threads \
        --localmem=$mem \
        --uiport=3003
    """
}
