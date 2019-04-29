#!/usr/bin/env nextflow

// Input parameters.
params.reference = null
params.dbsnp = null
params.targets = null
params.threads = null
params.mem = null
params.outdir = null
params.help = false

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
assert params.reference != null, 'Input parameter "reference" cannot be unasigned.'
assert params.dbsnp != null, 'Input parameter "dbsnp" cannot be unasigned.'
assert params.targets != null, 'Input parameter "targets" cannot be unasigned.'
assert params.threads != null, 'Input parameter "threads" cannot be unasigned.'
assert params.mem != null, 'Input parameter "mem" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

println "P I P E L I N E     I P U T S    "
println "================================="
println "reference          : ${params.reference}"
println "dbsnp              : ${params.dbsnp}"
println "targets            : ${params.targets}"
println "threads            : ${params.threads}"
println "mem                : ${params.mem}"
println "outdir             : ${params.outdir}"

// Get file handlers for input files.
//fastq_paths = file(params.fastq_paths)
reference = file(params.reference)  // Directory of 10x reference.
reference_fa = file(params.reference + '/fasta/genome.fa')  // Reference fasta file.
dbsnp = file(params.dbsnp)
targets = file(params.targets)

process consolidate_gvcf {

    script:
    """
    export TILEDB_DISABLE_FILE_LOCKING=1
    gatk GenomicsDBImport \
        -V $vcf
        -L $targets \
        --genomicsdb-workspace-path $genomicsdb \
        --merge-input-intervals \
        --tmp-dir=tmp \
        --java-options "-Xmx50g -Xms50g"
    """
}

process joint_genotyping {

    script:
    """
    gatk GenotypeGVCFs \
        -V gendb://$genomicsdb \
        -R $reference \
        -O $vcf \
        --tmp-dir=tmp \
        --java-options "-Xmx50g -Xms50g"
    """
}

