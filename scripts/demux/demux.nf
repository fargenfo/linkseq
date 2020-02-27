#!/usr/bin/env nextflow

params.rundir = null
params.outdir = null
params.samplesheet = null
params.help = false

helpMessage = """
    Input:
    rundir:         Path to FASTQ run directory.
    outdir:         Where to store output data.
    samplesheet:    Path to sample sheet.
    """.stripIndent()

if (params.help) {
    log.info helpMessage
    exit 0
}

assert params.rundir != null, 'Input parameter "rundir" cannot be unassigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unassigned.'
assert params.samplesheet != null, 'Input parameter "samplesheet" cannot be unassigned.'

rundir = file(params.rundir)
outdir = file(params.outdir)
samplesheet = file(params.samplesheet)
interop_dir = file(outdir + "InterOp")

println "D E M U X    L I N K    "
println "================================="
println "rundir              : ${rundir}"
println "outdir              : ${outdir}"
println "samplesheet         : ${samplesheet}"
println "interop_dir         : ${interop_dir}"

// Call bcl2fastq, performing simultaneous basecalling and demultiplexing.
// --use-bases-mask will use RunInfo.xml (in the run directory) to determine the length of read 1 and 2
// and of the index.
// Adapter sequences (read 1 and read2) should be contained in the sample sheet.
process bcl2fastq {
    publishDir "$outdir", mode: 'copy', pattern: '.command.log', saveAs: {filename -> 'bcl2fastq.log'}

    output:
    file "outs/*fastq.gz" into fastq_samplenames_ch
    file '.command.log'

    script:
    if(task.cpus > 20) {
        p_threads = task.cpus - 8
        w_threads = 4
        r_threads = 4
    } else if(task.cpus > 10) {
        p_threads = task.cpus - 2
        w_threads = 1
        r_threads = 1
    } else if(task.cpus > 3) {
        p_threads = task.cpus - 2
        w_threads = 1
        r_threads = 1
    } else {
        p_threads = 1
        w_threads = 1
        r_threads = 1
    }
    """
    bcl2fastq \
        -R $rundir \
        -o outs \
        --interop-dir interop \
        --sample-sheet $samplesheet \
        --use-bases-mask Y*,I*,Y* \
        --minimum-trimmed-read-length 8 \
        --mask-short-adapter-reads 8 \
        --ignore-missing-positions \
        --ignore-missing-filter \
        --ignore-missing-bcls \
        -p $p_threads -r $r_threads -w $w_threads
    """
}

// Prepare the channel for the merge process.

// Get (key, FASTQ files) tuples, where key is a (sample names, lane, read) tuple.
// First flatten the channel because each instance of process "bcl2fastq" outputs a tuple.
// Then map the channel to (key, FASTQ path) tuples. This channel has one record per file.
// Then group by key to get (key, FASTQ list) tuples.
fastq_ch = fastq_samplenames_ch.flatten()
    .map { file ->
        def sample = file.name.toString().split('_')[0]
        def lane = file.name.toString().split('_')[2]
        def read = file.name.toString().split('_')[3]
        def key = tuple(sample, lane, read)
        return tuple(key, file)}
    .groupTuple()

// Since 10x samples have multiple indexes per sample, we merge these.
// Since this process is only concatenating (cat) and zipping files, it doesn't need much memory or many cores.
// We use the "when" directive to avoid processing the "Undetermined" sample.
process merge {
    input:
    set key, file(fastqs) from fastq_ch

    output:
    set sample, file("*merged.fastq.gz") into fastq_sync_ch

    when:
    key[0] != "Undetermined"

    script:
    sample = key[0]
    lane = key[1]
    read = key[2]
    // Sort the FASTQ names so that they are merged in the proper order.
    fastqs = (fastqs as List)
    fastqs.sort()
    fastqs = fastqs.join(' ')
    """
    # Note: Piping the zcat output to gzip causes "unexpected end of file" errors sporadically.
    # Therefore, the zcat and gzip are done in separate steps.
    zcat $fastqs > $sample\\_$lane\\_$read\\_merged.fastq
    gzip -k $sample\\_$lane\\_$read\\_merged.fastq
    """
}

// Prepare channel for sync_reads.

// Get (key, FASTQ files) tuples, where key is a (samplename, lane) tuple.
// Then map the channel to (key, FASTQ path) tuples. This channel has one record per file.
// Then group by key to get (key, FASTQ list) tuples.
fastq_sync_ch = fastq_sync_ch.map { it ->
        // The record is a (sample, FASTQ file) tuple.
        def sample = it[0]
        def file = it[1]
        // Process the file string to get the lane number.
        def lane = file.name.toString().split('_')[1]
        def key = tuple(sample, lane)
        return tuple(key, file)}
    .groupTuple()

// Sort the file names, so that the list is always in order (read1, read2).
fastq_sync_ch = fastq_sync_ch.map { it ->
    def key = it[0]
    def fastqs = it[1]
    return tuple(key, fastqs.sort())}

// Synchronize reads, in case the reads got out of order in the merge process.
process sync_reads {
    publishDir "$outdir/$sample/fastqs", mode: 'copy'

    input:
    set key, file(fastqs) from fastq_sync_ch

    output:
    set key, file("*synced.fastq.gz") into fastq_qc_ch

    script:
    sample = key[0]
    lane = key[1]
    read1 = fastqs[0]
    read2 = fastqs[1]
    """
    # We use ziplevel=1 to get fast but low-level compression.
    # NOTE: singletons.fastq.gz should be empty.
    repair.sh -Xmx${task.memory.toGiga()}g ziplevel=1 in1=$read1 in2=$read2 out1=$sample\\_$lane\\_R1\\_synced.fastq.gz out2=$sample\\_$lane\\_R2\\_synced.fastq.gz outs=singletons.fastq.gz repair
    """
}

// Run FastQC for QC metrics of raw data.
process fastqc_analysis {
    publishDir "$outdir/$sample/fastqc", mode: 'copy', pattern: '{*.zip,*.html}',
        saveAs: {filename -> filename.indexOf('.zip') > 0 ? "zips/$filename" : "$filename"}
    publishDir "$outdir/$sample/fastqc", mode: 'copy', pattern: '.command.log',
        saveAs: {filename -> 'fastqc.log'}

    input:
    set key, file(fastqs) from fastq_qc_ch

    output:
    set sample, file('*.{zip,html}') into fastqc_report_ch
    set sample, file('.command.log') into fastqc_stdout_ch

    script:
    fastq_list = (fastqs as List).join(' ')
    sample = key[0]
    lane = key[1] // Unused.
    """
    # We unset the DISPLAY variable to avoid having FastQC try to open the GUI.
    unset DISPLAY
    mkdir tmp
    fastqc -q --dir tmp --threads ${task.cpus} --outdir . $fastq_list
    """
}


