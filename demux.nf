#!/usr/bin/env nextflow

params.rundir = null
params.outdir = null
params.samplesheet = null
params.whitelist = null
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
assert params.whitelist != null, 'Input parameter "whitelist" cannot be unassigned.'

rundir = file(params.rundir)
outdir = file(params.outdir)
samplesheet = file(params.samplesheet)
whitelist = file(params.whitelist)
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
    set sample, file("*merged.fastq.gz") into fastq_check_sync_ch

    when:
    key[0] != "Undetermined"

    script:
    sample = key[0]
    lane = key[1]
    read = key[2]
    // If there are multiple FASTQs in the input, sort the names so that they are merged in the proper order.
    if(fastqs instanceof List) {
        fastqs = (fastqs as List)
        fastqs.sort()
        fastqs = fastqs.join(' ')
    }
    """
    # Note: Piping the zcat output to gzip causes "unexpected end of file" errors sporadically.
    # Therefore, the zcat and gzip are done in separate steps.
    zcat $fastqs > $sample\\_$lane\\_$read\\_merged.fastq
    gzip -k $sample\\_$lane\\_$read\\_merged.fastq
    """
}

process mess_up_sync_test {
   input:
   

// Prepare input channel for check_sync.
fastq_ch = fastq_check_sync_ch.map { it ->
        def sample = it[0]
        def file = it[1]
        def lane = file.name.toString().split('_')[2]
        def key = tuple(sample, lane)
        return tuple(key, file)}
    .groupTuple()

process check_sync {
    input:
    set key, file(fastqs) from fastq_sync_ch

    output:
    file check_sync.log into check_sync_log_ch

    script:
    sample = key[0]
    lane = key[1]
    read1 = fastqs[0]
    read2 = fastqs[1]
    """
    # Check if reads are synchronized.
    reformat.sh in="*R1*.fastq.gz" in2="*R2*.fastq.gz" out="check_sync.log" vpair
    """
}

// NOTE: the code below synchronizes read 1 and 2 in FASTQ.
//
//// Prepare channel for sync_reads.
//
//// Get (key, FASTQ files) tuples, where key is a (samplename, lane) tuple.
//// Then map the channel to (key, FASTQ path) tuples. This channel has one record per file.
//// Then group by key to get (key, FASTQ list) tuples.
//fastq_sync_ch = fastq_sync_ch.map { it ->
//        // The record is a (sample, FASTQ file) tuple.
//        def sample = it[0]
//        def file = it[1]
//        // Process the file string to get the lane number.
//        def lane = file.name.toString().split('_')[1]
//        def key = tuple(sample, lane)
//        return tuple(key, file)}
//    .groupTuple()
//
//// Sort the file names, so that the list is always in order (read1, read2).
//fastq_sync_ch = fastq_sync_ch.map { it ->
//    def key = it[0]
//    def fastqs = it[1]
//    return tuple(key, fastqs.sort())}
//// Synchronize reads, in case the reads got out of order in the merge process.
//process sync_reads {
//    publishDir "$outdir/$sample/fastqs", mode: 'copy'
//
//    input:
//    set key, file(fastqs) from fastq_sync_ch
//
//    output:
//    set key, file("*R1*synced.fastq.gz"), file("*R2*synced.fastq.gz") into fastq_trim_adapters_ch
//
//    script:
//    sample = key[0]
//    lane = key[1]
//    read1 = fastqs[0]
//    read2 = fastqs[1]
//    """
//    # We use ziplevel=1 to get fast but low-level compression.
//    # NOTE: singletons.fastq.gz should be empty.
//    repair.sh -Xmx${task.memory.toGiga()}g ziplevel=1 in1=$read1 in2=$read2 out1=$sample\\_$lane\\_R1\\_synced.fastq.gz out2=$sample\\_$lane\\_R2\\_synced.fastq.gz outs=singletons.fastq.gz repair
//    """
//}

// Parses samplesheet and saves adapter in a FASTA file.
process extract_adapter {
    output:
    file "adapter.fasta" into adapter_fasta_ch

    script:
    """
    samplesheet_extract_adapter.py $samplesheet > adapter.fasta
    """
}

// Combine the FASTQs with the adapter FASTA file to get (key, read1 FASTQ, read2 FASTQ, adapter FASTA) tuples.
trim_adapters_data_ch = fastq_trim_adapters_ch.combine(adapter_fasta_ch)

// FIXME: output log
// Trim adapters.
process trim_adapters {
    input:
    set key, file(read1), file(read2), file(adapter_fasta) from trim_adapters_data_ch

    output:
    set key, file("*R1*adapter_trimmed.fastq.gz"), file("*R2*adapter_trimmed.fastq.gz") into fastq_bctrim_ch

    script:
    sample = key[0]
    lane = key[1]
    """
    # Trim adapters from 3' end (ktrim=r) with up to 2 mismatches (hdist=2).
    # k-mer size 21, and 11 at the end of the read (mink=11).
    # Use pair overlap detection (tbo), and trim both reads to the same length (tpe).
    bbduk.sh in1=$read1 in2=$read2 out1=$sample\\_$lane\\_R1\\_adapter_trimmed.fastq.gz out2=$sample\\_$lane\\_R2\\_adapter_trimmed.fastq.gz ref=$adapter_fasta ktrim=r k=21 mink=11 hdist=2 tbo tpe 2> bbduk.log
    """
}

//// FIXME:
//// Check if reads are synchronized. If they are not, exit with an error.
//
//// FIXME: 
//// Trim 10x barcode from read 2.
//// The barcode is taken from the first 16 bases of read 1.
//// If the barcode does not match any in the list of known barcodes (whitelist), we do not trim.
//process bctrim {
//    input:
//    output:
//    script:
//    """
//    trimR2bc.py $read1 $read2 $whitelist [OUTPUT R2] 1> bctrim_stats.log
//    """
//}
//
//// FIXME: remember to take read2 from bctrim, and read1 from... the previous process...
//
//// Trim poly-G tail.
//process polyG_trim {
//    input:
//    output:
//    script:
//    """
//    # Trim poly G of reads
//    # Q: Disable quality filter, -L: Disable length filter, -A: Disable adapter filtering, -g: Enable polyG trim
//    fastp -i $read1 -I $read2 -o $read1_out -O $read2_out -Q -L -A -g -h "polyG_trim_log.html" -j "polyG_trim_log.json" 2> polyG_trim.log
//    # NOTE: is does not seem like HTML and JSON reports (-h and -j) can be disabled.
//    """
//}
//
//// FIXME: in polyG_trim, output to two channels (read 1 and read 2).
//
//// Do quality trimming and minimum length filtering (read 1).
//process quality_trim_read1 {
//    input:
//    output:
//    script:
//    """
//    # We don't trim from 5' end, because this would trim the barcode (-x).
//    sickle se -f $file -t sanger -g -o $out -x -q 20 -l 58 1> sickle.log
//    """
//}
//
//// Do quality trimming and minimum length filtering (read 2).
//process quality_trim_read2 {
//    input:
//    output:
//    script:
//    """
//    sickle se -f $file -t sanger -g -o $out -q 20 -l 35 1> sickle.log
//    """
//}
//
//// FIXME: combing read 1 and 2 channels.
//
//process sync_reads_qtrim {
//    input:
//    output:
//    script:
//    """
//    # We use ziplevel=1 to get fast but low-level compression.
//    # NOTE: singletons.fastq.gz should be empty.
//    # FIXME: output names????
//    repair.sh -Xmx${task.memory.toGiga()}g ziplevel=1 in1=$read1 in2=$read2 out1=$sample\\_$lane\\_R1\\_synced.fastq.gz out2=$sample\\_$lane\\_R2\\_synced.fastq.gz outs=singletons.fastq.gz repair
//    """
//}
//
//// Run FastQC for QC metrics of raw data.
//process fastqc_analysis {
//    publishDir "$outdir/$sample/fastqc", mode: 'copy', pattern: '{*.zip,*.html}',
//        saveAs: {filename -> filename.indexOf('.zip') > 0 ? "zips/$filename" : "$filename"}
//    publishDir "$outdir/$sample/fastqc", mode: 'copy', pattern: '.command.log',
//        saveAs: {filename -> 'fastqc.log'}
//
//    input:
//    set key, file(fastqs) from fastq_qc_ch
//
//    output:
//    set sample, file('*.{zip,html}') into fastqc_report_ch
//    set sample, file('.command.log') into fastqc_stdout_ch
//
//    script:
//    fastq_list = (fastqs as List).join(' ')
//    sample = key[0]
//    lane = key[1] // Unused.
//    """
//    # We unset the DISPLAY variable to avoid having FastQC try to open the GUI.
//    unset DISPLAY
//    mkdir tmp
//    fastqc -q --dir tmp --threads ${task.cpus} --outdir . $fastq_list
//    """
//}


