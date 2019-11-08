#!/usr/bin/env nextflow
/*
Author: Ólavur Mortensen <olavur@fargen.fo>
*/

// TODO:
// memory and cpu specifications for processes.

// Input parameters.
params.fastq_path = null
params.reference = null
params.targets = null
params.whitelist = null
params.bcbins = null
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
assert params.fastq_path != null, 'Input parameter "fastq_path" cannot be unasigned.'
assert params.reference != null, 'Input parameter "reference" cannot be unasigned.'
assert params.targets != null, 'Input parameter "targets" cannot be unasigned.'
assert params.whitelist != null, 'Input parameter "whitelist" cannot be unasigned.'
assert params.bcbins != null, 'Input parameter "bcbins" cannot be unasigned.'
assert params.threads != null, 'Input parameter "threads" cannot be unasigned.'
assert params.mem != null, 'Input parameter "mem" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

println "P I P E L I N E     I P U T S    "
println "================================="
println "fastq_path         : ${params.fastq_path}"
println "reference          : ${params.reference}"
println "targets            : ${params.targets}"
println "whitelist          : ${params.whitelist}"
println "bcbins             : ${params.bcbins}"
println "threads            : ${params.threads}"
println "mem                : ${params.mem}"
println "outdir             : ${params.outdir}"

// Get file handlers for input files.
reference = file(params.reference)
targets = file(params.targets)
whitelist = file(params.whitelist)

// Get lists of the read 1 and 2 FASTQ files.
// We assume two read pairs, R1 and R2, and that it is compressed. Multiple lanes work.
file_r1 = file(params.fastq_path + '/*R1*.gz')
file_r2 = file(params.fastq_path + '/*R2*.gz')

// Merge all lanes in read 1 and 2.
process merge_lanes {
    output:
    file 'R1.fastq' into fastq_r1_ch
    file 'R2.fastq' into fastq_r2_ch

    script:
    r1_list = file_r1.join(' ')
    r2_list = file_r2.join(' ')
    """
    zcat $r1_list > 'R1.fastq'
    zcat $r2_list > 'R2.fastq'
    """
}

process interleave_fastq {
    input:
    file r1 from fastq_r1_ch
    file r2 from fastq_r2_ch

    output:
    file 'interleaved.fastq' into fastq_count_ch, fastq_preproc_ch

    script:
    """
    interleave_fastq.sh $r1 $r2 > 'interleaved.fastq'
    """
}

process bc_count {
    memory = "${params.mem}GB"
    cpus = params.threads

    input:
    val fastq from fastq_count_ch

    output:
    set file('*.ema-fcnt'), file('*.ema-ncnt') into bc_count_ch

    script:
    """
    cat $fastq | ema count -w $whitelist -o bc_count
    """
}

// TODO:
// How many bins to use? There should be a params.bin
// Number of bins seems to have a large effect on how many reads end up in the "non-barcode" (nobc) bin.
// The ema GitHub recomments 500 bins.
process preproc {
    memory = "${params.mem}GB"
    cpus = params.threads

    input:
    file fastq from fastq_preproc_ch
    set file(fcnt), file(ncnt) from bc_count_ch

    output:
    file 'preproc_dir' into preproc_ema_ch, preproc_bwa_ch

    script:
    """
    cat $fastq | ema preproc -h -w $whitelist -n ${params.bcbins} -t ${params.threads} -o 'preproc_dir' $ncnt
    """
}

//rg = "@RG\tID:rg1\tSM:sample1"

// TODO:
// read group
// More threads and memory for "samtools sort"?
// This writes all the BAMs to the "preproc_dir" directory, and then merges them to a new BAM in the working directory. This is messy.
// FIXME: This creates the same number of read groups as there are bins.
process ema_align {
    memory = "${params.mem}GB"
    cpus = params.threads

    input:
    file preproc_dir from preproc_ema_ch

    output:
    file "ema_aligned.bam" into ema_bam_ch

    script:
    """
    #ema align -t ${params.threads} -d -r $reference -R "@RG\tID:ema\tSM:sample1" -s $preproc_dir/ema-bin-* | \
    #    samtools sort -@ 1 -O bam -l 0 -m 4G -o "ema_aligned.bam" -
    parallel --bar -j${params.threads} "ema align -t 1 -d -r $reference -s {} |\
        samtools sort -@ 1 -O bam -l 0 -m 4G -o {}.bam -" ::: $preproc_dir/ema-bin-???
    # The -c -p arguments make it so that the duplicate @PG and @RG tags are merged.
    samtools merge -c -p "ema_aligned.bam" $preproc_dir/ema-bin*.bam
    """
}

// TODO:
// Add read group. picard MarkDuplicates gets confused because there are \t characters in the @PG line.
process map_nobc {
    memory = "${params.mem}GB"
    cpus = params.threads

    input:
    file preproc_dir from preproc_bwa_ch

    output:
    file "bwa_nobc.bam" into nobc_bam_ch

    script:
    """
    bwa mem -p -t ${params.threads} -M -R "@RG\tID:bwa\tSM:sample1" $reference $preproc_dir/ema-nobc | \
      samtools sort -@ 1 -O bam -l 0 -m 4G -o bwa_nobc.bam
    #bwa mem -p -t ${params.threads} -M -R "@RG\tID:rg1\tSM:sample1" $reference $preproc_dir/ema-nobc | \
    #  samtools sort -@ 1 -O bam -l 0 -m 4G -o bwa_nobc.bam
    """
}

//process add_rg_nobc {
//    input:
//    file bam from nobc_bam_ch
//
//    output:
//    file "nobc_rg.bam" into nobc_rg_bam_ch
//
//    script:
//    """
//    samtools addreplacerg -r "@RG\tID:rg1\tSM:sample1" -o "nobc_rg.bam" $bam
//    """
//}


process mark_dup_ema {
    memory = "${params.mem}GB"
    cpus = params.threads

    input:
    file bam from ema_bam_ch

    output:
    file "marked_dup_ema.bam" into marked_dup_ema_ch

    script:
    """
    gatk MarkDuplicates -I $bam -O "marked_dup_ema.bam" -M "marked_dup_metrics.txt"
    """
}

process mark_dup_nobc {
    memory = "${params.mem}GB"
    cpus = params.threads

    input:
    file bam from nobc_bam_ch

    output:
    file "marked_dup_nobc.bam" into marked_dup_nobc_ch

    script:
    """
    gatk MarkDuplicates -I $bam -O "marked_dup_nobc.bam" -M "marked_dup_metrics.txt"
    """
}

process merge_bams {
    memory = "${params.mem}GB"
    cpus = params.threads

    input:
    file bam_ema from marked_dup_ema_ch
    file bam_nobc from marked_dup_nobc_ch

    output:
    file 'merged.bam' into merged_bam_ch

    script:
    """
    samtools merge "merged.bam" $bam_ema $bam_nobc
    """
}

// Run Qualimap for QC metrics of aligned and recalibrated BAM.
process qualimap_analysis {
    memory = "${params.mem}GB"
    cpus = params.threads

    input:
    file bam from merged_bam_ch

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
        -nt ${params.threads} \
        --java-mem-size=${params.mem}G
    """
}















