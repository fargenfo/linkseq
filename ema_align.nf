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
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

println "P I P E L I N E     I P U T S    "
println "================================="
println "fastq_path         : ${params.fastq_path}"
println "reference          : ${params.reference}"
println "targets            : ${params.targets}"
println "whitelist          : ${params.whitelist}"
println "bcbins             : ${params.bcbins}"
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
    input:
    file fastq from fastq_count_ch

    output:
    set file('*.ema-fcnt'), file('*.ema-ncnt') into bc_count_ch

    script:
    """
    cat $fastq | ema count -w $whitelist -o bc_count
    """
}

// TODO:
// How many bins to use?
// Number of bins seems to have a large effect on how many reads end up in the "non-barcode" (nobc) bin.
// The ema GitHub recomments 500 bins.
process preproc {
    input:
    file fastq from fastq_preproc_ch
    set file(fcnt), file(ncnt) from bc_count_ch

    output:
    file "preproc_dir/ema-bin-*" into bins_ema_ch mode flatten
    file "preproc_dir/ema-nobc" into nobc_bin_bwa_ch

    script:
    """
    cat $fastq | ema preproc -h -w $whitelist -n ${params.bcbins} -t ${task.cpus} -o 'preproc_dir' $ncnt
    """
}

// It is important that the input to EMA (and BWA further down) contains literal "\t" strings, not actual tabs.
// Using the "/" character here and enclosing it in quotes, "$READGROUP_EMA" in the processes achieves this.
READGROUP_EMA = /@RG\tID:ema\tSM:sample1/


// TODO:
// read group
// More threads and memory for "samtools sort"?
// Maybe better to sort the BAM in a separate process.
// Why do we sort the BAM at this stage? Maybe it would be better to sort once everything is merged.
process ema_align {
    input:
    file bin from bins_ema_ch

    output:
    file "${bin}.bam" into ema_bam_ch

    script:
    """
    ema align -t ${task.cpus} -d -r $reference -R '$READGROUP_EMA' -s $bin | \
        samtools view -b -o ${bin}.bam
    """
}

READGROUP_BWA = /@RG\tID:bwa\tSM:sample1/

process map_nobc {
    input:
    file nobc_bin from nobc_bin_bwa_ch

    output:
    file "nobc.bam" into nobc_bam_ch

    script:
    """
    bwa mem -p -t ${task.cpus} -M -R "$READGROUP_BWA" $reference $nobc_bin | \
        samtools view -b -o nobc.bam
    """
}

// Merge BAMs from both EMA and BWA.
// All BAMs from EMA bins have the same readgroup, so the RG and PG headers wil be combined.
process merge_bams {
    input:
    file ema_bams from ema_bam_ch.collect()
    file bwa_bam from nobc_bam_ch

    output:
    file "merged.bam" into merged_bam_sort_ch

    script:
    ema_bam_list = (ema_bams as List).join(' ')
    """
    samtools merge -@ ${task.cpus} -O bam -l 0 -c -p "merged.bam" $ema_bam_list $bwa_bam
    """
}

process sort_bams {
    input:
    file bam from merged_bam_sort_ch

    output:
    file "sorted.bam" into sorted_bam_markdup_ch

    script:
    """
    samtools sort -@ ${task.cpus} -O bam -l 0 -m 4G -o sorted.bam $bam
    """
}


// NOTE:
// MarkDuplicates has the following option, I wonder why:
// --BARCODE_TAG:String          Barcode SAM tag (ex. BC for 10X Genomics)  Default value: null.                          
process mark_dup {
    input:
    file bam from sorted_bam_markdup_ch

    output:
    file "marked_dup.bam" into marked_bam_qc_ch

    script:
    """
    gatk MarkDuplicates -I $bam -O "marked_dup.bam" -M "marked_dup_metrics.txt"
    """
}

/*


// FIXME:
// Qualimap gives some Java error related to fonts.
// Run Qualimap for QC metrics of aligned and recalibrated BAM.
process qualimap_analysis {
    input:
    file bam from marked_bam_qc_ch

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

*/








