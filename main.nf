#!/usr/bin/env nextflow
/*
Author: Ólavur Mortensen <olavur@fargen.fo>
*/


/*
TODO:

If one sample fails, and I remove it from the CSV, is Nextflow able to get the other runned samples from
the cache? Or will all samples run from the start?

*/

// Input parameters.
params.fastq_csv = null
params.reference = null
params.targets = null
params.whitelist = null
params.bcbins = null
params.dbsnp = null
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
assert params.fastq_csv != null, 'Input parameter "fastq_csv" cannot be unasigned.'
assert params.reference != null, 'Input parameter "reference" cannot be unasigned.'
assert params.targets != null, 'Input parameter "targets" cannot be unasigned.'
assert params.whitelist != null, 'Input parameter "whitelist" cannot be unasigned.'
assert params.bcbins != null, 'Input parameter "bcbins" cannot be unasigned.'
assert params.dbsnp != null, 'Input parameter "dbsnp" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

println "P I P E L I N E     I P U T S    "
println "================================="
println "fastq_csv          : ${params.fastq_csv}"
println "reference          : ${params.reference}"
println "targets            : ${params.targets}"
println "whitelist          : ${params.whitelist}"
println "bcbins             : ${params.bcbins}"
println "dbsnp              : ${params.dbsnp}"
println "outdir             : ${params.outdir}"
println '=================================='

// Get file handlers for input files.
reference = file(params.reference, checkIfExists: true)
targets = file(params.targets, checkIfExists: true)
whitelist = file(params.whitelist, checkIfExists: true)
dbsnp = file(params.dbsnp, checkIfExists: true)
outdir = file(params.outdir)

// Read FASTQ read 1 and 2, as well as sample IDs, from input CSV file.
Channel.fromPath(params.fastq_csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sample, file(row.read1), file(row.read2)) }
    .into { fastq_ch; fastq_print_ch }

println 'Sample ID\tFASTQ read 1 files\tFASTQ read 2 files'
fastq_print_ch.subscribe onNext: { row ->
    sample = row[0]
    read1 = row[1].name.join(',')
    read2 = row[2].name.join(',')
    println sample + "\t"  + read1 + "\t" + read2
    }, onComplete: {
    println '==================================' }

/*
First, we align the data to reference with EMA. In order to do so, we need to do some pre-processing, including,
but not limited to, merging lanes, counting barcodes, and binning reads.
*/

// Merge all lanes in read 1 and 2.
// If there is only one lane, all this process does is decompress the files.
process merge_lanes {
    input:
    set sample, file(read1), file(read2) from fastq_ch

    output:
    set sample, file('R1.fastq.gz'), file('R2.fastq.gz') into merged_fastq_ch

    script:
    // If there are multiple input FASTQs, sort the list of FASTQ by file name, so that they are
    // concatenated in the same order. Join the names in a single string also.
    if(read1 instanceof List) {
        read1 = read1.sort{ it.name }
        read2 = read2.sort{ it.name }
        read1 = read1.join(' ')
        read2 = read2.join(' ')
    }
    script:
    """
    zcat $read1 | gzip -c > 'R1.fastq.gz'
    zcat $read2 | gzip -c > 'R2.fastq.gz'
    """
}

// Interleave reads 1 and 2.
process interleave_fastq {
    input:
    set sample, file(read1), file(read2) from merged_fastq_ch

    output:
    set sample, file('interleaved.fastq.gz') into fastq_count_ch, fastq_preproc_ch, fastq_readgroup_ch, fastq_check_sync_ch

    script:
    """
    reformat.sh in=$read1 in2=$read2 out=interleaved.fastq.gz
    """
}

// In the unlikely event that either the merging or interleaving procedures went wrong, this process
// will see that the reads are out of sync and throw an error.
process check_sync {
    input:
    set sample, file(fastq) from fastq_check_sync_ch

    output:
    set sample, val('done') into check_sync_status_ch

    script:
    """
    reformat.sh in=$fastq vint
    """
}

fastq_count_ch.join(check_sync_status_ch).set{data_bc_count_ch}

// Count barcodes in FASTQ.
process bc_count {
    input:
    set sample, file(fastq), status from data_bc_count_ch

    output:
    set sample, file('*.ema-fcnt'), file('*.ema-ncnt') into bc_count_ch

    script:
    """
    zcat $fastq | ema count -w $whitelist -o bc_count
    """
}

fastq_preproc_ch.join(bc_count_ch).set{data_preproc_ch}

// Statistical binning of reads, splitting the reads into bins.
// TODO:
// How many bins to use?
// Number of bins seems to have a large effect on how many reads end up in the "non-barcode" (nobc) bin.
// The ema GitHub recomments 500 bins.
// Barcode correction report?
process preproc {
    input:
    set sample, file(fastq), file(fcnt), file(ncnt) from data_preproc_ch

    output:
    set sample, file("preproc_dir/ema-bin-*") into bins_ema_ch mode flatten
    set sample, file("preproc_dir/ema-nobc") into nobc_bin_bwa_ch

    script:
    """
    zcat $fastq | ema preproc -h -w $whitelist -n ${params.bcbins} -t ${task.cpus} -o 'preproc_dir' $ncnt
    """
}

// Construct a readgroup from the sequence identifier in one of the input FASTQ files.
process get_readgroup {
    input:
    set sample, file(fastq) from fastq_readgroup_ch

    output:
    set val(sample), stdout into readgroup_ema_ch, readgroup_bwa_ch

    script:
    """
    get_readgroups.py $fastq $sample
    """
}

// Combine the readgroup channel with the EMA bins channel so that each instance of the ema_align process gets
// a readgroup object.
// First combine all (sample, rg) pairs with all (sample, bin) pairs.
readgroup_ema_ch.combine(bins_ema_ch).set{data_ema_ch}
// Filter out all pairs where the sample ID doesn't match.
data_ema_ch.filter { it[0] == it[2] }.set{data_ema_ch}
// Map from (sample, rg, sample, bin) to (sample, rg, bin).
data_ema_ch.map { tuple(it[0], it[1], it[3]) }.set{data_ema_ch}

// Align reads from each bin with EMA.
process ema_align {
    input:
    set sample, rg, file(bin) from data_ema_ch

    output:
    set sample, file("${bin}.bam") into ema_bam_ch

    script:
    """
    ema align -t ${task.cpus} -d -r $reference -R '$rg' -s $bin | \
        samtools view -b -o ${bin}.bam
    """
}

readgroup_bwa_ch.join(nobc_bin_bwa_ch).set{data_bwa_ch}

// Align the no-barcode bin. These reads had barcodes that didn't match the whitelist.
process map_nobc {
    input:
    set sample, rg, file(nobc_bin) from data_bwa_ch

    output:
    set sample, file("nobc.bam") into nobc_bam_ch

    script:
    """
    bwa mem -p -t ${task.cpus} -M -R '$rg' $reference $nobc_bin | \
        samtools view -b -o nobc.bam
    """
}

// Group EMA bin BAMs by key, such that we have (sample, BAM list) tuples.
ema_bam_ch.groupTuple().set{ema_bam_grouped_ch}
// Combine this with the no-barcode bin, yielding a (sample, no-bc BAM, EMA bin BAM list) tuples.
nobc_bam_ch.join(ema_bam_grouped_ch).set{data_merge_bams_ch}

// Merge BAMs from both EMA and BWA.
// All BAMs have the same readgroup, so the RG and PG headers wil be combined.
process merge_bams {
    input:
    set sample, nobc_bam, ema_bams from data_merge_bams_ch

    output:
    set sample, file("merged.bam") into merged_bam_sort_ch

    script:
    ema_bam_list = (ema_bams as List).join(' ')
    """
    samtools merge -@ ${task.cpus} -O bam -l 0 -c -p "merged.bam" $nobc_bam $ema_bam_list
    """
}

// Coordinate sort BAM.
process sort_bam {
    input:
    set sample, file(bam) from merged_bam_sort_ch

    output:
    set sample, file("sorted.bam") into sorted_bam_markdup_ch

    script:
    """
    samtools sort -@ ${task.cpus} -O bam -l 0 -m 4G -o "sorted.bam" $bam
    """
}

// Mark duplicates in BAM.
// NOTE:
// MarkDuplicates has the following option, I wonder why:
// --BARCODE_TAG:String          Barcode SAM tag (ex. BC for 10X Genomics)  Default value: null.                          
process mark_dup {
    input:
    set sample, file(bam) from sorted_bam_markdup_ch

    output:
    set sample, file("marked_dup.bam") into marked_bam_index_ch

    script:
    """
    gatk MarkDuplicates -I $bam -O "marked_dup.bam" -M "marked_dup_metrics.txt"
    """
}

// Index the BAM.
process index_bam {
    input:
    set sample, file(bam) from marked_bam_index_ch

    output:
    set sample, file("$bam"), file("${bam}.bai") into indexed_bam_prepare_ch, indexed_bam_apply_ch

    script:
    """
    gatk BuildBamIndex -I $bam -O "${bam}.bai"
    """
}

/*
The next three processes, prepare_bqsr_table, analyze_covariates, and apply_bqsr, deal with base quality score
recalibration, in preparation for GATK best practices.
BQSR: https://software.broadinstitute.org/gatk/documentation/article?id=44
*/

// Generate recalibration table for BQSR.
process prepare_bqsr_table {
    publishDir "$outdir/$sample/bam/bqsr", mode: 'copy', overwrite: true

    input:
    set sample, file(bam), file(bai) from indexed_bam_prepare_ch

    output:
    set sample, file('bqsr.table') into bqsr_table_analyze_ch, bqsr_table_apply_ch

    script:
    """
    mkdir tmp
    gatk BaseRecalibrator \
            -I $bam \
            -R $reference \
            -L $targets \
            --known-sites $dbsnp \
            -O 'bqsr.table' \
            --tmp-dir=tmp \
            --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g"
    """
}

indexed_bam_apply_ch.join(bqsr_table_apply_ch).set{data_apply_ch}

// Apply recalibration to BAM file.
// NOTE: this BAM will be phased at a later stage.
process apply_bqsr {
    input:
    set sample, file(bam), file(bai), bqsr_table from data_apply_ch

    output:
    set sample, file("recalibrated.bam"), file("recalibrated.bam.bai") into recalibrated_bam_call_ch, recalibrated_bam_second_pass_ch, bam_phase_vcf_ch, bam_phase_bam_ch

    script:
    """
    mkdir tmp
    gatk ApplyBQSR \
        -R $reference \
        -I $bam \
        --bqsr-recal-file $bqsr_table \
        -L $targets \
        -O "recalibrated.bam" \
        --tmp-dir=tmp \
        --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g"
    mv "recalibrated.bai" "recalibrated.bam.bai"
    """
}

// Second pass of BQSR, giving a "before and after" picture of BQSR.
process bqsr_second_pass {
    publishDir "$outdir/$sample/bam/bqsr", mode: 'copy', overwrite: true

    input:
    set sample, file(bam), file(bai) from recalibrated_bam_second_pass_ch

    output:
    set sample, file('bqsr_second_pass.table') into bqsr_second_pass_table_ch

    script:
    """
    mkdir tmp
    gatk BaseRecalibrator \
            -I $bam \
            -R $reference \
            -L $targets \
            --known-sites $dbsnp \
            -O 'bqsr_second_pass.table' \
            --tmp-dir=tmp \
            --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g"
    """
}

bqsr_table_analyze_ch.join(bqsr_second_pass_table_ch).set{data_analyze_covariates_ch}

// Evaluate BAM before and after recalibration, by comparing the BQSR tables of the first and second
// pass.
process analyze_covariates {
    publishDir "$outdir/$sample/bam/bqsr", mode: 'copy', overwrite: true

    input:
    set sample, file(bqsr_table), file(bqsr_table_second_pass) from data_analyze_covariates_ch

    output:
    file 'AnalyzeCovariates.pdf'

    script:
    """
    gatk AnalyzeCovariates \
        -before $bqsr_table \
        -after $bqsr_table_second_pass \
        -plots 'AnalyzeCovariates.pdf'
    """
}

/*
Call, annotate, and filter variants.
*/

// Call variants in sample with HapltypeCaller, yielding a GVCF.
process call_sample {
    publishDir "$outdir/$sample/gvcf", mode: 'copy', overwrite: true

    input:
    set sample, file(bam), file(bai) from recalibrated_bam_call_ch

    output:
    set sample, file("gvcf.g.vcf"), file("gvcf.g.vcf.idx") into gvcf_ch

    script:
    """
    mkdir tmp
    gatk HaplotypeCaller  \
        -I $bam \
        -O "gvcf.g.vcf" \
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
        --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g"
    """
}

// Genotype the GVCF in the previous process, yielding a VCF.
process genotyping {
    input:
    set sample, file(gvcf), file(idx) from gvcf_ch

    output:
    set sample, file("genotyped.vcf"), file("genotyped.vcf.idx") into genotyped_vcf_ch

    script:
    """
    mkdir tmp
    gatk GenotypeGVCFs \
        -V $gvcf \
        -R $reference \
        -O "genotyped.vcf" \
        --tmp-dir=tmp \
        --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g"
    """
}

// Add rsid from dbSNP
// NOTE: VariantAnnotator is still in beta (as of 20th of March 2019).
process annotate_rsid {
    input:
    set sample, file(vcf), file(idx) from genotyped_vcf_ch

    output:
    set sample, file("rsid_ann.vcf"), file("rsid_ann.vcf.idx") into rsid_annotated_vcf_ch

    script:
    """
    gatk VariantAnnotator \
        -R $reference \
        -V $vcf \
        --dbsnp $dbsnp \
        -O "rsid_ann.vcf" \
        --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g"
    """
}

// Filter variants, adding various filter tags to the "FILTER" field of the VCF.
// FIXME: I'm getting warnings that MQRankSum and ReadPosRankSum don't exist.
process filter_variants {
    input:
    set sample, file(vcf), file(idx) from rsid_annotated_vcf_ch

    output:
    set sample, file("filtered.vcf"), file("filtered.vcf.idx") into filtered_vcf_ch

    script:
    """
    gatk VariantFiltration \
        -V $vcf \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O "filtered.vcf" \
        --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g"
    """
}

// Annotate the VCF with effect prediction. Output some summary stats from the effect prediction as well.
process annotate_effect {
    publishDir "$outdir/$sample/vcf", pattern: "snpEff_stats.csv", mode: 'copy', overwrite: true

    input:
    set sample, file(vcf), file(idx) from filtered_vcf_ch

    output:
    set sample, file("effect_annotated.vcf") into variants_phase_ch
    file "snpEff_stats.csv"

    script:
    """
    snpEff -Xmx${task.memory.toGiga()}g \
         -i vcf \
         -o vcf \
         -csvStats "snpEff_stats.csv" \
         hg38 \
         -v \
         $vcf > "effect_annotated.vcf"
    """
}

/*
Phase the haplotypes in the VCF using HapCUT2, and then phase the BAM using WhatsHap.
*/

// Prepare input data for processes.
variants_phase_ch.join(bam_phase_vcf_ch).set { data_extract_hairs_ch }


// Convert BAM file to the compact fragment file format containing only haplotype-relevant information.
process extract_hairs {
    input:
    set sample, file(vcf), file(bam), file(bai) from data_extract_hairs_ch

    output:
    set sample, file(vcf), file(bam), file(bai), file("unlinked_fragment") into data_link_fragments_ch

    script:
    """
    extractHAIRS --10X 1 --bam $bam --VCF $vcf --out "unlinked_fragment"
    """
}

// Use LinkFragments to link fragments into barcoded molecules.
process link_fragments {
    input:
    set sample, file(vcf), file(bam), file(bai), file(unlinked_fragments) from data_link_fragments_ch

    output:
    set sample, file(vcf), file(bam), file(bai), file('linked_fragments') into data_phase_vcf_ch

    script:
    """
    LinkFragments.py --bam $bam --VCF $vcf --fragments $unlinked_fragments --out "linked_fragments"
    """
}

// Use HAPCUT2 to assemble fragment file into haplotype blocks.
process phase_vcf {
    input:
    set sample, file(vcf), file(bam), file(bai), file(linked_fragments) from data_phase_vcf_ch

    output:
    set sample, file("haplotypes.phased.VCF") into phased_vcf_ch

    script:
    """
    HAPCUT2 --nf 1 --outvcf 1 --fragments $linked_fragments --VCF $vcf --output "haplotypes"
    """
}

// Compress and index the phased VCF.
process zip_and_index_vcf {
    publishDir "$outdir/$sample/vcf", mode: 'copy', pattern: '*.vcf.gz', overwrite: true,
        saveAs: { filename -> "${sample}.vcf.gz" }
    publishDir "$outdir/$sample/vcf", mode: 'copy', pattern: '*.vcf.gz.tbi', overwrite: true,
        saveAs: { filename -> "${sample}.vcf.gz.tbi" }

    input:
    set sample, file(vcf) from phased_vcf_ch

    output:
    set sample, file("variants.vcf.gz"), file("variants.vcf.gz.tbi") into variants_phase_bam_ch, variants_evaluate_ch, variants_phasing_stats_ch

    script:
    """
    cat $vcf | bgzip -c > "variants.vcf.gz"
    tabix "variants.vcf.gz"
    """
}

variants_phase_bam_ch.join(bam_phase_bam_ch).set{data_haplotag_bam_ch}

// Add haplotype information to BAM, tagging each read with a haplotype (when possible), using
// the haplotype information from the phased VCF.
process haplotag_bam {
    input:
    set sample, file(vcf), file(idx), file(bam), file(bai) from data_haplotag_bam_ch

    output:
    set sample, file("phased.bam") into phased_bam_ch

    script:
    """
    whatshap haplotag --ignore-read-groups --reference $reference -o "phased.bam" $vcf $bam
    """
}

process index_phased_bam {
    publishDir "$outdir/$sample/bam", mode: 'copy', pattern: '*.bam', overwrite: true,
        saveAs: { filename -> "${sample}.bam" }
    publishDir "$outdir/$sample/bam", mode: 'copy', pattern: '*.bam.bai', overwrite: true,
        saveAs: { filename -> "${sample}.bam.bai" }

    input:
    set sample, file(bam) from phased_bam_ch

    output:
    set sample, file("phased.bam"), file("phased.bam.bai") into indexed_phased_bam_qualimap_ch, indexed_phased_bam_bx_ch

    script:
    """
    samtools index "phased.bam"
    """
}

/*
Below we perform QC of data.
*/

// The GATK variant evaluation module counts variants stratified w.r.t. filters, compares
// overlap with DBSNP, and more.
process variant_evaluation {
    publishDir "$outdir/$sample/vcf", mode: 'copy', overwrite: true

    input:
    set sample, file(vcf), file(idx) from variants_evaluate_ch

    output:
    file "variant_eval.table"
    set sample, val('done') into variants_status_ch

    script:
    """
    # VariantEval fails if the output file doesn't already exist. NOTE: this should be fixed in a newer version of GATK, as of the 19th of February 2019.
    echo -n > "variant_eval.table"

    gatk VariantEval \
        -R $reference \
        --eval $vcf \
        --output "variant_eval.table" \
        --dbsnp $dbsnp \
        -L $targets \
        -no-ev -no-st \
        --eval-module TiTvVariantEvaluator \
        --eval-module CountVariants \
        --eval-module CompOverlap \
        --eval-module ValidationReport \
        --stratification-module Filter
    """
}

// Run Qualimap for QC metrics of recalibrated BAM.
process qualimap_analysis {
    publishDir "$outdir/$sample/bam", mode: 'copy', overwrite: true

    input:
    set sample, file(bam), file(bai) from indexed_phased_bam_qualimap_ch

    output:
    file "qualimap_results"
    set sample, val('done') into qualimap_status_ch

    script:
    """
    # Make sure QualiMap doesn't attemt to open a display server.
    unset DISPLAY
    # Run QualiMap.
    qualimap bamqc \
        -gd HUMAN \
        -bam $bam \
        -gff $targets \
        -outdir "qualimap_results" \
        --skip-duplicated \
        --collect-overlap-pairs \
        -nt ${task.cpus} \
        --java-mem-size=${task.memory.toGiga()}G
    """
}

// Get basic statistics about haplotype phasing blocks.
// NOTE: provide a list of reference chromosome sizes to get N50.
// NOTE: this does not produce a MultiQC report, and neither does bx_stats. If I do implement
// that, there must be a status value emitted from these, to create dependency.
process phasing_stats {
    publishDir "$outdir/$sample/vcf/phasing", mode: 'copy', overwrite: true

    input:
    set sample, file(vcf), file(idx) from variants_phasing_stats_ch

    output:
    file 'phase_blocks.gtf'
    file 'phasing_stats.tsv'

    script:
    """
    whatshap stats $vcf --gtf phase_blocks.gtf --tsv phasing_stats.tsv
    """
}

// Get basic statistics about linked-read barcodes in the BAM.
process bx_stats {
    publishDir "$outdir/$sample/bam", mode: 'copy', overwrite: true

    input:
    set sample, file(bam), file(bai) from indexed_phased_bam_bx_ch

    output:
    file 'bx_stats.csv'
    file 'bx_summary.txt'

    script:
    """
    echo -e 'BX\tread count\tmedian insert size\tmedian mapq\tmedian AS' > bx_stats.csv
    bxtools stats $bam >> bx_stats.csv
    bx_summary.py bx_stats.csv > bx_summary.txt
    """
}

// FIXME: it seems like this process is imported as cached even though the qualimap process is not.
// This means that the multiqc output folder is not updated.
process multiqc {
    publishDir "$outdir/multiqc", mode: 'copy', overwrite: true

    input:
    val qstatus from qualimap_status_ch
    val vstatus from variants_status_ch

    output:
    file "multiqc_report.html" into multiqc_report_ch
    file "multiqc_data" into multiqc_data_ch

    script:
    """
    multiqc -f $outdir
    """
}

