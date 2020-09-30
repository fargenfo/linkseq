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
params.sample = null
params.fastq_r1 = null
params.fastq_r2 = null
params.reference = null
params.targets = null
params.whitelist = null
params.bcbins = null
params.dbsnp = null
params.snpeff_datadir = null
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
assert params.reference != null, 'Input parameter "reference" cannot be unasigned.'
assert params.targets != null, 'Input parameter "targets" cannot be unasigned.'
assert params.whitelist != null, 'Input parameter "whitelist" cannot be unasigned.'
assert params.bcbins != null, 'Input parameter "bcbins" cannot be unasigned.'
assert params.dbsnp != null, 'Input parameter "dbsnp" cannot be unasigned.'
assert params.snpeff_datadir != null, 'Input parameter "snpeff_datadir" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

if(params.fastq_csv == null & params.sample == null) {
    assert false, "Either --fastq_csv must be provided, or --sample, --fastq_r1 and --fastq_r2 must be provided."
}

println "L I N K S E Q    "
println "================================="
println "reference          : ${params.reference}"
println "targets            : ${params.targets}"
println "whitelist          : ${params.whitelist}"
println "bcbins             : ${params.bcbins}"
println "dbsnp              : ${params.dbsnp}"
println "outdir             : ${params.outdir}"
println "================================="
println "Command line        : ${workflow.commandLine}"
println "Profile             : ${workflow.profile}"
println "Project dir         : ${workflow.projectDir}"
println "Launch dir          : ${workflow.launchDir}"
println "Work dir            : ${workflow.workDir}"
println "Container engine    : ${workflow.containerEngine}"
println "================================="
println "Project             : $workflow.projectDir"
println "Git info            : $workflow.repository - $workflow.revision [$workflow.commitId]"
println "Cmd line            : $workflow.commandLine"
println "Manifest version    : $workflow.manifest.version"
println "================================="

// Get file handlers for input files.
reference = file(params.reference, checkIfExists: true)
targets = file(params.targets, checkIfExists: true)
whitelist = file(params.whitelist, checkIfExists: true)
dbsnp = file(params.dbsnp, checkIfExists: true)
snpeff_datadir = file(params.snpeff_datadir, checkIfExists: true)
outdir = file(params.outdir)

if(params.fastq_csv != null) {
    // Read FASTQ read 1 and 2, as well as sample IDs, from input CSV file.
    Channel.fromPath(params.fastq_csv)
        .splitCsv(header:true)
        .map{ row-> tuple(row.sample, file(row.read1), file(row.read2)) }
        .into { fastq_ch; fastq_print_ch }
} else {
    // Put sample name and FASTQ read 1 and 2 in a channel.
    Channel.of( tuple(params.sample, file(params.fastq_r1), file(params.fastq_r2)) )
        .into { fastq_ch; fastq_print_ch }
}

// Print the samples to be processed.
println 'Sample ID\tFASTQ read 1 files\tFASTQ read 2 files'
fastq_print_ch.subscribe onNext: { row ->
    sample = row[0]
    // Get the filenames from the file objects.
    read1 = row[1].name
    read2 = row[2].name
    // read1 and read2 may be either a single string or a list of strings.
    if(read1 instanceof List) {
        read1 = read1.join(',')
        read2 = read2.join(',')
    }
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
    set sample, file('interleaved.fastq.gz') into fastq_count_ch, fastq_preproc_ch, fastq_readgroup_ch, fastq_check_sync_ch, fastq_qc_ch

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
// First group EMA bins by key (sample), yielding a channel with (sample, bin tuple) records.
bins_ema_ch.groupTuple().set{grouped_bins_ema_ch}
// Join this channel with the readgroups, obtaining (sample, reagroup, bin tuple) records.
readgroup_ema_ch.join(grouped_bins_ema_ch).set{data_ema_ch}
// Transpose this channel (w.r.t. the tuples), obtaining (sample, reagroup, bin) channel.
data_ema_ch.transpose().set{data_ema_ch}

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
    publishDir "$outdir/multiqc_logs/bqsr_before", mode: 'copy', overwrite: true,
        saveAs: { filename -> "$sample" }

    input:
    set sample, file(bam), file(bai) from indexed_bam_prepare_ch

    output:
    set sample, file('bqsr.table') into bqsr_table_analyze_ch, bqsr_table_apply_ch, bqsr_table_multiqc_ch

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
    publishDir "$outdir/multiqc_logs/bqsr_after", mode: 'copy', overwrite: true,
        saveAs: { filename -> "$sample" }

    input:
    set sample, file(bam), file(bai) from recalibrated_bam_second_pass_ch

    output:
    set sample, file('bqsr_second_pass.table') into bqsr_second_pass_table_ch, bqsr_second_pass_multiqc_ch

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
    publishDir "$outdir/multiqc_logs/AnalyzeCovariates", mode: 'copy', overwrite: true,
        saveAs: { filename -> "$sample" }

    input:
    set sample, file(bqsr_table), file(bqsr_table_second_pass) from data_analyze_covariates_ch

    output:
    set sample, file('AnalyzeCovariates.pdf') into analyze_covariates_multiqc_ch

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
    set sample, file("rsid_ann.vcf"), file("rsid_ann.vcf.idx") into rsid_annotated_vcf_snp_ch, rsid_annotated_vcf_indel_ch

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


// Splitting VCF into SNPs and indels, because they have to be filtered seperately
process subset_snps {

    input:
    set sample, file(vcf), file(idx) from rsid_annotated_vcf_snp_ch

    output:
    set sample, file("snp.vcf"), file("snp.vcf.idx") into snpsubset_filter_ch

    script:
    """
    gatk SelectVariants \
    -V $vcf \
    -select-type SNP \
    -O "snp.vcf" \
    --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g"
    """
}

process subset_indels {

    input:
    set sample, file(vcf), file(idx) from rsid_annotated_vcf_indel_ch

    output:
    set sample, file("indel.vcf"), file("indel.vcf.idx") into indelsubset_filter_ch

    script:
    """
    gatk SelectVariants \
    -V $vcf \
    -select-type INDEL \
    -O "indel.vcf" \
    --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g"
    """
}

// Hard filter variants, adding various filter tags to the "FILTER" field of the VCF.
// FIXME: I'm getting warnings that MQRankSum and ReadPosRankSum don't exist.
process hard_filter_snps {
    input:
    set sample, file(vcf), file(idx) from snpsubset_filter_ch

    output:
    set sample, file("filtered.vcf"), file("filtered.vcf.idx") into filtered_snp_vcf_ch

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

process hard_filter_indels {

    input:
    set sample, file(vcf), file(idx) from indelsubset_filter_ch

    output:
    set sample, file("filtered_indel.vcf"), file("filtered_indel.vcf.idx") into filtered_indel_vcf_ch

    script:
    """
    gatk VariantFiltration \
    -V $vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O "filtered_indel.vcf" \
    --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g"
    """
 }

// Merge the SNP and INDEL vcfs before continuing.
// There will only be one joint vcf with all samples, so we don't need to set sample.
process join_snps_indels {

    input:
    set sample, file(vcf_snp), file(idx_snp) from filtered_snp_vcf_ch
    set sample, file(vcf_indel), file(idx_indel) from filtered_indel_vcf_ch

    output:
    set sample, file("joined_snp_indel.vcf"), file("joined_snp_indel.vcf.idx") into joined_snp_indel_vcf_ch

    script:
    """
    gatk MergeVcfs \
    -I $vcf_snp \
    -I $vcf_indel \
    -O "joined_snp_indel.vcf" \
    --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g"
    """
}

// Annotate the VCF with effect prediction. Output some summary stats from the effect prediction as well.
// We tell SnpEff not to attempt to download the reference data, and supply the reference data directory
// path explicitly instead. Otherwise, SnpEff will download these data for every new environment and for
// every new container.
process annotate_effect {
    publishDir "$outdir/multiqc_logs/SnpEff", pattern: "snpEff_stats.csv", mode: 'copy', overwrite: true,
        saveAs: { filename -> "$sample" }

    input:
    set sample, file(vcf), file(idx) from joined_snp_indel_vcf_ch

    output:
    set sample, file("effect_annotated.vcf") into variants_phase_ch
    set sample, file("snpEff_stats.csv") into snpeff_multiqc_ch

    script:
    """
    snpEff -Xmx${task.memory.toGiga()}g \
         -i vcf \
         -o vcf \
         -csvStats "snpEff_stats.csv" \
         -nodownload \
         -dataDir $snpeff_datadir \
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

// Compress the phased VCF.
process zip_vcf {
    publishDir "$outdir/$sample/vcf", mode: 'copy', pattern: '*.vcf.gz', overwrite: true,
        saveAs: { filename -> "${sample}.vcf.gz" }

    input:
    set sample, file(vcf) from phased_vcf_ch

    output:
    set sample, file("variants.vcf.gz") into variants_compressed_index_ch

    script:
    """
    cat $vcf | bgzip -c > "variants.vcf.gz"
    """
}

// Index the phased VCF.
process index_vcf {
    publishDir "$outdir/$sample/vcf", mode: 'copy', pattern: '*.vcf.gz.tbi', overwrite: true,
        saveAs: { filename -> "${sample}.vcf.gz.tbi" }

    input:
    set sample, file(vcf) from variants_compressed_index_ch

    output:
    set sample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into variants_phase_bam_ch, variants_evaluate_ch, variants_phasing_stats_ch

    script:
    """
    # Rename VCF, so that we may put both the VCF and the index into the output channel.
    mv $vcf indexed.vcf.gz
    tabix indexed.vcf.gz
    """
}

variants_phase_bam_ch.join(bam_phase_bam_ch).set{data_haplotag_bam_ch}

// Add haplotype information to BAM, tagging each read with a haplotype (when possible), using
// the haplotype information from the phased VCF.
process haplotag_bam {
    input:
    set sample, file(vcf), file(idx), file(bam), file(bai) from data_haplotag_bam_ch

    output:
    set sample, file("${sample}.bam") into phased_bam_ch

    script:
    """
    whatshap haplotag --ignore-read-groups --reference $reference -o "${sample}.bam" $vcf $bam
    """
}

process index_phased_bam {
    publishDir "$outdir/$sample/bam", mode: 'copy', pattern: '*.bam', overwrite: true
    publishDir "$outdir/$sample/bam", mode: 'copy', pattern: '*.bam.bai', overwrite: true

    input:
    set sample, file(bam) from phased_bam_ch

    output:
    set sample, file(bam), file("${bam}.bai") into indexed_phased_bam_qualimap_ch, indexed_phased_bam_bx_ch

    script:
    """
    samtools index $bam
    """
}

/*
Below we perform QC of data.
*/

// Do QG of interleaved FASTQ.
// While having stats on lanes and reads would be beneficial, doing QC of the merged interleaved FASTQ
// makes the MultiQC easier to view, as the FastQC report is associated with one sample. When we combine
// MultiQC reports of several samples, we can easily compare samples.
// NOTE: FastQC claims 250 MB of memory for every thread that is allocated to it.
// FIXME: control that the files are published in such a way that the MultiQC report has data
// with the correct sample name.
process fastqc_analysis {
    memory { 250.MB * task.cpus }

    publishDir "$outdir/multiqc_logs/fastqc", mode: 'copy', pattern: '*.zip', overwrite: true

	input:
	set sample, file(fastq) from fastq_qc_ch

    output:
    set sample, file('*.zip') into fastqc_multiqc_ch

    script:
    """
    # We unset the DISPLAY variable to avoid having FastQC try to open the GUI.
    unset DISPLAY
    mkdir tmp
    # Rename the FASTQ such that the sample name is correctly displayed in MultiQC
    mv $fastq ${sample}.fastq.gz
    fastqc -q --dir tmp --threads ${task.cpus} --outdir . ${sample}.fastq.gz
    """
}

// The GATK variant evaluation module counts variants stratified w.r.t. filters, compares
// overlap with DBSNP, and more.
process variant_evaluation {
    publishDir "$outdir/multiqc_logs/VariantEval", mode: 'copy', overwrite: true,
        saveAs: { filename -> "$sample" }

    input:
    set sample, file(vcf), file(idx) from variants_evaluate_ch

    output:
    set sample, file("variant_eval.table") into varianteval_multiqc_ch

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
    publishDir "$outdir/multiqc_logs/qualimap", mode: 'copy', overwrite: true,
        saveAs: { filename -> "$sample" }

    input:
    set sample, file(bam), file(bai) from indexed_phased_bam_qualimap_ch

    output:
    set sample, file("qualimap_results") into qualimap_multiqc_ch

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
// NOTE: this does not produce a MultiQC report, and neither does bx_stats.
process phasing_stats {
    publishDir "$outdir/$sample/vcf/phasing/", pattern: "*.gtf", mode: 'copy', overwrite: true
    publishDir "$outdir/multiqc_logs/WhatsHap", pattern: "*.tsv", mode: 'copy', overwrite: true,
        saveAs: { filename -> "$sample" }

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
// NOTE: this does not produce a MultiQC report.
process bx_stats {
    publishDir "$outdir/$sample/bam", pattern: "*.csv", mode: 'copy', overwrite: true
    publishDir "$outdir/multiqc_logs/bx_stats", pattern: ".txt", mode: 'copy', overwrite: true,
        saveAs: { filename -> "$sample" }

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

// To ensure MultiQC is run only when all reports have been created, and to ensure that a new MultiQC report is created if
// any file is created, we need to create a dependency between the QC processes and Multiqc.

// Get all QC reports in one channel.
multiqc_ch = bqsr_table_multiqc_ch.mix(  bqsr_second_pass_multiqc_ch, analyze_covariates_multiqc_ch, snpeff_multiqc_ch, fastqc_multiqc_ch, varianteval_multiqc_ch, qualimap_multiqc_ch )

process multiqc {
    publishDir "$outdir/multiqc", mode: 'copy', overwrite: true

    input:
    // This process depends on all the QC processes.
    val temp from multiqc_ch.collect()

    output:
    file "multiqc_report.html" into multiqc_report_ch
    file "multiqc_data" into multiqc_data_ch

    script:
    """
    multiqc -f $outdir/multiqc_logs
    """
}


workflow.onComplete {
    log.info "L I N K S E Q   "
    log.info "================================="
    log.info "reference          : ${params.reference}"
    log.info "targets            : ${params.targets}"
    log.info "whitelist          : ${params.whitelist}"
    log.info "bcbins             : ${params.bcbins}"
    log.info "dbsnp              : ${params.dbsnp}"
    log.info "outdir             : ${params.outdir}"
    log.info "================================="
    log.info "Command line        : ${workflow.commandLine}"
    log.info "Profile             : ${workflow.profile}"
    log.info "Project dir         : ${workflow.projectDir}"
    log.info "Launch dir          : ${workflow.launchDir}"
    log.info "Work dir            : ${workflow.workDir}"
    log.info "Container engine    : ${workflow.containerEngine}"
    log.info "================================="
    log.info "Project             : $workflow.projectDir"
    log.info "Git info            : $workflow.repository - $workflow.revision [$workflow.commitId]"
    log.info "Cmd line            : $workflow.commandLine"
    log.info "Manifest version    : $workflow.manifest.version"
    log.info "================================="
    log.info "Completed at        : ${workflow.complete}"
    log.info "Duration            : ${workflow.duration}"
    log.info "Success             : ${workflow.success}"
    log.info "Exit status         : ${workflow.exitStatus}"
    log.info "Error report        : ${(workflow.errorReport ?: '-')}"
    log.info "================================="
}

