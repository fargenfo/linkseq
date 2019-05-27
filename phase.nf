#!/usr/bin/env nextflow
/*
Author: Ólavur Mortensen <olavur@fargen.fo>
*/

// Input parameters.
//params.sample = null
params.bam_paths = null
params.vcf = null
params.reference = null
params.targets = null
params.threads = null
params.mem = null
params.HapCUT2 = null
params.interval = null
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
//assert params.sample != null, 'Input parameter "sample" cannot be unasigned.'
//assert params.bam != null, 'Input parameter "bam" cannot be unasigned.'
assert params.bam_paths != null, 'Input parameter "bam_paths" cannot be unasigned.'
assert params.vcf != null, 'Input parameter "vcf" cannot be unasigned.'
assert params.HapCUT2!= null, 'Input parameter "HapCUT2" cannot be unasigned.'
assert params.interval != null, 'Input parameter "interval" cannot be unasigned.'
assert params.reference != null, 'Input parameter "reference" cannot be unasigned.'
assert params.targets != null, 'Input parameter "targets" cannot be unasigned.'
assert params.threads != null, 'Input parameter "threads" cannot be unasigned.'
assert params.mem != null, 'Input parameter "mem" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

println "P I P E L I N E     I P U T S    "
println "================================="
//println "sample             : ${params.sample}"
//println "bam                : ${params.bam}"
println "bam_paths          : ${params.bam_paths}"
println "vcf                : ${params.vcf}"
println "HapCUT2            : ${params.HapCUT2}"
println "interval           : ${params.interval}"
println "reference          : ${params.reference}"
println "targets            : ${params.targets}"
println "threads            : ${params.threads}"
println "mem                : ${params.mem}"
println "outdir             : ${params.outdir}"

// Turn the file with FASTQ paths into a channel with [sample, path] tuples.
bam_paths_ch = Channel.fromPath(params.bam_paths)
bam_paths_ch = bam_paths_ch.splitCsv(header: true).map { it -> [it.sample, file(it.bam_path), file(it.bai_path)] }

// Get file handlers for input files.
//bam = file(params.bam)
vcf = file(params.vcf)
reference = file(params.reference)  // Directory of 10x reference.
reference_fa = file(params.reference + '/fasta/genome.fa')  // Reference fasta file.
targets = file(params.targets)

HAPCUT2 = params.HapCUT2 + "/build/HAPCUT2"
extractHAIRS = params.HapCUT2 + "/build/extractHAIRS"
LinkFragments = params.HapCUT2 + "/utilities/LinkFragments.py"

process make_small_bam {
    input:
    set sample, file(bam), file(bai) from bam_paths_ch

    output:
    //set file("small.bam"), file("small.bam.bai") into small_bam_genotyping_ch, small_bam_extract_ch, small_bam_link_ch, small_bam_phase_bam_ch
    set sample, file("small.bam"), file("small.bam.bai") into small_bam_ch, small_bam_haplotag_bam_ch

    script:
    """
    samtools view -b -o "small.bam" -T $reference_fa $bam $params.interval
    samtools index -b "small.bam"
    """
}

process get_sample_names {
    output:
    file "samples.txt" into sample_names_ch

    script:
    """
    # Get the line of the VCF containing the sample names.
    grep -m 1 "^#CHROM" $vcf > "header_line.txt"
    # Get the sample names in one line.
    cat "header_line.txt" | awk '{print substr(\$0, index(\$0,\$10))}' > "samples_line.txt"
    # Make a file with one line per sample name.
    cat "samples_line.txt" | tr -s '\t' '\n' > "samples.txt"
    """
}

// Convert this channel from a file with sample names, to a channel with one item per sample name.
sample_names_ch = sample_names_ch.splitText().map { it.trim() }

// FIXME: -L option only for testing
process get_sample_vcf {
    input:
    val sample from sample_names_ch

    output:
    //set sample, file("sample.vcf") into vcf_extract_ch, vcf_link_ch, vcf_phase_ch
    set sample, file("sample.vcf"), file("sample.vcf.idx") into vcf_ch

    script:
    """
    mkdir tmp
    gatk SelectVariants \
        -V $vcf \
        -R $reference_fa \
        -sn $sample \
        -L ${params.interval} \
        -O "sample.vcf" \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    """
}

// Make a channel with (sample ID, VCF, VCF index, BAM, BAM index) tuples.
vcf_ch.join(small_bam_ch).into { data_extract_ch; data_link_ch; data_phase_ch }

process extract_hairs {
    input:
    //set sample, file(vcf) from vcf_extract_ch
    //set file(bam), file(bai) from small_bam_extract_ch
    set sample, file(vcf), file(idx), file(bam), file(bai) from data_extract_ch

    output:
    //file "unlinked_fragment" into unlinked_fragments_ch
    set sample, file("unlinked_fragment") into unlinked_fragments_ch

    script:
    """
    $extractHAIRS --10X 1 --region $params.interval --bam $bam --VCF $vcf --out "unlinked_fragment"
    """
}

// Add unlinked fragment file to tuple for next process.
data_link_ch = data_link_ch.join(unlinked_fragments_ch)

process link_fragments {
    input:
    //file vcf from vcf_link_ch
    //set file(bam), file(bai) from small_bam_link_ch
    //file unlinked_fragments from unlinked_fragments_ch
    set sample, file(vcf), file(idx), file(bam), file(bai), file(unlinked_fragments) from data_link_ch

    output:
    set sample, file("linked_fragments") into linked_fragments_ch

    script:
    """
    python3 $LinkFragments --bam $bam --VCF $vcf --fragments $unlinked_fragments --out "linked_fragments"
    """
}

// Add linked fragment file to tuple for next process.
data_phase_ch = data_phase_ch.join(linked_fragments_ch)

process phase_vcf {
    input:
    //file vcf from vcf_phase_ch
    //file linked_fragments from linked_fragments_ch
    set sample, file(vcf), file(idx), file(bam), file(bai), file(linked_fragments) from data_phase_ch

    output:
    file "haplotypes" into haplotypes_ch
    //set sample, file("haplotypes.phased.VCF") into phased_vcf_get_header_ch, phased_vcf_reheader_ch
    set sample, file("haplotypes.phased.VCF") into phased_vcf_ch

    script:
    """
    $HAPCUT2 --nf 1 --outvcf 1 --fragments $linked_fragments --VCF $vcf --output "haplotypes"
    """
}

process index_and_zip_vcf {
    input:
    //set sample, file(vcf) from vcf_reheaded_ch
    set sample, file(vcf) from phased_vcf_ch

    output:
    set sample, file("phased.vcf.gz"), file("phased.vcf.gz.tbi") into phased_vcf_haplotag_ch
    file "phased.vcf.gz" into phased_vcf_merge_ch

    script:
    """
    bgzip -c $vcf > "phased.vcf.gz"
    tabix "phased.vcf.gz"
    """
}


process merge_phased_vcf {
    input:
    val vcf_list from phased_vcf_merge_ch.toList()

    output:
    set file("phased_merged.vcf.gz"), file("phased_merged.vcf.gz.tbi") into phased_merged_ch

    script:
    // Prepare input string for vcf-merge.
    vcf_list_str = (vcf_list as List)
        .join(' ')  // Join paths in single string.
    """
    #vcf-merge $vcf_list_str > "phased_merged.vcf"
    vcf-merge $vcf_list_str | bgzip -c > "phased_merged.vcf.gz"
    tabix "phased_merged.vcf.gz"
    """
}


//process merge_phased_vcf {
//    echo true
//    input:
//    val vcf_list from phased_vcf_merge_ch.toList()
//
//    script:
//    // Prepare input string for MergeVcfs.
//    vcf_list_str = (vcf_list as List)
//        .collect({ "-I " + it })  // Prepend path with "-I".
//        .join(' ')  // Join paths in single string.
//    """
//    gatk MergeVcfs \
//        $vcf_list_str \
//        -O phased_multisample.vcf \
//        --TMP_DIR=tmp \
//        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
//    """
//}




// Take the compressed and indexed VCF from above and join by sample name with BAM files.
data_haplotag_ch = phased_vcf_haplotag_ch.join(small_bam_haplotag_bam_ch)

process haplotag_bam {
    input:
    //set file(phased_vcf), file(phased_vcf_idx) from phased_vcf_idx_gz_ch
    //set file(bam), file(bai) from small_bam_phase_bam_ch
    set sample, file(phased_vcf), file(idx), file(bam), file(bai) from data_haplotag_ch

    output:
    file "phased.bam" into phased_bam_ch

    script:
    """
    whatshap haplotag --ignore-read-groups --reference $reference_fa -o "phased.bam" $phased_vcf $bam
    """
}

// NOTE: When we have a VCF phased by HapCUT2, we can apply this haplotype information to
// a separate VCF using WhatsHap. I do not know how much loss of information there is in
// this process.
//process whatshap_phase_vcf {
//    input:
//    file vcf from ???
//    file phased_vcf from ???
//
//    output:
//
//    script:
//    """
//    whatshap phase --ignore-read-groups --reference $reference_fa --indels -o $vcf $phased_vcf
//    """
//}

//process genotyping {
//    input:
//    set file(bam), file(bai) from small_bam_genotyping_ch
//
//    output:
//    file "genotyped.vcf" into vcf_extract_ch, vcf_link_ch, vcf_phase_ch
//
//    script:
//    """
//    mkdir tmp
//    gatk GenotypeGVCFs \
//        -V $gvcf \
//        -R $reference_fa \
//        --intervals $params.interval \
//        -O "genotyped.vcf" \
//        --tmp-dir=tmp \
//        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
//    """
//}

// NOTE: I can make a very simply VCF to phase with mpileup
// use either -T or -r to restrict analysis to taget regions.
// -x could make this faster
//process call_simple_vcf {
//    output:
//    file "simple.vcf" into vcf_ch
//
//    script:
//    """
//    bcftools mpileup -f $reference_fa -o "simple.vcf" -O "v" $bam
//    """
//}
