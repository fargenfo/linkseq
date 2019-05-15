#!/usr/bin/env nextflow
/*
Author: Ólavur Mortensen <olavur@fargen.fo>
*/

// Input parameters.
params.sample = null
params.bam = null
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
assert params.sample != null, 'Input parameter "sample" cannot be unasigned.'
assert params.bam != null, 'Input parameter "bam" cannot be unasigned.'
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
println "sample             : ${params.sample}"
println "bam                : ${params.bam}"
println "vcf                : ${params.vcf}"
println "HapCUT2            : ${params.HapCUT2}"
println "interval           : ${params.interval}"
println "reference          : ${params.reference}"
println "targets            : ${params.targets}"
println "threads            : ${params.threads}"
println "mem                : ${params.mem}"
println "outdir             : ${params.outdir}"

// Get file handlers for input files.
bam = file(params.bam)
vcf = file(params.vcf)
reference = file(params.reference)  // Directory of 10x reference.
reference_fa = file(params.reference + '/fasta/genome.fa')  // Reference fasta file.
targets = file(params.targets)

HAPCUT2 = params.HapCUT2 + "/build/HAPCUT2"
extractHAIRS = params.HapCUT2 + "/build/extractHAIRS"
LinkFragments = params.HapCUT2 + "/utilities/LinkFragments.py"

process make_small_bam {
    output:
    set file("small.bam"), file("small.bam.bai") into small_bam_genotyping_ch, small_bam_extract_ch, small_bam_link_ch, small_bam_phase_bam_ch

    script:
    """
    samtools view -b -o "small.bam" -T $reference_fa $bam $params.interval
    samtools index -b "small.bam"
    """
}

// FIXME: -L option only for testing
process get_sample_vcf {
    output:
    file "sample.vcf" into vcf_extract_ch, vcf_link_ch, vcf_phase_ch

    script:
    """
    mkdir tmp
    gatk SelectVariants \
        -V $vcf \
        -R $reference_fa \
        -sn ${params.sample} \
        -L ${params.interval} \
        -O "sample.vcf" \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    """
}


process extract_hairs {
    input:
    file vcf from vcf_extract_ch
    set file(bam), file(bai) from small_bam_extract_ch

    output:
    file "unlinked_fragment" into unlinked_fragments_ch

    script:
    """
    $extractHAIRS --10X 1 --region $params.interval --bam $bam --VCF $vcf --out "unlinked_fragment"
    """
}


process link_fragments {
    input:
    file vcf from vcf_link_ch
    set file(bam), file(bai) from small_bam_link_ch
    file unlinked_fragments from unlinked_fragments_ch

    output:
    file "linked_fragments" into linked_fragments_ch

    script:
    """
    python3 $LinkFragments --bam $bam --VCF $vcf --fragments $unlinked_fragments --out "linked_fragments"
    """
}


process phase_vcf {
    input:
    file vcf from vcf_phase_ch
    file linked_fragments from linked_fragments_ch

    output:
    file "haplotypes" into haplotypes_ch
    file "haplotypes.phased.vcf" into phased_vcf_ch

    script:
    """
    $HAPCUT2 --nf 1 --outvcf 1 --fragments $linked_fragments --VCF $vcf --output "haplotypes"
    """
}

process index_and_zip_vcf {
    input:
    file vcf from phased_vcf_ch

    output:
    set file("phased.vcf.gz"), file("phased.vcf.gz.tbi") into phased_vcf_idx_gz_ch

    script:
    """
    bgzip -c $vcf > "phased.vcf.gz"
    tabix "phased.vcf.gz"
    """
}

process haplotag_bam {
    input:
    set file(phased_vcf), file(phased_vcf_idx) from phased_vcf_idx_gz_ch
    set file(bam), file(bai) from small_bam_phase_bam_ch

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
