#!/usr/bin/env nextflow

// Input parameters.
params.gvcf_path = null
params.genomicsdb = null
params.reference = null
params.dbsnp = null
params.targets = null
params.mills = null
params.kGphase1 = null
params.kGphase3 = null
params.omni = null
params.hapmap = null
params.threads = null
params.mem = null
params.outdir = null
params.help = false

# Resources. Reference sequence and target regions.
mills=$resources_dir/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
kG=$resources_dir/gatk_bundle/1000G_phase1.snps.high_confidence.hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
omni=$resources_dir/gatk_bundle/1000G_omni2.5.hg38/1000G_omni2.5.hg38.vcf.gz
hapmap=$resources_dir/gatk_bundle/hapmap_3.3.hg38.vcf.gz/hapmap_3.3.hg38.vcf.gz

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
assert params.gvcf_path != null, 'Input parameter "gvcf_path" cannot be unasigned.'
assert params.genomicsdb != null, 'Input parameter "genomicsdb" cannot be unasigned.'
assert params.reference != null, 'Input parameter "reference" cannot be unasigned.'
assert params.dbsnp != null, 'Input parameter "dbsnp" cannot be unasigned.'
assert params.mills != null, 'Input parameter "mills" cannot be unasigned.'
assert params.kGphase1 != null, 'Input parameter "kGphase1" cannot be unasigned.'
assert params.kGphase3 != null, 'Input parameter "kGphase3" cannot be unasigned.'
assert params.omni != null, 'Input parameter "omni" cannot be unasigned.'
assert params.hapmap != null, 'Input parameter "hapmap" cannot be unasigned.'
assert params.targets != null, 'Input parameter "targets" cannot be unasigned.'
assert params.threads != null, 'Input parameter "threads" cannot be unasigned.'
assert params.mem != null, 'Input parameter "mem" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

println "P I P E L I N E     I P U T S    "
println "================================="
println "gvcf_path          : ${params.gvcf_path}"
println "genomicsdb         : ${params.genomicsdb}"
println "reference          : ${params.reference}"
println "dbsnp              : ${params.dbsnp}"
println "mills              : ${params.mills}"
println "kGphase1           : ${params.kGphase1}"
println "kGphase3           : ${params.kGphase3}"
println "omni               : ${params.omni}"
println "hapmap             : ${params.hapmap}"
println "targets            : ${params.targets}"
println "threads            : ${params.threads}"
println "mem                : ${params.mem}"
println "outdir             : ${params.outdir}"

// Get file handlers for input files.
genomicsdb = file(params.genomicsdb)
reference = file(params.reference)  // Directory of 10x reference.
reference_fa = file(params.reference + '/fasta/genome.fa')  // Reference fasta file.
dbsnp = file(params.dbsnp)
mill = file(params.mills)
kGphase1 = file(params.kGphase1)
kGphase3 = file(params.kGphase3)
omni = file(params.omni)
hapmap = file(params.hapmap)
targets = file(params.targets)

gvcf_paths_ch = Channel.fromPath(params.gvcf_path + "/*.gvcf")
gvcf_arg_ch = gvcf_paths_ch.map { "-V " + it }

process consolidate_gvcf {
    echo true

    input:
    val gvcf_arg from gvcf_arg_ch.toList()

    output:
    val "done" into consolidate_gvcf_status_ch

    script:
    gvcf_arg_str = (gvcf_arg as List).join(' ')
    """
    mkdir tmp
    export TILEDB_DISABLE_FILE_LOCKING=1
    echo gatk GenomicsDBImport \
        $gvcf_arg_str \
        -L $targets \
        --genomicsdb-workspace-path $genomicsdb \
        --merge-input-intervals \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    """
}

process joint_genotyping {
    echo true

    input:
    val "done" from consolidate_gvcf_status_ch

    output:
    file "variants.vcf" into genotyped_vcf_ch
    file "variants.vcf" into genotyped_vcf_copy_ch

    script:
    """
    mkdir tmp
    echo gatk GenotypeGVCFs \
        -V gendb://$genomicsdb \
        -R $reference \
        -O "variants.vcf" \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    """
}


// TODO:
// Increasing --max-gaussians may work for larger sample sizes. For two samples, --max-gaussians=4 failed. exoseq uses 4.
// Filtering with VQSR based on DP is not recommended for exome data. Don't know if I'm currently doing this.



snps=snps.vcf
snps_recal=snps_recal.vcf
tranches=snps.tranches
recal=snps.recal
plots=snps.plots.R

process get_snps {
    input:
    file vcf from genotyped_vcf_ch

    output:
    file "snps.vcf" into snps_ch
    file "snps.vcf" into snps_copy_ch

    script:
    """
    gatk SelectVariants \
        -R $reference \
        -V $vcf \
        -O "snps.vcf" \
        --select-type-to-include SNP
    """
}

// TODO: do I want the plots from "snps.plots.R" (and "snps.plots.R.pdf")?
process recalibrate_snps {
    input:
    file vcf from snps_ch

    output:
    file "recal.table" into snps_recal_table_ch
    file "trances.table" into snps_trances_table_ch

    script:
    """
    gatk VariantRecalibrator \
        -R $reference \
        -V $vcf \
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
        -resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 $kGphase1 \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP \
        --max-gaussians 4 \
        -O "recal.table" \
        --tranches-file "tranches.table" \
        --rscript-file "snps.plots.R"
    """
}

process apply_vqsr_snps {
    input:
    file vcf from snps_copy_ch
    file recal_table from snps_recal_table_ch
    file trances_table from snps_trances_table_ch

    output:
    file "snps_recal.vcf" into recalibrated_snps_ch

    script:
    """
    gatk ApplyVQSR \
        -R $reference \
        -V $vcf \
        -O "snps_recal.vcf" \
        --truth-sensitivity-filter-level 99.0 \
        --tranches-file $tranches_table \
        --recal-file $recal_table \
        -mode SNP
    """
}

indels=indels.vcf
indels_recal=indels_recal.vcf
tranches=indels.tranches
recal=indels.recal
plots=indels.plots.R


process get_indels {
    input:
    file vcf from genotyped_vcf_copy_ch

    output:
    file "indels.vcf" into indels_ch
    file "indels.vcf" into indels_copy_ch

    script:
    """
    gatk SelectVariants \
        -R $reference \
        -V $vcf \
        -O "indels.vcf" \
        --select-type-to-include INDEL \
        --select-type-to-include MIXED \
        --select-type-to-include MNP \
        --select-type-to-include SYMBOLIC \
        --select-type-to-include NO_VARIATION
    """
}

process recalibrate_indels {
    input:
    file vcf from indels_ch

    output:
    file "recal.table" into indels_recal_table_ch
    file "trances.table" into indels_trances_table_ch

    script:
    """
    gatk VariantRecalibrator \
        -R $reference \
        -V $vcf \
        -resource:mills,known=false,training=true,truth=true,prior=12.0 $mills \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode INDEL \
        --max-gaussians 4 \
        -O "recal.table" \
        --tranches-file "tranches.table" \
        --rscript-file "plots.plots.R"
    """
}

process apply_vqsr_indels {
    input:
    file vcf from indels_copy_ch

    output:
    file "indels_recal.vcf" into recalibrated_indels_ch
    file recal_table from indels_recal_table_ch
    file trances_table from indels_trances_table_ch

    script:
    """
    gatk ApplyVQSR \
        -R $reference \
        -V vcf \
        -O "indels_recal.vcf" \
        --truth-sensitivity-filter-level 99.0 \
        --tranches-file $tranches_table \
        --recal-file $recal_table \
        -mode INDEL
    """
}

// FIXME: Merge variants some other way. CombineVariants may not be included in GATK4 in the near future.
//java -jar $gatk3 -T CombineVariants \
//    -R $reference \
//    -V:snps $snps_recal \
//    -V:indels $indels_recal \
//    -o recal.vcf \
//    --genotypemergeoption PRIORITIZE \
//    --rod_priority_list snps,indels

