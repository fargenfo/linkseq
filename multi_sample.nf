#!/usr/bin/env nextflow
/*
Author: Ólavur Mortensen <olavur@fargen.fo>
*/

// Input parameters.
params.gvcf_path = null
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

// Help message.
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
reference = file(params.reference)  // Directory of 10x reference.
reference_fa = file(params.reference + '/fasta/genome.fa')  // Reference fasta file.
dbsnp = file(params.dbsnp)
mills = file(params.mills)
kGphase1 = file(params.kGphase1)
// NOTE: 1000 Genomes phase 3 is not used anywhere at the moment, but could be used
// in CalculateGenotypePosteriors.
kGphase3 = file(params.kGphase3)
omni = file(params.omni)
hapmap = file(params.hapmap)
targets = file(params.targets)

// Get a list of GVCFs, and parse them to a format that GenomicsDBImport understands.
gvcf_paths = file(params.gvcf_path + "/*.g.vcf")
    .collect {"-V " + it}
    .join(' ')


// FIXME:
// remove this when done testing
targets = "chr17"
// FIXME



// Consolidate the GVCFs with a "genomicsdb" database, so that we are ready for joint genotyping.
process consolidate_gvcf {
    output:
    file "genomicsdb" into genomicsdb_ch

    script:
    """
    mkdir tmp
    export TILEDB_DISABLE_FILE_LOCKING=1
    gatk GenomicsDBImport \
        $gvcf_paths \
        -L $targets \
        --genomicsdb-workspace-path "genomicsdb" \
        --merge-input-intervals \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    """
}

// Multisample variant calling via GenotypeGVCFs.
process joint_genotyping {
    input:
    file genomicsdb from genomicsdb_ch

    output:
    set file("genotyped.vcf"), file("genotyped.vcf.idx") into genotyped_snprecal_ch, genotyped_indelrecal_ch, genotyped_applyrecal_ch

    script:
    """
    mkdir tmp
    gatk GenotypeGVCFs \
        -V gendb://$genomicsdb \
        -R $reference_fa \
        -O "genotyped.vcf" \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    """
}

/*
The next few processes do variant recalibration. SNPs and indels are recalibrated separately.
*/

// TODO:
// Increasing --max-gaussians may work for larger sample sizes. For two samples, --max-gaussians=4 failed. exoseq uses 4.
// Filtering with VQSR based on DP is not recommended for exome data. Don't know if I'm currently doing this.

// TODO: do I want the plots from "snps.plots.R" (and "snps.plots.R.pdf")?
// Generate recalibration and tranches tables for recalibrating the SNP variants in the next step.
process recalibrate_snps {
    input:
    //file vcf from snps_recalibrate_ch
    set file(vcf), file(idx) from genotyped_snprecal_ch

    output:
    set file("recal.table"), file("recal.table.idx") into snps_recal_table_ch
    file "tranches.table" into snps_trances_table_ch

    script:
    """
    mkdir tmp
    gatk VariantRecalibrator \
        -R $reference_fa \
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
        --rscript-file "snps.plots.R" \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    """
}

// TODO: use --trust-all-polymorphic?
// Generate recalibration and tranches tables for recalibrating the indel variants in the next step.
process recalibrate_indels {
    input:
    //file vcf from indels_recalibrate_ch
    set file(vcf), file(idx) from genotyped_indelrecal_ch

    output:
    set file("recal.table"), file("recal.table.idx") into indels_recal_table_ch
    file "tranches.table" into indels_trances_table_ch

    script:
    """
    mkdir tmp
    gatk VariantRecalibrator \
        -R $reference_fa \
        -V $vcf \
        -resource:mills,known=false,training=true,truth=true,prior=12.0 $mills \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode INDEL \
        --max-gaussians 4 \
        -O "recal.table" \
        --tranches-file "tranches.table" \
        --rscript-file "plots.plots.R" \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    """
}

// Realibrate indels.
process apply_vqsr_indels {
    input:
    //file vcf from indels_apply_ch
    set file(vcf), file(idx) from genotyped_applyrecal_ch
    set file(recal_table), file(recal_table_idx) from indels_recal_table_ch
    file tranches_table from indels_trances_table_ch

    output:
    set file("indels_recal.vcf"), file("indels_recal.vcf.idx") into recalibrated_indels_ch

    script:
    """
    mkdir tmp
    gatk ApplyVQSR \
        -R $reference_fa \
        -V $vcf \
        -O "indels_recal.vcf" \
        --truth-sensitivity-filter-level 99.0 \
        --tranches-file $tranches_table \
        --recal-file $recal_table \
        -mode INDEL \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    """
}

// Recalibrate SNPs.
process apply_vqsr_snps {
    input:
    //file vcf from snps_apply_ch
    set file(vcf), file(idx) from recalibrated_indels_ch
    set file(recal_table), file(recal_table_idx) from snps_recal_table_ch
    file tranches_table from snps_trances_table_ch

    output:
    set file("recalibrated.vcf"), file("recalibrated.vcf.idx")into recalibrated_vcf_ch

    script:
    """
    mkdir tmp
    gatk ApplyVQSR \
        -R $reference_fa \
        -V $vcf \
        -O "recalibrated.vcf" \
        --truth-sensitivity-filter-level 99.0 \
        --tranches-file $tranches_table \
        --recal-file $recal_table \
        -mode SNP \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    """
}

// TODO:
// Consider whether to use a supporting dataset. I commented out the "-supporting" argument,
// because it biases the data toward the population the supporting dataset is based on. If
// only few samples are in the dataset, a supporting dataset can be quite useful.
// Supply pedigree information (does this use more than just trios?)
// Consider which filters, if any, to apply.
// Calculate the genotype posteriors based on all the samples in the VCF.
process refine_genotypes {
    input:
    set file(vcf), file(idx) from recalibrated_vcf_ch

    output:
    set file("refined.vcf"), file("refined.vcf.idx") into refined_vcf_ch

    script:
    """
    gatk CalculateGenotypePosteriors \
        -V $vcf \
        -O "refined.vcf"
    #    -supporting $kGphase3
    """
}

// Add rsid from dbSNP
// NOTE: VariantAnnotator is still in beta (as of 20th of March 2019).
process annotate_rsid {
    input:
    set file(vcf), file(idx) from refined_vcf_ch

    output:
    set file("rsid_ann.vcf"), file("rsid_ann.vcf.idx") into rsid_annotated_vcf_ch

    script:
    """
    gatk VariantAnnotator \
        -R $reference_fa \
        -V $vcf \
        --dbsnp $dbsnp \
        -O "rsid_ann.vcf"
    """
}

// TODO: memory specification?
// Annotate the VCF with effect prediction. Output some summary stats from the effect prediction as well.
process annotate_effect {
    publishDir "${params.outdir}/variants", pattern: "snpEff_stats.csv", mode: 'copy', overwrite: true

    input:
    set file(vcf), file(idx) from rsid_annotated_vcf_ch

    output:
    file "effect_annotated.vcf" into effect_vcf_annotate_ch
    file "snpEff_stats.csv" into snpeff_stats_ch

    script:
    """
    snpEff -Xmx${params.mem}g \
         -i vcf \
         -o vcf \
         -csvStats "snpEff_stats.csv" \
         hg38 \
         -v \
         $vcf > "effect_annotated.vcf"
    """
}

process zip_and_index_vcf {
    publishDir "${params.outdir}/variants", mode: 'copy', overwrite: true

    input:
    file vcf from effect_vcf_annotate_ch

    output:
    set file("variants.vcf.gz"), file("variants.vcf.gz.tbi") into variants_evaluate_ch, variants_validate_ch

    script:
    """
    cat $vcf | bgzip -c > "variants.vcf.gz"
    tabix "variants.vcf.gz"
    """
}

// TODO: is this step necessary? Perhaps put it at the very end, to make sure the VCF is
// correctly formatted and such.
// Validate the VCF sice we used a non-GAKT tool.
process validate_vcf {
    publishDir "${params.outdir}/variants", mode: 'copy', overwrite: true, saveAs: { filename -> "validation.log" }
    input:
    set file(vcf), file(idx) from variants_validate_ch

    output:
    file ".command.log" into validation_log_ch

    script:
    """
    gatk ValidateVariants \
        -V $vcf \
        -R $reference_fa \
        --dbsnp $dbsnp
    """
}

process variant_evaluation {
    publishDir "${params.outdir}/variants", mode: 'copy', overwrite: true

    input:
    //file vcf from rsid_annotated_vcf_ch
    set file(vcf), file(idx) from variants_evaluate_ch

    output:
    file "variant_eval.table" into variant_eval_table_ch
    val "done" into status_ch

    script:
    """
    # VariantEval fails if the output file doesn't already exist. NOTE: this should be fixed in a newer version of GATK, as of the 19th of February 2019.
    echo -n > "variant_eval.table"

    gatk VariantEval \
        -R $reference_fa \
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


process multiqc {
    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

    input:
    val status from status_ch

    output:
    file "multiqc_report.html" into multiqc_report_ch
    file "multiqc_data" into multiqc_data_ch

    script:
    """
    multiqc -f ${params.outdir} --config ${params.multiqc_config}
    """
}

