
# GATK best-practices pipeline adapted to linked-reads


**TODO:**

* When done testing fix:
    * single_sample.nf
        * remove make_small_bam
        * remove 'targets = "chr17"'
        * fix qualimap_analysis
    * multi_sample.nf
        * remove 'targets = "chr17"'
* Make sure all the sample names match. For example, if the sample names in the genotyped VCF don't match the sample names of the BAM (either folder name or as defined in the read group), I may have problems joining channels by sample name.
    * If the sample name in the input FASTQ path CSV doesn't match the read group in the aligned BAM, then rename it in the BAM.
    * Make this an option, i.e. `param.rename_samples = false` by default.
    * Issue a warning as well, that the name doesn't match.
    * If you get this warning, you check that it's not an error but that you do want your samples to have different names.
    * The workflow prints which samples it has to rename.
    * The workflow creates a list of the samples it renamed, and includes it in the output directory.
* Help strings
* Each process should be allocated only as much memory as it needs. Profile the workflow to see the memory consumption
* Dockerize
* Phase VCF with HapCut2
    * Pretty much done. Check results.
* Use CRAM. Most likely just use BAM all the way through and convert to CRAM at the end
* Align reads with Lariat
* Parameter checking:
    * Check that the FASTQ path exists
    * Check that the FASTQ path has data
    * Check that "sample" is in the FASTQ path
    * Check that files such as dbsnp exist
    * Check that "gvcf_path" has at least one GVCF. Report number of GVCFs
* Unit tests and Continuous integration
    * Simulate small test data with LRSIM
        * https://github.com/aquaskyline/LRSIM
    * Make a tiny reference fasta. Maybe this can help:
        * https://github.com/SciLifeLab/Sarek-data/blob/9087faa53d25fca90c1a84a48cfaf7cbed496317/scripts/makeShortContigs.py
        * Also have to package the reference with longranger mkref
    * Make tiny dbsnp, target BED, and so on

