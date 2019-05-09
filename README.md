
# GATK best-practices pipeline adapted to linked-reads



* Help strings
* Each process should be allocated only as much memory as it needs. Profile the workflow to see the memory consumption.
* Always output index files a la 'set file("1.vcf"), file("1.vcf.idx") into variants_ch'
* Dockerize
* Phase VCF with HapCut2
* Phase BAM with WhatsHap using phased VCF
* Use CRAM. Most likely just use BAM all the way through, convert to CRAM at the end, and delete all intermediate files.

