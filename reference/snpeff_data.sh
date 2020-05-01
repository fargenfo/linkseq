#!/bin/bash
#
# Download the reference data for SnpEff. The file downloaded below corresponds to hg38, and the link
# can be found by running the following command:
# snpEff databases | grep hg38
# This data can be used using the -dataDir parameter in snpEff.

wget http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_hg38.zip
unzip snpEff_v4_3_hg38.zip
rm snpEff_v4_3_hg38.zip
mv data snpeff_data
