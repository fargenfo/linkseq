#!/usr/bin/python
#
# This scripts reads sample names in [Data] section of an Illumina sample sheet, and checks that there are no
# underscores. If there are underscores in any sample names, this script exits fatally, causing an error.

import sys

samplesheet = sys.argv[1]

fid = open(samplesheet)

while 'Lane,' not in fid.readline():
    continue

lines = fid.readlines()

for line in lines:
    line = line.strip()
    samplename = line.split(',')[1]
    assert '_' not in samplename, 'Error: there must not be underscores in sample names in sample sheet: %s' % samplesheet

