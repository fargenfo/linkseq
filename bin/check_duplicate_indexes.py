#!/usr/bin/python
#
# This scripts reads indexes in [Data] section of an Illumina sample sheet, and checks that they are all unique.
# If there are duplicate indexes, this script exits fatally, causing an error.

import sys

samplesheet = sys.argv[1]

fid = open(samplesheet)

while 'Lane,' not in fid.readline():
    continue

lines = fid.readlines()

indexes = []
for line in lines:
    line = line.strip()
    idx = line.split(',')[2]
    indexes.append(idx)

unique_indexes = list(set(indexes))

assert len(indexes) == len(unique_indexes), 'Error: there are duplicate indexes in samplesheet: %s' % samplesheet
