#!/usr/bin/python
#
# Check that the sample sheet is proper. We assume the following about the sample sheet:
# * It is comma-separated.
# * The "data" section contains minimum the following fields: Lane, Sample_ID, and index.
# * The indexes are unique.
# * The sample IDs do not contain underscores.
# * The "data" section is the last section of the samplesheet.
#
# Usage:
# python check_duplicate_indexes.py SampleSheet.csv
# Where "SampleSheet.csv" is an Illumina sample sheet (see project README.md).

import sys

samplesheet = sys.argv[1]

fid = open(samplesheet)

# Check the "Data" section of the file.

# Discard all lines to the header of the [Data] section (inclusive).
while '[Data]' not in fid.readline():
    continue

# Skip empty lines until the header is reached.
for line in fid:
    # Check whether the line is empty or only whitespace.
    if line.strip() != '':
        table_header = line
        break

# Get the list of fields in the header.
table_header = table_header.split(',')

# Check if the header is comma-separated.
# Assuming that if it contains commas, then it is comma-separated.
# This seems like a fair assumption.
assert len(table_header) > 1, 'Error: sample sheet must be comma separated.'

# Check that the header has the correct fields.

# Remove potential whitespace from fields.
table_header = [field.strip() for field in table_header]

# FIXME: get index of fields of interest.
# FIXME: check that fields of interested are IN list, other fields may be there, and the order is not important.

# The table should have these fields.
header_fields = ['Lane', 'Sample_ID', 'index']

for field in header_fields:
    assert field in table_header, 'Error: field "%s" not found in sample sheet. Data section must contain the following fields: Lane, Sample_ID, and index.' % field

field_indexes = { field : table_header.index(field) for field in header_fields }

# Read the remaining lines, which is a table that should contain lanes, sample names, indexes, and so on.
lines = fid.readlines()

# Check that the indexes are all unique.
indexes = []
for line in lines:
    line = line.strip()
    idx = line.split(',')[field_indexes['index']]
    indexes.append(idx)


unique_indexes = list(set(indexes))

assert len(indexes) == len(unique_indexes), 'Error: there are duplicate indexes in samplesheet: %s' % samplesheet

# Check that no of the sample names have underscores in them.
for line in lines:
    line = line.strip()
    samplename = line.split(',')[field_indexes['Sample_ID']]
    assert '_' not in samplename, 'Error: there must not be underscores in sample names in sample sheet: %s' % samplesheet

