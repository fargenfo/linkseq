#!/usr/bin/env python3
#
# This script assumes the FASTQ sequence identifier is formatted as follows:
# @Instrument:RunID:FlowCellID:Lane:Tile:X:Y:UMIRead:Filter:0:IndexSequence
#
# And that the gzip compressed fastq_path is formatted as:
# [sample name]_[lane number]_[read number]_[].fastq.gz
#
# The readgroup will be formatted as:
# @RG\tID:[sample name]\tPL:Illumina\tPM:[platform model]\tPU:[flow cell ID]\tSM:[sample name]


import sys, gzip, re

fastq_path = sys.argv[1]
sample = sys.argv[2]

# Platform is always Illumina.
platform = 'Illumina'

# Open gzipped FASTQ file in "read text" mode.
with gzip.open(fastq_path) as fid:
    # Read the first line and remove the newline.
    line = fid.readline().strip()

# Get the instrument and the flow cell ID. Remove the "@" from the instrument ID.
instrument = re.split(b':', line)[0][1:]
flowcell = re.split(b':', line)[2]

instrument = instrument.decode()
flowcell = flowcell.decode()

# Read group format.
rg_format = "@RG\\tID:%s\\tPL:%s\\tPM:%s\\tPU:%s\\tSM:%s"

print(rg_format %(sample, platform, instrument, flowcell, sample), end='')

