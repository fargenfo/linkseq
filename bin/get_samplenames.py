#!/usr/bin/env python3


import sys, re

fastq_path = sys.argv[1]

# Get the filename from the path.
filename = re.split('/', fastq_path)[-1]

# Get the sample name from the filename.
sample = re.split('_', filename)[0]

print(sample, end='')

