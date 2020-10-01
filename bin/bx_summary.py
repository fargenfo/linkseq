#!/usr/bin/env python3
'''
Summarize output from `bxtools stats` (https://github.com/walaj/bxtools).
'''

import sys
import numpy as np

csv = sys.argv[1]

bx_counts = []
with open(csv) as fid:
    # Discard the header.
    temp = fid.readline()
    for line in fid:
        # Split the tab separated line into fields, and get the count (the second field).
        count = line.split('\t')[1]
        # Convert it to an integer and append to list.
        bx_counts.append(int(count))

    # Convert to a numpy array.
    bx_counts = np.array(bx_counts)

# Number of unique barcodes.
n_bx = len(bx_counts)

# Maximum, minimum, median and mean of barcode count.
n_max = max(bx_counts)
n_min = min(bx_counts)
median = np.median(bx_counts)
mean = np.mean(bx_counts)

# Number of barcodes with only one read.
singletons = sum(bx_counts == 1)

print_format = "Unique barcodes: %d\nPer-barcode read counts:\nMax: %d\nMin: %d\nMedian: %d\nMean: %f\nNumber of singletons: %d"

print(print_format %(n_bx, n_max, n_min, median, mean, singletons))

