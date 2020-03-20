#!/usr/bin/env python3

"""
Author: Elisabet Thomsen
Date: 18-03-2020
Description: Program that trims barcodes off R2.
It checks if the 16 first bases in R1 are barcodes in the whitelist.
If positive, it reverse compliments the barcode, and checks if the rev_barcode is in the corresponding R2. If negative, it searches for rev_bc with one mismatch.
It searches from the 5' end and cuts the first match it finds.
If cuts line 4 to same length as line 2.
If R1 or R2 are shorter than 16 bp it will not trim.
Checks if line 2 contains nothing (''). If true, replaces '' with 'N' and replaces line 4 with '!'.
Throughput: Approx. 6381 read-pairs/second
Output: A barcode trimmed R2 fastq file (not gzipped).
Input: R1 and R2 must be in two seperate fastq files. If they are gzipped, they must end with '.gz'.
Usage: ./trimR2bc.py <R1fastqfile> <R2fastqfile> <whitelist> <R2outfile> 1>>bctrim_stats.txt
"""


import sys, gzip

# Check commandline
if len(sys.argv) < 5:
    print('Usage: ./trimR2bc.py <R1fastqfile> <R2fastqfile> <whitelist> <R2outfile>');
    sys.exit(1)

# Get filename from commandline
readfile = sys.argv[1]
readfile2 = sys.argv[2]
whitelistfile = sys.argv[3]
outfilename = sys.argv[4]

# Open files
try:
    # Check if file1 is gziped or not
    if readfile.lower().endswith('.gz'):
        infile = gzip.open(readfile,"rt")
    else:
        infile = open(readfile,"rt")
    # Check if file2 is gziped or not
    if readfile2.lower().endswith('.gz'):
        infile2 = gzip.open(readfile2, "rt")
    else:
        infile2 = open(readfile2, "rt")

    barcodefile = open(whitelistfile,"r")
    outfile = open(outfilename,'wt')
except IOError as err:
    print("Cant open file:", str(err));
    sys.exit(1)


# Define function
# Assumption: ksq and rsq have already been checked for 0 mismatch
def bi_search(ksq,rsq):
    flag = True
    l = len(ksq)
    # While there is not more than one mismatch
    while flag:
        l = len(ksq)
        # If this is not the last round
        if l != 1:
            # Calculate half of seq length
            half_l = int(len(ksq)/2)
            # Split the remaining seq into two
            k1 = ksq[0:half_l]
            k2 = ksq[half_l:l]
            r1 = rsq[0:half_l]
            r2 = rsq[half_l:l]
        # If this is the last round
        else:
            # If this last letter matches, then there is one mismatch
            if k1 == r1 or k2 == r2:
                onemismatch = True
                flag = False
            # If it doesnt match then stop this function
            else:
                flag = False
                onemismatch = False
        # If the first half matches, then prepare the other half for split
        if k1 == r1:
            ksq = k2
            rsq = r2
        # If the second half matches, then prepare the other half for split
        elif k2 == r2:
            ksq = k1
            rsq = r1
        # If none of the halfs matches, then there is more than one mismatch -> stop function
        else:
            flag = False
            onemismatch = False
    return onemismatch


# Make translation table for reverse complementation of barcode
complementtable = str.maketrans('ATCG', 'TAGC')

print("Infiles:", readfile, readfile2)
print("Whitelistfile:", whitelistfile)

# Make barcodelist (set)
barcodelist = set()

print("Making the barcodelist ...")

for line in barcodefile:
    barcodelist.add(line[:-1])

print("Done!")


# Initializing
count = 0
kmerlength = 16
chunksize = 200000
printthis = ''
readcount = 0
chunkreadcount = 0
trimfbccount = 0
trimhbccount = 0
totalbptrim = 0
totalbp = 0
Nflag = False

# Iterate through lines in both files
for line1, line2 in zip(infile, infile2):
    count += 1

    # Line two in fastq files contains the sequence
    if count == 2:
        readcount += 1
        # 16 first bases in R1 is the barcode
        bc = line1[0:16]
        R2sq = line2[:-1]
        R2len_pretrim = len(R2sq)

        # Only continue, if R1 and R2 are 16 bp or more and bc is in barcodelist
        if bc in barcodelist and len(bc) >= 16 and len(R2sq) >= 16:
            # Reverse complement the barcode
            rev_bc = bc.translate(complementtable)[::-1]

            ### Find rev_barcode in R2 with 0 or 1 mismatch ###
            bait1 = rev_bc[0:8]
            bait2 = rev_bc[8:16]

            pos1_fin = len(R2sq)
            pos2_fin = len(R2sq)
            fflag = False
            hflag = False

            # Look if first half of the bc can be found in R2seq
            pos1 = R2sq.find(bait1)
            pos2 = -1

            # While you still find the first half of bc in R2seq
            while pos1 != -1:
                ksq = bait2
                # Take the 8 bases after where bc matched
                rsq = R2sq[pos1+8:pos1+16]
                # Check if the second half has 0 mismatch
                if rsq == ksq:
                    R2sq = R2sq[:pos1]
                    trimfbccount += 1
                    fflag = True
                    break
                # If you find the other half with one mismatch
                elif bi_search(ksq,rsq) == True:
                    pos1_fin = pos1
                    hflag = True
                    break
                pos1 = R2sq.find(bait1, (pos1+1))

            # If first half AND second half not had 0 mismatch, then look for second half in R2seq
            if not fflag:
                # There is no need to look after a potential 1 mismatch from first half (stop looking at pos1_fin)
                pos2 = R2sq.find(bait2, 0, pos1_fin)

                # While you still find the second half bc in R2seq
                while pos2 != -1:
                    ksq = bait1
                    # Take the 8 bases infront of where bc matched
                    rsq = R2sq[pos2-8:pos2]
                    # Check for 1 mismatch
                    if bi_search(ksq,rsq) == True:
                        pos2_fin = (pos2-8)
                        hflag = True
                        break
                    pos2 = R2sq.find(bait2, pos2+1, pos1_fin)

                # Find the first 1 mismatch from 5' end
                if hflag:
                    if pos1_fin > pos2_fin:
                        pos_fin = pos2_fin
                    else:
                        pos_fin = pos1_fin
                    R2sq = R2sq[:pos_fin]
                    trimhbccount += 1

        # For stats (% bp trimmed off)
        R2len_posttrim = len(R2sq)
        bptrimlen = R2len_pretrim - R2len_posttrim
        totalbptrim += bptrimlen
        totalbp += R2len_pretrim

        # If R2sq = '', then replace with 'N'
        if len(R2sq) == 0:
            R2sq = 'N'
            Nflag = True

        printthis += R2sq + '\n'

    # Set the count to 0 again (fq files have 4 lines per read)
    elif count == 4:
        count = 0

        # If R2sq was '', then replace with '!'
        if Nflag == True:
            line2 = '!'
            Nflag = False

        # Trim the quality line to same length as R2
        printthis += line2[:len(R2sq)] + '\n'
        chunkreadcount += 1
        # Print in chunks of "chunksize" reads, to make algorithm faster
        if chunkreadcount == chunksize:
            outfile.write(printthis)
            printthis = ''
            chunkreadcount = 0

    # Line 1 and 3 need to be printed also
    else:
        printthis += line2


# Print the last chunk of reads
outfile.write(printthis)
print('BCtrimmedfile:', outfilename)


# Print stats
print('R2_Stats:')
print('Total_reads:', readcount)
perc = (trimfbccount/readcount) * 100
print('BCtrimmed(0_mismatch):', trimfbccount, '(' + '{:.2f}'.format(perc) + '%)')
perc = (trimhbccount/readcount) * 100
print('BCtrimmed(1_mismatch):', trimhbccount, '(' + '{:.2f}'.format(perc) + '%)')
print('Total_bases:', totalbp)
perc = (totalbptrim/totalbp) * 100
print('Total_bases_trimmed:', totalbptrim, '(' + '{:.2f}'.format(perc) + '%)')
print('\n')


# Close files
infile.close()
infile2.close()
barcodefile.close()
outfile.close()
