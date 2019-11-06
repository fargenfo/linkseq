
import re, sys

# Read standard input.
rg = sys.stdin.read()

# Extract the "SM" field of the read group line.
# Match any sample name that doesn't contain space ("\s").
# Assign the sample name itself to a group via parentheses.
matches = re.search('SM:([^\s]*)', rg)

# Group 1 of the match is the sample name. Print to stdout.
print(matches.group(1))

