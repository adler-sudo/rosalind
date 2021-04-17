# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 22:28:30 2021

@author: james
"""

# take a dna strand and list of introns
# remove the introns form the dna strand
# perform translation and return the
# aa sequence

# import the modules
# the protein function is a custom function used to convert
# nucleotide sequence to aa sequence
from dna import protein
import argparse
import re


# parse the arguments
parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

# initate list of introns and dna sequence
introns = []
seq = ''

# separate full strand and introns
with open(args.filename) as f:
    f.readline()
    while True:
        newline = f.readline().rsplit()
        if newline[0][0] == '>':
            break
        seq = seq + newline[0]
    for n, line in enumerate(f.read().splitlines()):
        if n % 2 == 0:
            introns.append(line)

# remove introns from the sequence
for intron in introns:
    loc = re.search(intron, seq)
    start, end = loc.start(), loc.end()
    seq = seq[:start] + seq[end:]

# convert T to U
seq = seq.replace("T", "U")

# convert string to aa sequence
aa_seq = protein(seq)

# print aa seq to console
print(aa_seq)