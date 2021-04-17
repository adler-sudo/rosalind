# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 17:30:52 2021

@author: james
"""

# impor tmodules
import argparse


# establish cmd parser
parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

# get sequence from file
with open(args.filename) as f:
    f.readline()
    seq = f.read().splitlines()
    
seq = ''.join(seq)
print(seq)

# establish reverse compliment tools
old_chars = "ACGT"
replace_chars = "TGCA"
tab = str.maketrans(old_chars, replace_chars)

# identify site indeces and length and print to console
for a in range(len(seq)):
    i = a + 3
    while i - a <= 12 and i < len(seq) + 1:
        fawd = seq[a:i]
        bawd = seq[a:i].translate(tab)[::-1]
        if fawd == bawd:
            print(a+1, len(fawd))
        i += 1
        