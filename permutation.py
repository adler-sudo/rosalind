# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 15:43:49 2021

@author: james
"""

# takes a number and returns the total possible permutations
# followed by each of said permutations

# import modules
import math
import argparse


# get the argument file
parser = argparse.ArgumentParser()
parser.add_argument('filename')
parser.add_argument('output_file')
args = parser.parse_args()

# open the file and read the range value
with open(args.filename) as f:
    content = f.readline()[0]
    content = int(content)
    
# initiate list
l = [*range(1, content+1)]
    
# write permutation function
def permutation(l):
    
    # if no values left in list, return
    if len(l) == 0:
        return []
    
    # if only one value left return that value
    if len(l) == 1:
        return [l]
    
    # initiate new empty list
    y = []
    
    # if more than a value left then loop through
    # we have kind of a nested function within itself here
    # almost like a circular reasoning?
    for i in range(len(l)):
        
        # pick your value
        m = l[i]
        
        # track remaining values
        reml = l[:i] + l[i+1:]
        
        # this is where we make that circular loop
        for p in permutation(reml):
            y.append([m] + p)
    
    return y
        
# caluclate total possible values
cnt_combos = math.factorial(content)

# pull output file from args
output_file = args.output_file

# write contents to file
with open(output_file, 'w') as w:
    w.write("%s" % cnt_combos)
    w.write("\n")
    
with open(output_file, 'a') as w:
    for p in permutation(l):
        for item in p:
            w.write("%s " % item)
        w.write("\n")


