# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 16:32:01 2021

@author: james
"""

# a function that takes a protein sequence and returns
# the molecular weight
import argparse


# establish argument parser
parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

# monoisotopic mass table
weights = dict(A = 71.03711,
C = 103.00919,
D = 115.02694,
E =  129.04259,
F = 147.06841,
G = 57.02146,
H = 137.05891,
I = 113.08406,
K = 128.09496,
L = 113.08406,
M = 131.04049,
N = 114.04293,
P = 97.05276,
Q = 128.05858,
R = 156.10111,
S = 87.03203,
T = 101.04768,
V = 99.06841,
W = 186.07931,
Y = 163.06333)

# protein sequence of aa
with open(args.filename) as f:
    protein_seq = f.readline().rsplit()

# initiate weight variable
total_weight = 0

# calculate sum of weight sof aas
for letter in list(protein_seq[0]):
    total_weight += weights[letter]

# print the weight total to the terminal
print(total_weight)
