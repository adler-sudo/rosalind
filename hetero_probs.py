# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 11:38:24 2021

@author: james
"""

import argparse
import math


parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()
with open(args.filename) as f:
    content = f.readline().split()
    generation = int(content[0])
    heteros = int(content[1])
    
def probs_hetero(generation, heteros):
	kin = 2**generation
	probs = 0
	for i in range(heteros, kin+1):
		fact = math.factorial(kin) / (math.factorial(i) * math.factorial(kin-i))
		yes = 0.25 ** i
		no = 0.75 ** (kin - i)
		probs += fact * yes * no
	return probs

probs = probs_hetero(generation, heteros)

print(generation)
print(heteros)
print(probs)
