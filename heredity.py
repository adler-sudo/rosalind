# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 21:37:05 2020

@author: james
"""

# functions for solving problems related to heredity

def domProb(k, m, n):
    """determines probability of offspring displaying
    dominant phenotype given:
        
        k: number of homozygous dominant individuals
        
        m: number of heterozygous individuals
        
        n: number of homozygous recessive individuals
        
        """
    
    # total number of individuals
    t = k + m + n
    
    # calculation of mating probability and subsequent recessive phenotype
    rec_mate_self = (n-1) / (t-1) * (4/4)
    
    rec_mate_het = m / (t-1) * (2/4)
    
    het_mate_self = (m-1) / (t-1) * (1/4)
    
    het_mate_rec = n / (t-1) * (2/4)
    
    # probability of recessive phenotype
    prob_rec = (rec_mate_self + rec_mate_het) * (n/t) + (het_mate_self + het_mate_rec) * (m/t)
    
    # probbility of dominant phenotype
    prob_dom = 1 - prob_rec
    
    return prob_dom
    

def randomVariable(a, b, c, d, e, f):

    a1 = a * 4 * .5
    
    b1 = b * 4 * .5
    
    c1 = c * 4 * .5
    
    d1 = d * 3 * .5
    d2 = d * 1 * .5
    
    e1 = e * 2 * .5
    e2 = e * 2 * .5
    
    f2 = f * 4 * .5
    
    total = a1 + b1 + c1 +d1 + e1
    
    return total
