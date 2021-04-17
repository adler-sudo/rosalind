# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 15:15:46 2020

@author: james
"""

# functions for dealing with population dynamics

def mortalRabbits(n, m):
    """returns the total number of rabbits in n-th month if
    each rabbit pair produces a pair of offspring each month
    following maturity (1 month) and each rabbit lives for m months
    
    n: month of rabbit count
    
    m: number of months each rabbit lives
    
    """
    
    rabbits = [1, 1]
    
    i = 2
    
    # sum rabbits for last two months and subtract rabbits born m
    # months ago (which is equal to rabbits in population the month
    # before m months ago(i-m-1))
    while i < n:
        # before any rabbits are dying
        if i < m:
            population = rabbits[i-1] + rabbits[i-2]
        
        # the first month that rabbits die, you only lose the first pair
        # this avoids the possibility of taking the -1 of the list, which
        # would give the last value in the list
        elif i == m:
            population = rabbits[i-1] + rabbits[i-2] - rabbits[0]
        
        # for all subsequent months, subtract the rabbit population born
        # i-m months ago, which is equal to the population of mature
        # rabbits the generation before (i-m-1)
        elif i > m:
            population = rabbits[i-1] + rabbits[i-2] - rabbits[i-m-1]

        rabbits.append(population)
        i += 1
    
    return rabbits[n-1]