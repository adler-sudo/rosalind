# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 14:51:11 2020

@author: james
"""

# functions for solving dna-related problems

def rev(strand):
    """returns the reverse of the string"""
    a = strand[::-1]
    return a


def compliment(strand):
    """returns the compliment of a dna strand"""
    
    # create the dictionary for variable assignment
    a = {'G':'C', 'C':'G', 'A':'T', 'T':'A'}
    
    # create the empty list
    h = []
    
    # create and return the compliment strand
    for i in strand:
        t = a[i]
        h.append(t)
        
    h = "".join(h)
    
    return h


def gc_content(strand):
    l = len(strand)
    gc = strand.count('G') + strand.count('C')
    gc_perc = gc / l * 100
    print(gc_perc)


def strand_compare(strand1, strand2):
    sum = 0
    for i in range(len(strand1)):
        if strand1[i] != strand2[i]:
            sum += 1
    return sum

    
codons = {'UUU': 'F',
 'CUU': 'L',
 'AUU': 'I',
 'GUU': 'V',
 'UUC': 'F',
 'CUC': 'L',
 'AUC': 'I',
 'GUC': 'V',
 'UUA': 'L',
 'CUA': 'L',
 'AUA': 'I',
 'GUA': 'V',
 'UUG': 'L',
 'CUG': 'L',
 'AUG': 'M',
 'GUG': 'V',
 'UCU': 'S',
 'CCU': 'P',
 'ACU': 'T',
 'GCU': 'A',
 'UCC': 'S',
 'CCC': 'P',
 'ACC': 'T',
 'GCC': 'A',
 'UCA': 'S',
 'CCA': 'P',
 'ACA': 'T',
 'GCA': 'A',
 'UCG': 'S',
 'CCG': 'P',
 'ACG': 'T',
 'GCG': 'A',
 'UAU': 'Y',
 'CAU': 'H',
 'AAU': 'N',
 'GAU': 'D',
 'UAC': 'Y',
 'CAC': 'H',
 'AAC': 'N',
 'GAC': 'D',
 'UAA': 'Stop',
 'CAA': 'Q',
 'AAA': 'K',
 'GAA': 'E',
 'UAG': 'Stop',
 'CAG': 'Q',
 'AAG': 'K',
 'GAG': 'E',
 'UGU': 'C',
 'CGU': 'R',
 'AGU': 'S',
 'GGU': 'G',
 'UGC': 'C',
 'CGC': 'R',
 'AGC': 'S',
 'GGC': 'G',
 'UGA': 'Stop',
 'CGA': 'R',
 'AGA': 'R',
 'GGA': 'G',
 'UGG': 'W',
 'CGG': 'R',
 'AGG': 'R',
 'GGG': 'G'}


rev_codons = {'F': ['UUU', 'UUC'],
 'L': ['CUU', 'CUC', 'UUA', 'CUA', 'UUG', 'CUG'],
 'I': ['AUU', 'AUC', 'AUA'],
 'V': ['GUU', 'GUC', 'GUA', 'GUG'],
 'M': ['AUG'],
 'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
 'P': ['CCU', 'CCC', 'CCA', 'CCG'],
 'T': ['ACU', 'ACC', 'ACA', 'ACG'],
 'A': ['GCU', 'GCC', 'GCA', 'GCG'],
 'Y': ['UAU', 'UAC'],
 'H': ['CAU', 'CAC'],
 'N': ['AAU', 'AAC'],
 'D': ['GAU', 'GAC'],
 'Stop': ['UAA', 'UAG', 'UGA'],
 'Q': ['CAA', 'CAG'],
 'K': ['AAA', 'AAG'],
 'E': ['GAA', 'GAG'],
 'C': ['UGU', 'UGC'],
 'R': ['CGU', 'CGC', 'CGA', 'AGA', 'CGG', 'AGG'],
 'G': ['GGU', 'GGC', 'GGA', 'GGG'],
 'W': ['UGG']}



def chunks(s, n):
    """converts string to specified number of 
    chunks. Yields a generator object.
    
    s: string to be split up
    
    n: size of each chunk
    """
    
    for start in range(0, len(s), n):
        yield s[start:start+n]


def protein(strand):
    """provides amino acid sequence of protein
    based on rna strand provided
    
    strand: rna strand to be transmitted
    """
    from dna import codons
    from dna import chunks
    
    aa_seq = []
    
    for chunk in chunks(strand, 3):
        aa = codons[chunk]
        if aa == 'Stop':
            break
        aa_seq.append(aa)
        
    aa_seq = ''.join(aa_seq)

    return aa_seq


def substringOccurrence(strand1, substring):
    """determine the number of times that a substring
    appears in a strand
    
    strand1: the strand to be analyzed
    
    substring: the sequence to look for
    
    """
    
    sum = 0
    
    for i in range(len(strand1)):
        u = strand1[i:(i+4)]
        if u == substring:
            sum += 1

    return sum


def substringLocation(strand1, substring):
    """determine the position of a substring within 
    a strand 
    
    strand1: the strand to be analyzed
    
    substring: the sequence to look for 
    """
    
    positions = []
    
    for i in range(len(strand1)):
        u = strand1[i:(i+len(substring))]
        if u == substring:
            pos = i+1
            positions.append(pos)
            
    return positions


def rosalindSplit(file):
    """returns list of raw rosalind strands from text file
    
    file: text file containing strand
    
    """
    
    # open file
    with open(file, 'r') as f:
        raw = f.read()
        
    # drop and split up each strand
    raw = raw.replace("\n", "")
    q = raw.split(">Rosalind_")
    q.pop(0)
    
    # drop all remaining numbers and return strands in list
    nums = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    
    for i in range(len(q)):
        for letter in q[i]:
            if letter in nums:
                q[i] = q[i].replace(letter, "")
    
    return q


def rosalindSplitNames(file):
    """returns list of rosalind names associated
    with strands returned from rosalindSplit()
    
    file: text file containing strand
    
    """
    
    # open file
    with open(file, 'r') as f:
        raw = f.read()
        
    # drop and split up each strand
    raw = raw.replace("\n", "")
    raw = raw.split(">")
    raw.pop(0)
    
    names = []
    
    for i in raw:
        names.append(i[:13])
        
    return names


def matrixStrand(strand_list, file):
    """returns matrix values of each nucleotide base in a list
    of strands and returns the most common nucleotide for 
    each position
    
    strand_list: list of strands to compare (formatted properly
    when coming from function 'rosalindSplit')
    
    """
    
    # estblish empty lists
    A = []
    C = []
    G = []
    T = []
    dom_list = []
    
    # calculate occurrences of each letter in each strand
    for i in range(len(strand_list[0])):
            A.append(0)
            C.append(0)
            G.append(0)
            T.append(0)
    
    for q in range(len(strand_list)):
        for s in range(len(strand_list[0])):
            if strand_list[q][s] == "A":
                A[s] += 1
            elif strand_list[q][s] == "C":
                C[s] += 1
            elif strand_list[q][s] == "G":
                G[s] += 1
            elif strand_list[q][s] == "T":
                T[s] += 1

    for i in range(len(A)):
        a = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        a['A'] = A[i]
        a['C'] = C[i]
        a['G'] = G[i]
        a['T'] = T[i]
        
        dom_list.append(max(a, key = a.get))
    
    # creaate dominant string
    dom_list = "".join(dom_list)
    
    # return in desired matrix format
    with open(file, 'w') as f:
        f.write(dom_list)

    # convert nucleotide lists to strings and print in matrix format
    A = ['A:'] + A
    C = ['C:'] + C
    G = ['G:'] + G
    T = ['T:'] + T
    
    def stringify(li):
        """turn a list of values into a list of string values
        
        li: list to alter
        
        """
        
        for o in range(len(li)):
            li[o] = str(li[o])
        li = " ".join(li)
        return li
    
    A = stringify(A)
    C = stringify(C)
    G = stringify(G)
    T = stringify(T)
    
    with open(file, 'w') as f:
        f.write(dom_list + "\n" + A + "\n" + C + "\n" + G + "\n" + T)


def adjacentList(file, k):
    """returns an adjacent list given a file and
    match of length k
    
    file: file to analyze (contains FASTA strands)
    
    k: length of match
    
    """
    
    # import functions for dealing with rosalind text files
    from dna import rosalindSplit, rosalindSplitNames
    
    # split the strands and names in the file
    strands = rosalindSplit(file)
    names = rosalindSplitNames(file)
    
    for i in range(len(strands)):
        for t in range(len(strands)):
            if strands[i] == strands[t]:
                pass
            else:
                if strands[i][-k:] == strands[t][:k]:
                    print(names[i], names[t])


def longestCommonSubstring(file):
    """returns longest common substring of
    a list of strands in a file
    
    file: FASTA format file with a list of
    strands
    
    """
    
    # split up the dna strands in the file
    from dna import rosalindSplit
    
    x = rosalindSplit(file)
    
    # prepare empty dictionary
    a = {}

    # create a dictionary using the first strand. 
    # set values to 0 so that when the first strand is looped through
    # it will only set values to 1
    for i in range(len(x[0])):
        for b in range(len(x[0])):
            q = x[0][i:b]
            if q == '':
                pass
            else:
                a[q] = 0

    # check if substrings are present within dictionary and set value
    # to i, eliminating strings that do not have a value of i in order
    # to continuously reduce the size of the dictionary and speed up
    # run time
    for i in range(len(x)):
        for r in range(len(x[i])):
            for b in range(len(x[i])):
                if b > r:
                    q = x[i][r:b]
                    if q in a.keys():
                        a[q] = i
        a = {k:v for k, v in a.items() if a[k] == i}
    
    # create new dictionary with sequences that survived and their lengths
    # returning the substring with the longest length
    c = {k:len(k) for k in a.keys()}
    t = max(c, key = c.get)
    return t

            
def uniprotReader(protein, motif):
    """DRAFT this function will take a file with
    uniprot sequences and locate the positions of
    a given motif
    
    file: file of uniprot proteins
    
    motif: motif to search for in the sequence
    of amino acids (the motif will be written in 
    regular expression)
    """
    
    # import necessary modules
    import re
    from urllib.request import urlopen
    

    # read in and convert bytes to string
    s = urlopen("https://www.uniprot.org/uniprot/{}.fasta".format(protein))
    t = s.read()
    m = t.decode("utf-8")
    
    # isolate aa sequence
    m = m.replace("\n", "")
    a = m.find("SV=")
    n = m[a:]
    n = re.sub("SV=\d+", "", n)


    
    # compile list of all strands
    a = []
    
    for o in range(len(n)):
        mo = n[o:o+4]
        a.append(mo)

    b = []

    for k in range(len(a)):
        if re.match(motif, a[k]):
            b.append(k)
        # if it matchess add i to list

    print(protein)
    print(" ".join([str(i) for i in b]))


def uniprotReaderStolen():
    from re import finditer
    from sys import argv
    from urllib.request import urlopen
    
    # the % sign will be used later in a for loop
    # quicker way of running through the whole problem
    uniprot = "http://www.uniprot.org/uniprot/%s.fasta"
    
    # quickly and neatly pulls the lines from the fasta file
    # doesn't work with the strand sequence rosalind files i was working
    # with earlier in this file
    for protein in open(argv[1], 'r').read().strip().splitlines():
    
        # Fetch the protein's fasta file and get rid of newlines.'
        f = urlopen(uniprot % protein).read().decode('utf-8')
        f = ''.join(f.splitlines()[1:])
    
        # Find all the positions of the N-glycosylation motif.
        # performing the regex statement like shown avoids overlap misses
        locs = [g.start()+1 for g in finditer(r'(?=N[^P][ST][^P])', f)]
    
        # Print them out, if any.
        if locs != []:
            print(protein)
            print(' '.join(map(str, locs)))


def mrnaInfer(protein):
    """returns the number of possible mRNA
    sequences that could have produced the given
    protein (amino acid sequence)
    
    protein: the amino acid sequence to be
    analyzed
    
    """
    
    # import modules
    from functools import reduce
    from operator import mul
    from dna import rev_codons

    revCodonLengths = {k: len(v) for k, v in rev_codons.items()}

    # split the amino acids
    protein = list(protein)
    
    # map amino acid potential inputs
    mapper = list(map(revCodonLengths.get, protein))

    # caluculate possible combos
    nums = reduce(mul, mapper) * 3
    
    # modulo 1000000
    mod = nums % 1000000

    # return
    return mod


def potentialProteins(input_strand):
    """find potential proteins given one
    strand of DNA double helix
    
    strand: strand provided
    
    """
    # import necessary modules
    import re
    from dna import compliment
    from dna import protein
    
    # establish forward and complimentary strands
    forward = input_strand
    backward = compliment(forward)[::-1]

    strands = [forward, backward]

    # establish stop codons
    stops = ['TAA', 'TAG', 'TGA']

    # find viable sequences
    sequences = []
    
    # loop through each strand
    for strand in strands:
        
        # locate start codons
        begins = []
        
        for i in re.finditer("(?=ATG)", strand):
            begins.append(i.start())
        
        # locate stop codons        
        ends = []

        for stop in stops:
            seq = "(?=%s)"
            for i in re.finditer(seq % stop, strand):
                ends.append(i.start())

        for begin in begins:
            viableSequences = [end for end in ends if end > begin and end % 3 == begin % 3]
            if len(viableSequences) > 0:
                sequences.append(strand[begin:min(viableSequences)])
    
    # compile list of proteins
    proteins = []
    
    for sequence in sequences:
        sequence = sequence.replace('T', 'U')
        pro = protein(sequence)
        if pro in proteins:
            pass
        else:
            proteins.append(pro)
                
    return proteins

   



