# -*- coding: utf-8 -*-
"""
Takes a piece of DNA from FASTA file thought to contain proteins,
and returns a list of amino acid strings that are likely protein candidates 
that can be identified on BLAST. The screaning process looks for start and 
stop codons, as well as a minimum length, determined by scrambling the input
DNA sequence and looking for the longest non-coding protein (anything longer)
than the non-coding protein is likely not by chance).

This code is currently set up to identify proteins from a piece of salmonella
DNA. 

BLAST search: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome

Notes on the code: The first 6 functions were tested using a python visualization
site since I didn't understand doctest implementation at the time. The final
functions were tested via doctests (in each function, modified final call to 
current function) and running the code with print statements in the terminal
and sublime build box (Ctrl+B).

@author: Elizabeth Sundsmo

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq #this is a function, give a FASTA file to read DNA seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return 'T'  
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G' 
    elif nucleotide == 'G':
        return 'C'
    else:
        return -1

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    i = -1
    reverse_complement = ''

    while i >= -len(dna):
        reverse_complement+=(get_complement(dna[i]))
        i += -1
    return reverse_complement

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    i = 0
    i2 = 3
        
    while i < len(dna):
    	if dna[i:i2] in ['TAA', 'TAG', 'TGA']:
            return dna[:i]
        i+=3
        i2+=3
            
        if i2> len(dna):
           return dna

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    
    all_ORFs = []
    count = 0
    i = 0
    i2 = 3

    while i < len(dna):
        if dna[i:i2] == 'ATG':
            dna_piece = dna[i:]
            all_ORFs.append(rest_of_ORF(dna_piece))
            
            i+=3+len(all_ORFs[count])
            i2= i+3
            count +=1
        else :    
            i +=3
            i2+=3

    return all_ORFs

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_ORFs = []

    all_ORFs.extend(find_all_ORFs_oneframe(dna[0:]))
    all_ORFs.extend(find_all_ORFs_oneframe(dna[1:]))
    all_ORFs.extend(find_all_ORFs_oneframe(dna[2:]))
    
    return all_ORFs

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    all_ORFs_both_strands = []
    all_ORFs_both_strands.extend(find_all_ORFs(dna))
    all_ORFs_both_strands.extend(find_all_ORFs(get_reverse_complement(dna)))

    return all_ORFs_both_strands

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """  
    all_strands = find_all_ORFs_both_strands(dna)
    big = ''
    i = 1
    big = all_strands[0]

    while i < len(all_strands):
        if len(all_strands[i]) > len(all_strands[i-1]):
            big = all_strands[i]
        i+=1
    
    return big 

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 
    """
    #print type(longest_ORF(shuffle_string(dna)))
    returns = ''
    t = 0
    while t < num_trials:
    	shuffled = longest_ORF(shuffle_string(dna))
        if len(shuffled) > len(returns):
    	    returns = shuffled
        t+=1
    return returns

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
        
        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA("ATGCGAATG")
        'MRM'
        >>> coding_strand_to_AA("ATGCTACATTCGCAT")
        'MLHSH'
    """
    returns = ''
    i = 0
    i2 = 3
        
    while i < len(dna):
    	returns += aa_table[dna[i:i2]]
        i+=3
        i2+=3         
    return returns

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
        
        >>> gene_finder("ATGCGAATGTAGCATCAAA")
        ['MRM', 'MLHSH']
    """
    print 'Finding threshold...'
    threshold = len(longest_ORF_noncoding(dna, 1500))
    print threshold
    
    all_ORFs = find_all_ORFs_both_strands(dna) 
    returns = []
    i=0

    print 'Entering while loop... '
    while i < len(all_ORFs):
        if len(all_ORFs[i]) >= threshold:
    	    returns.append(coding_strand_to_AA(all_ORFs[i]))
    	    print 'Added sequence'
        i+=1

    print str(len(returns)) + 'sequences added: \n'
    print returns


if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    #doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)

gene_finder(load_seq("./data/X73525.fa"))