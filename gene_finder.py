# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: ANA KRISHNAN

"""
from load import load_seq
dna = load_seq("./data/X73525.fa")


import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


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

        I  don't need more doctests because they are all basically the same.

    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # TODO: implement this
    if nucleotide == 'A':
        return('T')
    if nucleotide == 'T':
        return('A')
    if nucleotide == 'C':
        return('G')
    if nucleotide == 'G':
        return('C')
    pass


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
    # TODO: implement this
    length = len(dna)-1
    revcomp = ''
    while length >= 0:
        revcomp = revcomp + get_complement(dna[length])
        length = length - 1
    return revcomp
    pass


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
    # TODO: implement this
    x = 3
    while x < len(dna):
        if dna[x:x+3] == "TAG" or dna[x:x+3] == "TGA" or dna[x:x+3] == "TAA":
            return dna[:x]
        x = x + 3
    return dna
    pass


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
    # TODO: implement this
    x = 0
    ORFs = []
    while x < len(dna):
        if dna[x:x+3] == 'ATG':
            ORF = rest_of_ORF(dna[x:])
            ORFs.append(ORF)
            x = x + len(ORF)
        else:
            x = x + 3
    return ORFs
    pass


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
    # TODO: implement this
    frameOne = find_all_ORFs_oneframe(dna)
    frameTwo = find_all_ORFs_oneframe(dna[1:])
    frameThree = find_all_ORFs_oneframe(dna[2:])
    return frameOne + frameTwo + frameThree
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    ORFs = []
    ORFs = find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))
    return ORFs
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    ORFs = find_all_ORFs_both_strands(dna)
    longest = max(ORFs, key=len)
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    for x in range(0, num_trials):
        dna = shuffle_string(dna)
        x = x + 1
    longest = longest_ORF(dna)
    int = len(longest)
    return int


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
    """
    x = 0
    aminos = ''
    while x < len(dna):
        codon = dna[x:x+3]
        if len(codon) == 3:
            amino = aa_table[codon]
            aminos += amino
        x = x + 3
    return aminos
    # hacked


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    aminoSequence = []
    threshold = longest_ORF_noncoding(dna, 1500)
    ORFs = find_all_ORFs_both_strands(dna)
    for ORF in ORFs:
        if len(ORF) > threshold:
            aminoSequence.append(coding_strand_to_AA(ORF))
    return aminoSequence


gene_finder(dna)


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
