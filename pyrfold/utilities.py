"""collection of data and utility functions for general RNA and
DNA work
"""
from random import random, choice

# Dictionaries of complement sequences
RNA_complement_dict = {'A': 'U',
                       'U': ['A', 'G'],
                       'G': ['C', 'U'],
                       'C': 'G'}

DNA_complement_dict = {'A': 'T',
                       'T': 'A',
                       'G': 'C',
                       'C': 'G'}


def random_RNA_sequence(size, GC_range=None):
    """
    Simple Random RNA generating random sequence generator
    """
    for i in range(10000):
        out = []
        for placecounter in range(size):
            out.append(choice(['A', 'U', 'G', 'C']))
        if GC_range:
            if GC_range[0] <= GC_content(out) <= GC_range[1]:
                return ''.join(out)
        else:
            return ''.join(out)


def random_DNA_sequence(size, GC_range=None):
    """
    Simple Random DNA generating random sequence generator
    """
    for i in range(10000):
        out = []
        for placecounter in range(size):
            out.append(choice(['A', 'T', 'G', 'C']))
        if GC_range:
            if GC_range[0] <= GC_content(out) <= GC_range[1]:
                return ''.join(out)
        else:
            return ''.join(out)


def RNA_sequences_complementary(sequence1, sequence2):
    """
    Function to determine if two RNA sequences are complementary
    """
    # Make sure the sequences are RNA
    sequence1 = convert_to_RNA(sequence1)
    sequence2 = convert_to_RNA(sequence2)
    if len(sequence1) != len(sequence2):
        # write error script
        pass
    for position, base1 in enumerate(sequence1):
        if sequence2[-(position+1)] in RNA_complement_dict[base1]:
            continue
        else:
            return False
    return True


def randomly_substitue_u_for_c(sequence, probability=0.5):
    """Function for randomly_substituing Us for Cs for a given sequence
    """
    sequence = str(sequence)
    outsequence = []
    for base in sequence:
        if base == 'C':
            if random() <= probability:
                base = 'U'
        outsequence.append(base)
    return ''.join(outsequence)


def GC_content(sequence):
    sequence = str(sequence)
    GC_count = sequence.count('G')
    GC_count += sequence.count('C')
    return float(GC_count)/len(sequence)


def convert_to_RNA(sequence):
    sequence = sequence.upper()
    return sequence.replace('T', 'U')


def convert_to_DNA(sequence):
    """
    Converts RNA to DNA
    """
    sequence = str(sequence).upper()
    return sequence.replace('U', 'T')


def reverse_complement(sequence, strand_type='RNA'):
    """ This block simply creates the reverse complement of a strand of either
    DNA or RNA
    :param sequence: The sequence of DNA/RNA to be reversed
    :type sequence: str
    :param strand_type: The specific strand type that is being used
    :type strand_type: str
    """
    if strand_type == 'RNA':
        sequence = convert_to_RNA(sequence)
        tempseq = sequence.replace('A', 'x')
        tempseq = tempseq.replace('U', 'A')
        tempseq = tempseq.replace('x', 'U')
        tempseq = tempseq.replace('G', 'x')
        tempseq = tempseq.replace('C', 'G')
        tempseq = tempseq.replace('x', 'C')
        sequence = tempseq[::-1]
    if strand_type == 'DNA':
        sequence = convert_to_DNA(sequence)
        tempseq = sequence.replace('A', 'x')
        tempseq = tempseq.replace('T', 'A')
        tempseq = tempseq.replace('x', 'T')
        tempseq = tempseq.replace('G', 'x')
        tempseq = tempseq.replace('C', 'G')
        tempseq = tempseq.replace('x', 'C')
        sequence = tempseq[::-1]
    return sequence
