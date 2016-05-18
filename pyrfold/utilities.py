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

# UNITTESTED CODE

def index_of_sequence(sequence, sub_sequence):
    """this function will look through the parts and identify the start
    and stop location of the part if it exists in the part.
    NOTE: This is 1 indeded
    :param sequence: The full sequence to search
    :type sequence: str
    :param sub_sequence: The sub sequence to search for
    :type sub_sequence: str
    :returns: Start and stop index (1 indexed).
    :rtype: list
    """
    # First convert the string to RNA
    sequence = convert_to_RNA(sequence)
    sub_sequence = convert_to_RNA(sub_sequence)
    try:
        lowindex = sequence.index(sub_sequence) + 1
        highindex = lowindex - 1 + len(sub_sequence)
    except ValueError:
        return 'Target sequence is not in the complete sequence'
    return [lowindex, highindex]


def random_sequence(size, GC_range=None, strand_type='RNA'):
    """
    Simple random nucleotide sequence generator
    Nucleobases accepts RNA or DNA as input

    :type size: int
    :param size: Length of sequence to be returned
    :type GC_range: list of floats
    :param GC_range: Low and high bounds of GC content to be accepted
    :type nucleobase: str
    :param nucleobase: RNA or DNA as input for the type of sequence to return
    """
    # Error checking
    if strand_type not in ['RNA', 'DNA']:
        raise ValueError('Strand_type accepts RNA or DNA')

    if strand_type == 'RNA':
        base_set = ['A', 'U', 'G', 'C']
    elif strand_type == 'DNA':
        base_set = ['A', 'T', 'G', 'C']

    for i in range(100000):
        out = []
        for placecounter in range(size):
            out.append(choice(base_set))
        if GC_range:
            if GC_range[0] > GC_range[1]:
                raise ValueError('GR_range must go from low to high')
            if ((GC_content(''.join(out)) <= GC_range[1]) &
               (GC_content(''.join(out)) >= GC_range[0])):
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


def GC_content(sequence):
    sequence = str(sequence)
    GC_count = sequence.count('G')
    GC_count += sequence.count('C')
    return float(GC_count)/len(sequence)


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


# NOT UNITTESTED CODE
class RNAsequence(object):
    def __init__(self, sequence):
        self.sequence = convert_to_RNA(sequence)

    def reverse_complement(self):
        self.sequence = reverse_complement(self.sequence, strand_type='RNA')

    def __str__(self):
        return self.sequence


def randomly_substitue_u_for_c(sequence, probability=0.5):
    """
    Function for randomly_substituing Us for Cs for a given sequence
    """
    sequence = str(sequence)
    outsequence = []
    for base in sequence:
        if base == 'C':
            if random() <= probability:
                base = 'U'
        outsequence.append(base)
    return ''.join(outsequence)


def convert_to_RNA(sequence):
    """
    Concerts DNA sequence to RNA Sequence
    """
    sequence = sequence.upper()
    return sequence.replace('T', 'U')


def convert_to_DNA(sequence):
    """
    Converts RNA to DNA
    """
    sequence = str(sequence).upper()
    return sequence.replace('U', 'T')


def random_RNA_sequence(size, GC_range=None):
    """
    Simple Random RNA generating random sequence generator
    """
    return random_sequence(size, GC_range, strand_type='RNA')


def random_DNA_sequence(size, GC_range=None):
    """
    Simple Random DNA generating random sequence generator
    """
    return random_sequence(size, GC_range, strand_type='DNA')
